let assert_true msg cond =
  if not cond then failwith msg

let count_pdb_records filename =
  In_channel.with_open_text filename @@ fun ic ->
  let rec loop acc =
    match In_channel.input_line ic with
    | None -> acc
    | Some line ->
        let is_atom =
          String.length line >= 6 &&
          (String.sub line 0 6 = "ATOM  " || String.sub line 0 6 = "HETATM")
        in
        loop (if is_atom then acc + 1 else acc)
  in
  loop 0

let read_xyz_count filename =
  In_channel.with_open_text filename @@ fun ic ->
  match In_channel.input_line ic with
  | None -> 0
  | Some line ->
      match int_of_string_opt (String.trim line) with
      | None -> 0
      | Some n -> n

let find_project_root () =
  let rec loop dir =
    let marker = Filename.concat dir "debug_objs" in
    if Sys.file_exists marker then dir
    else
      let parent = Filename.dirname dir in
      if parent = dir then failwith "Could not find project root"
      else loop parent
  in
  loop (Sys.getcwd ())

let () =
  let root = find_project_root () in
  let pdb_path = Filename.concat root "debug_objs/test.pdb" in
  let xyz_path = Filename.concat root "debug_objs/coords.xyz" in
  let dcd_path = Filename.concat root "debug_objs/test.dcd" in

  let raw_count = count_pdb_records pdb_path in
  let atoms = Bio_pdb.read_pdb_atoms pdb_path in
  assert_true "PDB atom count mismatch" (List.length atoms = raw_count);

  let xyz_text = Bio_pdb.atoms_to_xyz atoms in
  let tmp_xyz = "/tmp/biocaml_test.xyz" in
  Out_channel.with_open_text tmp_xyz (fun oc -> output_string oc xyz_text);
  let xyz_coords = Bio_xyz.read_xyz tmp_xyz in
  assert_true "XYZ length mismatch from PDB" (List.length xyz_coords = List.length atoms);

  let crd_atoms = Makecrd.atoms_to_crd_atoms atoms in
  let tmp_crd = "/tmp/biocaml_test.crd" in
  Bio_crd.write_crd tmp_crd crd_atoms;
  let crd = Bio_crd.read_crd tmp_crd in
  assert_true "CRD atom count mismatch" (List.length crd.atoms = List.length atoms);

  let _xyz_n = read_xyz_count xyz_path in
  let xyz_coords2 = Bio_xyz.read_xyz xyz_path in
  assert_true "coords.xyz should parse" (List.length xyz_coords2 > 0);

  let solvate_n = Bio_solvate.estimate_count_from_xyz xyz_path in
  assert_true "solvate result should be positive" (solvate_n > 0);

  let coords_arr = Array.of_list xyz_coords2 in
  let fit_result = Bio_fit.fit coords_arr coords_arr in
  begin match fit_result with
  | None -> failwith "fit failed on identical coords"
  | Some (aligned, _, _, _) ->
      let rmsd = Bio_fit.rmsd coords_arr aligned in
      assert_true "fit RMSD too large" (rmsd < 1e-6)
  end;

  begin match Bio_analysis.rmsd xyz_coords2 xyz_coords2 with
  | None -> failwith "RMSD failed on identical coords"
  | Some r -> assert_true "RMSD too large" (r < 1e-9)
  end;

  let header, frames = Bio_dcd.read_dcd dcd_path in
  assert_true "DCD nset should be > 0" (header.Bio_dcd.nset > 0);
  assert_true "DCD natoms should be > 0" (header.Bio_dcd.natoms > 0);
  assert_true "DCD frame count mismatch" (List.length frames = header.Bio_dcd.nset);
  begin match frames with
  | [] -> failwith "DCD frames empty"
  | f :: _ ->
      assert_true "DCD frame natoms mismatch" (Array.length f = header.Bio_dcd.natoms)
  end;

  let frame1 = Array.of_list xyz_coords2 in
  let frame2 = Array.copy frame1 in
  if Array.length frame2 > 0 then
    frame2.(0) <- { frame2.(0) with Bio_geom.x = frame2.(0).Bio_geom.x +. 0.1 };
  let mean, eigs, scores = Bio_pca.pca ~k:2 [frame1; frame2] in
  assert_true "PCA mean size" (Array.length mean = Array.length frame1 * 3);
  assert_true "PCA eig count" (List.length eigs = 2);
  assert_true "PCA scores rows" (Array.length scores = 2);
  assert_true "PCA scores cols" (Array.length scores.(0) = 2);

  print_endline "All tests passed."
