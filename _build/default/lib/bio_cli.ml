open Cmdliner

let cmd_pdb_info =
  let pdb =
    let doc = "Input PDB file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"PDB" ~doc)
  in
  let term =
    Term.(const (fun pdb ->
      let count, first = Pdb_parse.summary pdb in
      Printf.printf "atoms: %d\n%s\n" count first
    ) $ pdb)
  in
  let info = Cmd.info "pdb-info" ~doc:"Summarize a PDB file." in
  Cmd.v info term

let cmd_pdb_to_crd =
  let pdb =
    let doc = "Input PDB file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"PDB" ~doc)
  in
  let out =
    let doc = "Output CRD file." in
    Arg.(required & opt (some string) None & info ["o"; "out"] ~docv:"CRD" ~doc)
  in
  let term =
    Term.(const (fun pdb out ->
      Makecrd.pdb_to_crd pdb out
    ) $ pdb $ out)
  in
  let info = Cmd.info "pdb-to-crd" ~doc:"Convert PDB to CRD." in
  Cmd.v info term

let cmd_pdb_to_xyz =
  let pdb =
    let doc = "Input PDB file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"PDB" ~doc)
  in
  let out =
    let doc = "Output XYZ file." in
    Arg.(required & opt (some string) None & info ["o"; "out"] ~docv:"XYZ" ~doc)
  in
  let term =
    Term.(const (fun pdb out ->
      Makecrd.pdb_to_xyz pdb out
    ) $ pdb $ out)
  in
  let info = Cmd.info "pdb-to-xyz" ~doc:"Convert PDB to XYZ." in
  Cmd.v info term

let cmd_crd_info =
  let crd =
    let doc = "Input CRD file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"CRD" ~doc)
  in
  let term =
    Term.(const (fun crd ->
      let title_lines, count, first = Crd_parse.summary crd in
      Printf.printf "title lines: %d\natoms: %d\n%s\n" title_lines count first
    ) $ crd)
  in
  let info = Cmd.info "crd-info" ~doc:"Summarize a CRD file." in
  Cmd.v info term

let cmd_solvate =
  let xyz =
    let doc = "Input XYZ file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"XYZ" ~doc)
  in
  let density =
    let doc = "Target density (g/cm^3-style scalar). Default: 1.0." in
    Arg.(value & opt float 1.0 & info ["d"; "density"] ~docv:"DENSITY" ~doc)
  in
  let term =
    Term.(const (fun xyz density ->
      let n = Bio_solvate.estimate_count_from_xyz ~density xyz in
      Printf.printf "%d\n" n
    ) $ xyz $ density)
  in
  let info = Cmd.info "solvate" ~doc:"Estimate solvent count from XYZ bounds." in
  Cmd.v info term

let cmd_pdb_bounds =
  let pdb =
    let doc = "Input PDB file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"PDB" ~doc)
  in
  let term =
    Term.(const (fun pdb ->
      let atoms = Bio_pdb.read_pdb_atoms pdb in
      let coords = Bio_analysis.atoms_to_coords atoms in
      match Bio_geom.bounds coords with
      | None -> Printf.printf "n/a\n"
      | Some (min_v, max_v) ->
          let diameter = Bio_geom.sub max_v min_v in
          let side = Bio_geom.norm diameter in
          Printf.printf "min: %f %f %f\n" min_v.x min_v.y min_v.z;
          Printf.printf "max: %f %f %f\n" max_v.x max_v.y max_v.z;
          Printf.printf "diameter: %f\n" side
    ) $ pdb)
  in
  let info = Cmd.info "pdb-bounds" ~doc:"Report bounding box and diameter." in
  Cmd.v info term

let cmd_pdb_center =
  let pdb =
    let doc = "Input PDB file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"PDB" ~doc)
  in
  let term =
    Term.(const (fun pdb ->
      let atoms = Bio_pdb.read_pdb_atoms pdb in
      let coords = Bio_analysis.atoms_to_coords atoms in
      match Bio_analysis.center coords with
      | None -> Printf.printf "n/a\n"
      | Some c -> Printf.printf "%f %f %f\n" c.x c.y c.z
    ) $ pdb)
  in
  let info = Cmd.info "pdb-center" ~doc:"Report geometric center." in
  Cmd.v info term

let cmd_pdb_rmsd =
  let pdb_a =
    let doc = "First PDB file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"PDB_A" ~doc)
  in
  let pdb_b =
    let doc = "Second PDB file." in
    Arg.(required & pos 1 (some string) None & info [] ~docv:"PDB_B" ~doc)
  in
  let term =
    Term.(const (fun pdb_a pdb_b ->
      let a = Bio_pdb.read_pdb_atoms pdb_a |> Bio_analysis.atoms_to_coords in
      let b = Bio_pdb.read_pdb_atoms pdb_b |> Bio_analysis.atoms_to_coords in
      match Bio_analysis.rmsd a b with
      | None -> Printf.eprintf "rmsd: incompatible lengths\n"
      | Some r -> Printf.printf "%f\n" r
    ) $ pdb_a $ pdb_b)
  in
  let info = Cmd.info "pdb-rmsd" ~doc:"RMSD between two PDBs (same atom count)." in
  Cmd.v info term

let cmd_xyz_rmsd =
  let xyz_a =
    let doc = "First XYZ file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"XYZ_A" ~doc)
  in
  let xyz_b =
    let doc = "Second XYZ file." in
    Arg.(required & pos 1 (some string) None & info [] ~docv:"XYZ_B" ~doc)
  in
  let term =
    Term.(const (fun xyz_a xyz_b ->
      let a = Bio_xyz.read_xyz xyz_a in
      let b = Bio_xyz.read_xyz xyz_b in
      match Bio_analysis.rmsd a b with
      | None -> Printf.eprintf "rmsd: incompatible lengths\n"
      | Some r -> Printf.printf "%f\n" r
    ) $ xyz_a $ xyz_b)
  in
  let info = Cmd.info "xyz-rmsd" ~doc:"RMSD between two XYZ coordinate sets." in
  Cmd.v info term

let cmd_pdb_select =
  let pdb =
    let doc = "Input PDB file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"PDB" ~doc)
  in
  let out =
    let doc = "Output PDB file." in
    Arg.(required & opt (some string) None & info ["o"; "out"] ~docv:"PDB_OUT" ~doc)
  in
  let chain =
    let doc = "Chain identifier (single character)." in
    Arg.(value & opt (some string) None & info ["chain"] ~docv:"CHAIN" ~doc)
  in
  let resname =
    let doc = "Residue name (e.g., ALA)." in
    Arg.(value & opt (some string) None & info ["resname"] ~docv:"RESNAME" ~doc)
  in
  let resmin =
    let doc = "Minimum residue sequence number." in
    Arg.(value & opt (some int) None & info ["res-min"] ~docv:"RESMIN" ~doc)
  in
  let resmax =
    let doc = "Maximum residue sequence number." in
    Arg.(value & opt (some int) None & info ["res-max"] ~docv:"RESMAX" ~doc)
  in
  let element =
    let doc = "Element symbol (e.g., C, N, O)." in
    Arg.(value & opt (some string) None & info ["element"] ~docv:"ELEMENT" ~doc)
  in
  let term =
    Term.(const (fun pdb out chain resname resmin resmax element ->
      let atoms = Bio_pdb.read_pdb_atoms pdb in
      let chain_id =
        match chain with
        | None -> None
        | Some s when String.length s > 0 -> Some s.[0]
        | Some _ -> None
      in
      let criteria = Bio_select.{
        chain_id;
        res_name = resname;
        res_seq_min = resmin;
        res_seq_max = resmax;
        element = element;
      } in
      let sel = Bio_select.filter_atoms criteria atoms in
      Bio_pdb.write_pdb out sel
    ) $ pdb $ out $ chain $ resname $ resmin $ resmax $ element)
  in
  let info = Cmd.info "pdb-select" ~doc:"Filter PDB atoms and write a new PDB." in
  Cmd.v info term

let cmd_xyz_fit =
  let xyz =
    let doc = "Input XYZ file (multiple frames supported)." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"XYZ" ~doc)
  in
  let out =
    let doc = "Output XYZ file for aligned frames." in
    Arg.(required & opt (some string) None & info ["o"; "out"] ~docv:"XYZ_OUT" ~doc)
  in
  let ref_xyz =
    let doc = "Optional reference XYZ file (single frame). Defaults to first frame of input." in
    Arg.(value & opt (some string) None & info ["ref"] ~docv:"XYZ_REF" ~doc)
  in
  let term =
    Term.(const (fun xyz out ref_xyz ->
      let frames = Bio_xyz.read_xyz_frames xyz in
      let ref_frame =
        match ref_xyz with
        | Some r -> Bio_xyz.read_xyz r
        | None ->
            match frames with
            | f :: _ -> f
            | [] -> []
      in
      let ref_arr = Array.of_list ref_frame in
      let aligned =
        List.map (fun frame ->
          let mob = Array.of_list frame in
          match Bio_fit.fit ref_arr mob with
          | None -> frame
          | Some (al, _, _, _) -> Array.to_list al
        ) frames
      in
      Bio_xyz.write_xyz_frames out aligned
    ) $ xyz $ out $ ref_xyz)
  in
  let info = Cmd.info "xyz-fit" ~doc:"Align XYZ frames to a reference via Kabsch." in
  Cmd.v info term

let cmd_dcd_info =
  let dcd =
    let doc = "Input DCD file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"DCD" ~doc)
  in
  let term =
    Term.(const (fun dcd ->
      let header, frames = Bio_dcd.read_dcd dcd in
      Printf.printf "frames: %d\n" header.Bio_dcd.nset;
      Printf.printf "natoms: %d\n" header.Bio_dcd.natoms;
      Printf.printf "istart: %d\nnsavc: %d\n" header.Bio_dcd.istart header.Bio_dcd.nsavc;
      Printf.printf "loaded: %d\n" (List.length frames)
    ) $ dcd)
  in
  let info = Cmd.info "dcd-info" ~doc:"Summarize a DCD trajectory." in
  Cmd.v info term

let cmd_dcd_to_xyz =
  let dcd =
    let doc = "Input DCD file." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"DCD" ~doc)
  in
  let out =
    let doc = "Output XYZ file (multi-frame)." in
    Arg.(required & opt (some string) None & info ["o"; "out"] ~docv:"XYZ_OUT" ~doc)
  in
  let term =
    Term.(const (fun dcd out ->
      let _header, frames = Bio_dcd.read_dcd dcd in
      let frames_list = List.map Array.to_list frames in
      Bio_xyz.write_xyz_frames out frames_list
    ) $ dcd $ out)
  in
  let info = Cmd.info "dcd-to-xyz" ~doc:"Convert DCD trajectory to XYZ." in
  Cmd.v info term

let cmd_xyz_pca =
  let xyz =
    let doc = "Input XYZ file (multiple frames)." in
    Arg.(required & pos 0 (some string) None & info [] ~docv:"XYZ" ~doc)
  in
  let k =
    let doc = "Number of principal components." in
    Arg.(value & opt int 3 & info ["k"] ~docv:"K" ~doc)
  in
  let scores =
    let doc = "Output CSV file for scores." in
    Arg.(required & opt (some string) None & info ["scores"] ~docv:"CSV" ~doc)
  in
  let modes =
    let doc = "Optional output XYZ for mode visualization." in
    Arg.(value & opt (some string) None & info ["modes"] ~docv:"XYZ_MODES" ~doc)
  in
  let scale =
    let doc = "Scale factor for mode displacements (if --modes set)." in
    Arg.(value & opt float 1.0 & info ["scale"] ~docv:"SCALE" ~doc)
  in
  let term =
    Term.(const (fun xyz k scores modes scale ->
      let frames = Bio_xyz.read_xyz_frames xyz in
      let frames_arr = List.map Array.of_list frames in
      let mean, eigs, scores_mat = Bio_pca.pca ~k frames_arr in
      let oc = Out_channel.open_text scores in
      for i = 0 to Array.length scores_mat - 1 do
        let row = scores_mat.(i) in
        for j = 0 to Array.length row - 1 do
          if j > 0 then Out_channel.output_char oc ',';
          Out_channel.output_string oc (Printf.sprintf "%f" row.(j))
        done;
        Out_channel.output_char oc '\n'
      done;
      Out_channel.close oc;
      match modes with
      | None -> ()
      | Some out ->
          let natoms = if List.length frames_arr > 0 then Array.length (List.hd frames_arr) else 0 in
          let mean_frame = Bio_pca.unflatten_frame natoms mean in
          let mode_frames =
            List.concat (List.mapi (fun _idx (vec, _lambda) ->
              let pos = Array.map2 (fun m v -> m +. scale *. v) mean vec in
              let neg = Array.map2 (fun m v -> m -. scale *. v) mean vec in
              let fpos = Bio_pca.unflatten_frame natoms pos |> Array.to_list in
              let fneg = Bio_pca.unflatten_frame natoms neg |> Array.to_list in
              [fpos; fneg]
            ) eigs)
          in
          Bio_xyz.write_xyz_frames out (Array.to_list mean_frame :: mode_frames)
    ) $ xyz $ k $ scores $ modes $ scale)
  in
  let info = Cmd.info "xyz-pca" ~doc:"PCA over XYZ frames (scores CSV, optional modes)." in
  Cmd.v info term

let run () =
  let doc = "BIOCAML CLI utilities." in
  let info = Cmd.info "biocaml" ~doc in
  let _default = Cmd.v info Term.(ret (const (`Help (`Pager, None)))) in
  Cmd.eval (Cmd.group info [
    cmd_pdb_info;
    cmd_pdb_to_crd;
    cmd_pdb_to_xyz;
    cmd_crd_info;
    cmd_solvate;
    cmd_pdb_bounds;
    cmd_pdb_center;
    cmd_pdb_rmsd;
    cmd_xyz_rmsd;
    cmd_pdb_select;
    cmd_xyz_fit;
    cmd_dcd_info;
    cmd_dcd_to_xyz;
    cmd_xyz_pca;
  ])
