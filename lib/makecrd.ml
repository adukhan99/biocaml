open Bio_pdb
open Bio_crd

let atoms_to_crd_atoms (atoms : Bio_pdb.atom list) =
  List.map (fun (a : Bio_pdb.atom) ->
    {
      Bio_crd.serial = a.serial;
      resid = a.Bio_pdb.res_seq;
      resname = a.Bio_pdb.res_name;
      atomname = a.Bio_pdb.name;
      x = a.Bio_pdb.x;
      y = a.Bio_pdb.y;
      z = a.Bio_pdb.z;
      segid = Option.map (String.make 1) a.Bio_pdb.chain_id;
      resid2 = None;
      weight = None;
    }
  ) atoms

let pdb_to_crd pdb_path crd_path =
  let pdb = read_pdb_atoms pdb_path in
  let crd_atoms = atoms_to_crd_atoms pdb in
  write_crd crd_path crd_atoms

let pdb_to_xyz pdb_path xyz_path =
  let pdb = read_pdb_atoms pdb_path in
  let xyz = atoms_to_xyz pdb in
  Out_channel.with_open_text xyz_path (fun oc -> output_string oc xyz)
