open Bio_pdb

let summary pdb_path =
  let pdb = read_pdb_atoms pdb_path in
  let count = List.length pdb in
  let first =
    match pdb with
    | [] -> "first: n/a"
    | a :: _ ->
        Printf.sprintf "first: %s %s %d (%.3f, %.3f, %.3f)"
          (match a.record with Atom -> "ATOM" | Hetatm -> "HETATM")
          a.name a.serial a.x a.y a.z
  in
  count, first
