open Bio_crd

let summary crd_path =
  let crd = read_crd crd_path in
  let count = List.length crd.atoms in
  let first =
    match crd.atoms with
    | [] -> "first: n/a"
    | a :: _ ->
        Printf.sprintf "first: %s %s %d (%.3f, %.3f, %.3f)"
          a.resname a.atomname a.serial a.x a.y a.z
  in
  (List.length crd.title), count, first
