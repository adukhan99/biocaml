open Bio_pdb

type criteria = {
  chain_id : char option;
  res_name : string option;
  res_seq_min : int option;
  res_seq_max : int option;
  element : string option;
}

let matches (criteria : criteria) (atom : Bio_pdb.atom) =
  let ok_chain =
    match criteria.chain_id, atom.chain_id with
    | None, _ -> true
    | Some c, Some a -> c = a
    | Some _, None -> false
  in
  let ok_res_name =
    match criteria.res_name with
    | None -> true
    | Some rn -> rn = atom.res_name
  in
  let ok_res_min =
    match criteria.res_seq_min with
    | None -> true
    | Some n -> atom.res_seq >= n
  in
  let ok_res_max =
    match criteria.res_seq_max with
    | None -> true
    | Some n -> atom.res_seq <= n
  in
  let ok_element =
    match criteria.element with
    | None -> true
    | Some e -> (match atom.element with Some a -> a = e | None -> false)
  in
  ok_chain && ok_res_name && ok_res_min && ok_res_max && ok_element

let filter_atoms (criteria : criteria) (atoms : Bio_pdb.atom list) =
  List.filter (matches criteria) atoms
