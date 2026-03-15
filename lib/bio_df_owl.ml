open Bio_pdb

module DF = Owl_dataframe

let of_atoms atoms =
  let serials = Array.of_list (List.map (fun a -> a.serial) atoms) in
  let names = Array.of_list (List.map (fun a -> a.name) atoms) in
  let res_names = Array.of_list (List.map (fun a -> a.res_name) atoms) in
  let chains =
    Array.of_list (List.map (fun a ->
      match a.chain_id with None -> "" | Some c -> String.make 1 c
    ) atoms)
  in
  let res_seqs = Array.of_list (List.map (fun a -> a.res_seq) atoms) in
  let xs = Array.of_list (List.map (fun a -> a.x) atoms) in
  let ys = Array.of_list (List.map (fun a -> a.y) atoms) in
  let zs = Array.of_list (List.map (fun a -> a.z) atoms) in
  let occs =
    Array.of_list (List.map (fun a -> Option.value ~default:nan a.occupancy) atoms)
  in
  let temps =
    Array.of_list (List.map (fun a -> Option.value ~default:nan a.temp_factor) atoms)
  in
  let elements =
    Array.of_list (List.map (fun a -> Option.value ~default:"" a.element) atoms)
  in
  let data = [|
    DF.pack_int_series serials;
    DF.pack_string_series names;
    DF.pack_string_series res_names;
    DF.pack_string_series chains;
    DF.pack_int_series res_seqs;
    DF.pack_float_series xs;
    DF.pack_float_series ys;
    DF.pack_float_series zs;
    DF.pack_float_series occs;
    DF.pack_float_series temps;
    DF.pack_string_series elements;
  |] in
  let head = [|
    "serial"; "name"; "res_name"; "chain_id"; "res_seq";
    "x"; "y"; "z"; "occupancy"; "temp_factor"; "element"
  |] in
  DF.make ~data head

let to_csv ?sep df filename =
  DF.to_csv ?sep df filename
