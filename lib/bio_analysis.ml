open Base
open Bio_geom

let atoms_to_coords (atoms : Bio_pdb.atom list) =
  List.map atoms ~f:(fun (a : Bio_pdb.atom) ->
    Bio_geom.of_xyz a.Bio_pdb.x a.Bio_pdb.y a.Bio_pdb.z
  )

let center coords =
  match coords with
  | [] -> None
  | _ ->
      let sum = List.fold coords ~init:(Bio_geom.of_xyz 0.0 0.0 0.0) ~f:(fun acc c ->
        {
          Bio_geom.x = acc.Bio_geom.x +. c.Bio_geom.x;
          y = acc.Bio_geom.y +. c.Bio_geom.y;
          z = acc.Bio_geom.z +. c.Bio_geom.z;
        }
      ) in
      let n = Float.of_int (List.length coords) in
      Some (Bio_geom.of_xyz (sum.Bio_geom.x /. n) (sum.Bio_geom.y /. n) (sum.Bio_geom.z /. n))

let rmsd coords_a coords_b =
  let rec loop acc a b =
    match a, b with
    | [], [] -> Some acc
    | ca :: ta, cb :: tb ->
        let dx = ca.Bio_geom.x -. cb.Bio_geom.x in
        let dy = ca.Bio_geom.y -. cb.Bio_geom.y in
        let dz = ca.Bio_geom.z -. cb.Bio_geom.z in
        loop (acc +. (dx *. dx) +. (dy *. dy) +. (dz *. dz)) ta tb
    | _ -> None
  in
  match loop 0.0 coords_a coords_b with
  | None -> None
  | Some sum ->
      let n = Float.of_int (List.length coords_a) in
      if Float.equal n 0.0 then None else Some (Float.sqrt (sum /. n))
