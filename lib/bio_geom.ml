type coord = {
  x : float;
  y : float;
  z : float;
}

let of_xyz x y z = { x; y; z }

let sub a b = { x = a.x -. b.x; y = a.y -. b.y; z = a.z -. b.z }

let norm v = sqrt ((v.x *. v.x) +. (v.y *. v.y) +. (v.z *. v.z))

let bounds coords =
  match coords with
  | [] -> None
  | first :: rest ->
      let min_v = ref first in
      let max_v = ref first in
      List.iter (fun c ->
        min_v := {
          x = min !min_v.x c.x;
          y = min !min_v.y c.y;
          z = min !min_v.z c.z;
        };
        max_v := {
          x = max !max_v.x c.x;
          y = max !max_v.y c.y;
          z = max !max_v.z c.z;
        };
      ) rest;
      Some (!min_v, !max_v)
