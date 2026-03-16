open Base
open Bio_geom
open Bio_xyz

let inside_trunc_oct c side =
  Float.(c.x + c.y + c.z <= side) &&
  Float.(max c.x (max c.y c.z) <= side * 0.5)

let density_to_number density side =
  Int.of_float (density *. (Float.int_pow side 3) *. 0.03342)

let estimate_count_from_xyz ?(density=1.0) filename =
  let coords = read_xyz filename in
  match bounds coords with
  | None -> 0
  | Some (min_v, max_v) ->
      let diameter = sub max_v min_v in
      let side = norm diameter in
      density_to_number density side
