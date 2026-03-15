type dims =
  | Cube of float
  | Rect of float * float * float

let volume = function
  | Cube a -> a *. a *. a
  | Rect (a, b, c) -> a *. b *. c

let concentration vol = (vol /. 2.) *. 6.022 *. 0.0001 *. 0.150;;

print_endline (string_of_float (concentration ( volume (Rect (566.585, 259.966, 303.434) ) ) ) )