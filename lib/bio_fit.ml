open Base
open Bio_geom

let centroid coords =
  let n = Array.length coords in
  if n = 0 then of_xyz 0.0 0.0 0.0
  else
    let sum =
      Array.fold coords ~init:(of_xyz 0.0 0.0 0.0) ~f:(fun acc c ->
        {
          x = acc.x +. c.x;
          y = acc.y +. c.y;
          z = acc.z +. c.z;
        }
      )
    in
    let f_n = Float.of_int n in
    of_xyz (sum.x /. f_n) (sum.y /. f_n) (sum.z /. f_n)

let center coords =
  let c = centroid coords in
  Array.map coords ~f:(fun p -> { x = p.x -. c.x; y = p.y -. c.y; z = p.z -. c.z }), c

let cov_matrix a b =
  let sxx = ref 0.0 and sxy = ref 0.0 and sxz = ref 0.0 in
  let syx = ref 0.0 and syy = ref 0.0 and syz = ref 0.0 in
  let szx = ref 0.0 and szy = ref 0.0 and szz = ref 0.0 in
  for i = 0 to Array.length a - 1 do
    let ax = a.(i).x and ay = a.(i).y and az = a.(i).z in
    let bx = b.(i).x and by = b.(i).y and bz = b.(i).z in
    sxx := !sxx +. ax *. bx;
    sxy := !sxy +. ax *. by;
    sxz := !sxz +. ax *. bz;
    syx := !syx +. ay *. bx;
    syy := !syy +. ay *. by;
    syz := !syz +. ay *. bz;
    szx := !szx +. az *. bx;
    szy := !szy +. az *. by;
    szz := !szz +. az *. bz;
  done;
  (!sxx, !sxy, !sxz, !syx, !syy, !syz, !szx, !szy, !szz)

let quat_from_cov (sxx, sxy, sxz, syx, syy, syz, szx, szy, szz) =
  let trace = sxx +. syy +. szz in
  let k00 = trace in
  let k01 = syz -. szy in
  let k02 = szx -. sxz in
  let k03 = sxy -. syx in
  let k11 = sxx -. syy -. szz in
  let k12 = sxy +. syx in
  let k13 = szx +. sxz in
  let k22 = -.sxx +. syy -. szz in
  let k23 = syz +. szy in
  let k33 = -.sxx -. syy +. szz in
  let mat = [|
    [| k00; k01; k02; k03 |];
    [| k01; k11; k12; k13 |];
    [| k02; k12; k22; k23 |];
    [| k03; k13; k23; k33 |];
  |] in
  let v = ref [| 1.0; 0.0; 0.0; 0.0 |] in
  for _ = 1 to 60 do
    let w =
      Array.init 4 ~f:(fun i ->
        mat.(i).(0) *. (!v).(0)
        +. mat.(i).(1) *. (!v).(1)
        +. mat.(i).(2) *. (!v).(2)
        +. mat.(i).(3) *. (!v).(3))
    in
    let norm =
      Float.sqrt (Array.fold w ~init:0.0 ~f:(fun acc x -> acc +. x *. x))
    in
    if Float.(norm > 0.0) then
      v := Array.map w ~f:(fun x -> x /. norm)
  done;
  (!v).(0), (!v).(1), (!v).(2), (!v).(3)

let rot_from_quat (q0, q1, q2, q3) =
  let q00 = q0 *. q0 and q11 = q1 *. q1 and q22 = q2 *. q2 and q33 = q3 *. q3 in
  let q01 = q0 *. q1 and q02 = q0 *. q2 and q03 = q0 *. q3 in
  let q12 = q1 *. q2 and q13 = q1 *. q3 and q23 = q2 *. q3 in
  [|
    [| q00 +. q11 -. q22 -. q33; 2.0 *. (q12 -. q03); 2.0 *. (q13 +. q02) |];
    [| 2.0 *. (q12 +. q03); q00 -. q11 +. q22 -. q33; 2.0 *. (q23 -. q01) |];
    [| 2.0 *. (q13 -. q02); 2.0 *. (q23 +. q01); q00 -. q11 -. q22 +. q33 |];
  |]

let apply_rot rot p =
  {
    x = rot.(0).(0) *. p.x +. rot.(0).(1) *. p.y +. rot.(0).(2) *. p.z;
    y = rot.(1).(0) *. p.x +. rot.(1).(1) *. p.y +. rot.(1).(2) *. p.z;
    z = rot.(2).(0) *. p.x +. rot.(2).(1) *. p.y +. rot.(2).(2) *. p.z;
  }

let fit ref_coords mob_coords =
  if Array.length ref_coords <> Array.length mob_coords then
    None
  else
    let ref_centered, ref_centroid = center ref_coords in
    let mob_centered, mob_centroid = center mob_coords in
    let cov = cov_matrix mob_centered ref_centered in
    let q = quat_from_cov cov in
    let rot = rot_from_quat q in
    let aligned =
      Array.map mob_centered ~f:(fun p ->
        let pr = apply_rot rot p in
        { x = pr.x +. ref_centroid.x; y = pr.y +. ref_centroid.y; z = pr.z +. ref_centroid.z }
      )
    in
    Some (aligned, rot, ref_centroid, mob_centroid)

let rmsd ref_coords aligned =
  let n = Array.length ref_coords in
  if n = 0 then 0.0
  else
    let sum = ref 0.0 in
    for i = 0 to n - 1 do
      let dx = ref_coords.(i).x -. aligned.(i).x in
      let dy = ref_coords.(i).y -. aligned.(i).y in
      let dz = ref_coords.(i).z -. aligned.(i).z in
      sum := !sum +. dx *. dx +. dy *. dy +. dz *. dz
    done;
    Float.sqrt (!sum /. Float.of_int n)
