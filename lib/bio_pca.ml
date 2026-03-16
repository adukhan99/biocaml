open Base
open Bio_geom

let flatten_frame frame =
  let n = Array.length frame in
  let v = Array.create ~len:(n * 3) 0.0 in
  for i = 0 to n - 1 do
    let j = i * 3 in
    v.(j) <- frame.(i).x;
    v.(j + 1) <- frame.(i).y;
    v.(j + 2) <- frame.(i).z;
  done;
  v

let unflatten_frame natoms v =
  Array.init natoms ~f:(fun i ->
    let j = i * 3 in
    of_xyz v.(j) v.(j + 1) v.(j + 2)
  )

let mean_vector samples =
  let m = Array.length samples in
  if m = 0 then [||]
  else
    let d = Array.length samples.(0) in
    let mean = Array.create ~len:d 0.0 in
    Array.iter samples ~f:(fun v ->
      for i = 0 to d - 1 do
        mean.(i) <- mean.(i) +. v.(i)
      done
    );
    for i = 0 to d - 1 do
      mean.(i) <- mean.(i) /. Float.of_int m
    done;
    mean

let cov_matvec samples mean v =
  let m = Array.length samples in
  let d = Array.length v in
  let out = Array.create ~len:d 0.0 in
  if m <= 1 then out
  else
    let scale = 1.0 /. Float.of_int (m - 1) in
    Array.iter samples ~f:(fun x ->
      let dot_v = ref 0.0 in
      for i = 0 to d - 1 do
        dot_v := !dot_v +. (x.(i) -. mean.(i)) *. v.(i)
      done;
      for i = 0 to d - 1 do
        out.(i) <- out.(i) +. (x.(i) -. mean.(i)) *. !dot_v
      done
    );
    for i = 0 to d - 1 do
      out.(i) <- out.(i) *. scale
    done;
    out

let normalize v =
  let norm = Float.sqrt (Array.fold v ~init:0.0 ~f:(fun acc x -> acc +. x *. x)) in
  if Float.equal norm 0.0 then v else Array.map v ~f:(fun x -> x /. norm)

let dot a b =
  let acc = ref 0.0 in
  for i = 0 to Array.length a - 1 do
    acc := !acc +. a.(i) *. b.(i)
  done;
  !acc

let orthonormalize v basis =
  let v' = Array.copy v in
  List.iter basis ~f:(fun b ->
    let proj = dot v' b in
    for i = 0 to Array.length v' - 1 do
      v'.(i) <- v'.(i) -. proj *. b.(i)
    done
  );
  normalize v'

let power_iteration ?(iters=80) samples mean basis =
  let d = Array.length mean in
  let v0 = Array.init d ~f:(fun i -> if i % 2 = 0 then 1.0 else -1.0) in
  let v = ref (orthonormalize v0 basis) in
  for _ = 1 to iters do
    let w = cov_matvec samples mean !v in
    v := orthonormalize w basis
  done;
  let lambda = dot !v (cov_matvec samples mean !v) in
  !v, lambda

let pca ?(k=3) frames =
  let m = List.length frames in
  if m = 0 then ([||], [], [||])
  else
    let samples =
      Array.of_list (List.map frames ~f:flatten_frame)
    in
    let mean = mean_vector samples in
    let rec loop i basis eigs =
      if i = k then List.rev eigs, basis
      else
        let vec, lambda = power_iteration samples mean basis in
        loop (i + 1) (basis @ [vec]) ((vec, lambda) :: eigs)
    in
    let eigs, basis = loop 0 [] [] in
    let scores =
      Array.map samples ~f:(fun x ->
        Array.of_list (List.map basis ~f:(fun v ->
          dot (Array.map2_exn x mean ~f:(fun xi mi -> xi -. mi)) v
        ))
      )
    in
    mean, eigs, scores
