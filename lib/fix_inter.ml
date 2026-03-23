open Base
open Stdio
open Owl
open Cmdliner

module Vec = struct
  type t = Mat.mat
  let of_float x y z = Mat.of_array [|x; y; z|] 3 1
  let x v = Mat.get v 0 0
  let y v = Mat.get v 1 0
  let z v = Mat.get v 2 0
  let add = Mat.add
  let sub = Mat.sub
  let smult = Mat.mul_scalar
  let get_scalar m = 
    match Dense.Ndarray.D.num_dims m with
    | 0 -> Dense.Ndarray.D.get m [||]
    | 1 -> Dense.Ndarray.D.get m [|0|]
    | 2 -> Mat.get m 0 0
    | _ -> Mat.get m 0 0

  let dot v1 v2 = Mat.( (transpose v1 *@ v2) |> fun m -> get_scalar m)
  let norm v = Mat.l2norm v |> get_scalar
  let norm_sq v = let n = Mat.l2norm v |> get_scalar in n *. n
  let dist v1 v2 = Mat.l2norm (Mat.sub v1 v2) |> get_scalar
  let dist_sq v1 v2 = norm_sq (sub v1 v2)
  let normalize v = 
    let n = norm v in
    if Float.(n < 1e-9) then v else Mat.div_scalar v n
end

type atom = {
  id: int;
  name: string;
  resname: string;
  chain: char;
  resseq: int;
  coord: Vec.t;
}

type residue = {
  r_resname: string;
  r_chain: char;
  r_resseq: int;
  r_atoms: atom list;
}

let vdw_radius name =
  let name = String.strip name in
  let elem = 
    if String.length name = 0 then ""
    else if String.length name = 1 then name
    else if Char.is_digit name.[0] then String.sub name ~pos:1 ~len:1
    else String.sub name ~pos:0 ~len:1
  in
  match elem with
  | "H" -> 1.20 | "C" -> 1.70 | "N" -> 1.55 | "O" -> 1.52
  | "P" -> 1.80 | "S" -> 1.80 | "MG" -> 1.73 | "ZN" -> 1.39
  | "FE" -> 1.80 | _ -> 1.50

let is_dna_res r =
  let r = String.strip r in
  List.mem ["DA";"DT";"DG";"DC";"A";"T";"G";"C";"DA5";"DT5";"DG5";"DC5";"DA3";"DT3";"DG3";"DC3"] r ~equal:String.equal

let parse_atom line id =
  let substr i len = String.sub line ~pos:i ~len |> String.strip in
  {
    id;
    name = substr 12 4;
    resname = substr 17 3;
    chain = line.[21];
    resseq = Int.of_string (substr 22 4);
    coord = Vec.of_float 
      (Float.of_string (substr 30 8))
      (Float.of_string (substr 38 8))
      (Float.of_string (substr 46 8));
  }

let group_residues atoms =
  List.group atoms ~break:(fun a1 a2 ->
    not (Char.equal a1.chain a2.chain && Int.equal a1.resseq a2.resseq))
  |> List.map ~f:(fun r_atoms ->
      let a = List.hd_exn r_atoms in
      { r_resname = a.resname; r_chain = a.chain; r_resseq = a.resseq; r_atoms })

let centroid coords =
  let n = List.length coords |> Float.of_int in
  let sum = List.fold coords ~init:(Vec.of_float 0. 0. 0.) ~f:Vec.add in
  Vec.smult sum (1. /. n)

let principal_axis coords =
  let c = centroid coords in
  let x = List.map coords ~f:(fun v -> Vec.sub v c |> Mat.to_array) |> Array.of_list |> Mat.of_arrays in
  let _, _, v = Linalg.D.svd x in
  Mat.get_slice [[0];[]] v |> Mat.transpose

let quicksurf_density p_atom dna_atoms =
  List.fold dna_atoms ~init:0.0 ~f:(fun acc q ->
    let r_sum = vdw_radius p_atom.name +. vdw_radius q.name in
    let d2 = Vec.dist_sq p_atom.coord q.coord in
    acc +. Float.exp (-2.0 *. d2 /. (r_sum *. r_sum))
  )

let density_threshold = 0.2

let check_collision atoms dna_atoms =
  List.exists atoms ~f:(fun p ->
    Float.(quicksurf_density p dna_atoms > density_threshold))

let transform_residue dna_atoms res =
  if is_dna_res res.r_resname then res
  else
    let com = centroid (List.map res.r_atoms ~f:(fun a -> a.coord)) in
    let local_dna = List.filter dna_atoms ~f:(fun a -> Float.(Vec.dist a.coord com < 15.0)) in
    if List.length local_dna < 4 then res
    else
      let dna_coords = List.map local_dna ~f:(fun a -> a.coord) in
      let center = centroid dna_coords in
      let axis = principal_axis dna_coords in
      let op = Vec.sub com center in
      let radial_dir = Vec.normalize (Vec.sub op (Vec.smult axis (Vec.dot op axis))) in
      let rec find_push current_atoms total_shift =
        let max_d = List.fold current_atoms ~init:0.0 ~f:(fun acc p -> Float.max acc (quicksurf_density p local_dna)) in
        if Float.(max_d < density_threshold) then total_shift
        else if Float.(total_shift > 20.0) then total_shift
        else
          let step = 0.5 in
          let next_atoms = List.map current_atoms ~f:(fun a -> { a with coord = Vec.add a.coord (Vec.smult radial_dir step) }) in
          find_push next_atoms (total_shift +. step)
      in
      let shift_mag = find_push res.r_atoms 0.0 in
      if Float.(shift_mag <= 0.0) then res
      else
        let shift_vec = Vec.smult radial_dir shift_mag in
        { res with r_atoms = List.map res.r_atoms ~f:(fun a -> { a with coord = Vec.add a.coord shift_vec }) }

let find_components residues =
  let get_ca_coord res = List.find res.r_atoms ~f:(fun a -> String.equal (String.strip a.name) "CA") |> Option.map ~f:(fun a -> a.coord) in
  let n = List.length residues in
  let res_arr = Array.of_list residues in
  let adj = Array.init n ~f:(fun _ -> []) in
  for i = 0 to n - 1 do
    for j = i + 1 to n - 1 do
      match get_ca_coord res_arr.(i), get_ca_coord res_arr.(j) with
      | Some c1, Some c2 when Float.(Vec.dist c1 c2 <= 5.0) ->
          adj.(i) <- j :: adj.(i); adj.(j) <- i :: adj.(j)
      | _ -> ()
    done
  done;
  let visited = Array.create ~len:n false in
  let rec bfs q acc =
    match q with
    | [] -> acc
    | i :: tl ->
        if visited.(i) then bfs tl acc
        else (visited.(i) <- true; bfs (tl @ adj.(i)) (res_arr.(i) :: acc))
  in
  let comps = ref [] in
  for i = 0 to n - 1 do
    if not visited.(i) then comps := (bfs [i] []) :: !comps
  done;
  !comps

let regroup_protein residues dna_atoms =
  let protein = List.filter residues ~f:(fun r -> not (is_dna_res r.r_resname)) in
  let dna_only_residues = List.filter residues ~f:(fun r -> is_dna_res r.r_resname) in
  let comps = find_components protein in
  match comps with
  | [] | [_] -> residues
  | _ ->
      let sorted = List.sort comps ~compare:(fun c1 c2 -> Int.compare (List.length c2) (List.length c1)) in
      let main = List.hd_exn sorted in
      let others = List.tl_exn sorted in
      let get_ca_coord res = List.find res.r_atoms ~f:(fun a -> String.equal (String.strip a.name) "CA") |> Option.map ~f:(fun a -> a.coord) in
      let main_cas = List.filter_map main ~f:get_ca_coord in
      if List.is_empty main_cas then residues
      else
        let moved_others = List.map others ~f:(fun comp ->
          let comp_cas = List.filter_map comp ~f:get_ca_coord in
          if List.is_empty comp_cas then comp
          else
            let chains_in_comp = List.map comp ~f:(fun r -> r.r_chain) |> List.dedup_and_sort ~compare:Char.compare in
            let chains_in_main = List.map main ~f:(fun r -> r.r_chain) |> List.dedup_and_sort ~compare:Char.compare in
            if not (List.exists chains_in_comp ~f:(fun c -> List.mem chains_in_main c ~equal:Char.equal)) then comp
            else
              let best_dist = ref Float.infinity in
              let best_vec = ref (Vec.of_float 0. 0. 0.) in
              List.iter comp_cas ~f:(fun p1 -> List.iter main_cas ~f:(fun p2 ->
                let d = Vec.dist p1 p2 in if Float.(d < !best_dist) then (best_dist := d; best_vec := Vec.sub p2 p1)));
              let total_shift = Vec.smult !best_vec (1. -. 3.8 /. !best_dist) in
              let steps = Int.of_float (Float.round_up (Vec.norm total_shift /. 0.2)) |> Int.max 1 in
              let step_vec = Vec.smult total_shift (1. /. Float.of_int steps) in
              let rec try_move current_comp step =
                if step > steps then current_comp
                else
                  let next_comp = List.map current_comp ~f:(fun r ->
                    { r with r_atoms = List.map r.r_atoms ~f:(fun a -> { a with coord = Vec.add a.coord step_vec }) })
                  in
                  if check_collision (List.concat_map next_comp ~f:(fun r -> r.r_atoms)) dna_atoms then current_comp
                  else try_move next_comp (step + 1)
              in try_move comp 1
        ) in
        dna_only_residues @ main @ List.concat moved_others

let run input output =
  let lines = In_channel.read_lines input in
  let atom_lines = ref [] in
  List.iteri lines ~f:(fun i l ->
    if String.length l >= 4 && String.equal (String.sub l ~pos:0 ~len:4) "ATOM" then
      atom_lines := (parse_atom l i) :: !atom_lines);
  let atoms = List.rev !atom_lines in
  let residues = group_residues atoms in
  let dna_atoms = List.filter atoms ~f:(fun a -> is_dna_res a.resname) in
  let pushed = List.map residues ~f:(transform_residue dna_atoms) in
  let final_residues = regroup_protein pushed dna_atoms in
  
  let coord_map = Hashtbl.create (module Int) in
  List.iter final_residues ~f:(fun res ->
    List.iter res.r_atoms ~f:(fun a ->
      Hashtbl.set coord_map ~key:a.id ~data:a.coord));
  
  let outc = Out_channel.create output in
  List.iteri lines ~f:(fun i l ->
    match Hashtbl.find coord_map i with
    | Some c ->
        Out_channel.fprintf outc "%s%8.3f%8.3f%8.3f%s\n"
          (String.sub l ~pos:0 ~len:30) (Vec.x c) (Vec.y c) (Vec.z c)
          (if String.length l > 54 then String.sub l ~pos:54 ~len:(String.length l - 54) else "")
    | None -> Out_channel.fprintf outc "%s\n" l);
  Out_channel.close outc

let input = Arg.(required & pos 0 (some string) None & info [] ~docv:"INPUT")
let output = Arg.(required & pos 1 (some string) None & info [] ~docv:"OUTPUT")
let cmd = Cmd.v (Cmd.info "fix_inter") Term.(const run $ input $ output)
let () = Stdlib.exit (Cmd.eval cmd)
