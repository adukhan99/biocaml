open Bio_geom

type header = {
  nset : int;
  istart : int;
  nsavc : int;
  natoms : int;
  endian : [ `Little | `Big ];
}

let read_int32_endian ic endian =
  let b0 = input_byte ic in
  let b1 = input_byte ic in
  let b2 = input_byte ic in
  let b3 = input_byte ic in
  match endian with
  | `Little -> Int32.of_int (b0 lor (b1 lsl 8) lor (b2 lsl 16) lor (b3 lsl 24))
  | `Big -> Int32.of_int (b3 lor (b2 lsl 8) lor (b1 lsl 16) lor (b0 lsl 24))

let read_fortran_record ic endian =
  let len = read_int32_endian ic endian |> Int32.to_int in
  let buf = really_input_string ic len in
  let _ = read_int32_endian ic endian in
  buf

let int32_at_endian buf off endian =
  let b0 = Char.code buf.[off] in
  let b1 = Char.code buf.[off + 1] in
  let b2 = Char.code buf.[off + 2] in
  let b3 = Char.code buf.[off + 3] in
  match endian with
  | `Little -> Int32.of_int (b0 lor (b1 lsl 8) lor (b2 lsl 16) lor (b3 lsl 24))
  | `Big -> Int32.of_int (b3 lor (b2 lsl 8) lor (b1 lsl 16) lor (b0 lsl 24))

let float32_at_endian buf off endian =
  Int32.float_of_bits (int32_at_endian buf off endian)

let detect_endian ic =
  let b0 = input_byte ic in
  let b1 = input_byte ic in
  let b2 = input_byte ic in
  let b3 = input_byte ic in
  let le = Int32.of_int (b0 lor (b1 lsl 8) lor (b2 lsl 16) lor (b3 lsl 24)) |> Int32.to_int in
  let be = Int32.of_int (b3 lor (b2 lsl 8) lor (b1 lsl 16) lor (b0 lsl 24)) |> Int32.to_int in
  let is_plausible n = n >= 32 && n <= 2048 && n mod 4 = 0 in
  let endian =
    if is_plausible le then `Little
    else if is_plausible be then `Big
    else `Little
  in
  endian, (b0, b1, b2, b3)

let read_header ic =
  let endian, (b0, b1, b2, b3) = detect_endian ic in
  let len =
    match endian with
    | `Little -> Int32.of_int (b0 lor (b1 lsl 8) lor (b2 lsl 16) lor (b3 lsl 24)) |> Int32.to_int
    | `Big -> Int32.of_int (b3 lor (b2 lsl 8) lor (b1 lsl 16) lor (b0 lsl 24)) |> Int32.to_int
  in
  let rec1 = really_input_string ic len in
  let _ = read_int32_endian ic endian in
  let nset = int32_at_endian rec1 4 endian |> Int32.to_int in
  let istart = int32_at_endian rec1 8 endian |> Int32.to_int in
  let nsavc = int32_at_endian rec1 12 endian |> Int32.to_int in
  let _title = read_fortran_record ic endian in
  let rec3 = read_fortran_record ic endian in
  let natoms = int32_at_endian rec3 0 endian |> Int32.to_int in
  { nset; istart; nsavc; natoms; endian }

let read_dcd filename =
  In_channel.with_open_bin filename @@ fun ic ->
  let header = read_header ic in
  let natoms = header.natoms in
  let endian = header.endian in
  let read_coord () =
    let buf = read_fortran_record ic endian in
    Array.init natoms (fun i -> float32_at_endian buf (i * 4) endian)
  in
  let rec read_frames acc remaining =
    if remaining = 0 then List.rev acc
    else
      let maybe_uc = read_fortran_record ic endian in
      let xs, ys, zs =
        if String.length maybe_uc = natoms * 4 then
          let xs = Array.init natoms (fun i -> float32_at_endian maybe_uc (i * 4) endian) in
          let ys = read_coord () in
          let zs = read_coord () in
          xs, ys, zs
        else
          let xs = read_coord () in
          let ys = read_coord () in
          let zs = read_coord () in
          xs, ys, zs
      in
      let frame =
        Array.init natoms (fun i -> of_xyz xs.(i) ys.(i) zs.(i))
      in
      read_frames (frame :: acc) (remaining - 1)
  in
  header, read_frames [] header.nset
