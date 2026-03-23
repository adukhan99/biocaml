open Base
open Stdio

type atom = {
  name: string;
  chain: char;
  resseq: int;
  coord: float * float * float;
}

let parse_atom line =
  let substr i len = String.sub line ~pos:i ~len |> String.strip in
  {
    name = substr 12 4;
    chain = line.[21];
    resseq = Int.of_string (substr 22 4);
    coord = (
      Float.of_string (substr 30 8),
      Float.of_string (substr 38 8),
      Float.of_string (substr 46 8)
    );
  }

let dist (x1,y1,z1) (x2,y2,z2) =
  let dx = x1 -. x2 in
  let dy = y1 -. y2 in
  let dz = z1 -. z2 in
  Float.sqrt (dx *. dx +. dy *. dy +. dz *. dz)

let run input =
  let cas = In_channel.read_lines input
    |> List.filter ~f:(fun l -> String.length l >= 22 && String.equal (String.sub l ~pos:0 ~len:4) "ATOM")
    |> List.map ~f:parse_atom
    |> List.filter ~f:(fun a -> String.equal a.name "CA")
  in
  let chains = List.group cas ~break:(fun a1 a2 -> not (Char.equal a1.chain a2.chain)) in
  List.iter chains ~f:(fun chain ->
    let rec check = function
      | a1 :: a2 :: rest ->
          let d = dist a1.coord a2.coord in
          if Float.(d > 10.0) then
            printf "Chain %c: gap between %d and %d is %.2f A\n" a1.chain a1.resseq a2.resseq d;
          check (a2 :: rest)
      | _ -> ()
    in
    check chain
  )

let () =
  if Array.length (Stdlib.Sys.argv) < 2 then printf "Usage: check_bonds input.pdb\n"
  else run (Stdlib.Sys.argv.(1))
