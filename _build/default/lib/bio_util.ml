let safe_sub s start len =
  let n = String.length s in
  if start >= n || len <= 0 then
    ""
  else
    let len' = if start + len > n then n - start else len in
    String.sub s start len'

let trim s = String.trim s

let char_opt s idx =
  if idx < 0 || idx >= String.length s then None else Some s.[idx]

let int_opt s =
  let t = String.trim s in
  if t = "" then None else Some (int_of_string t)

let float_opt s =
  let t = String.trim s in
  if t = "" then None
  else
    try Some (float_of_string t) with
    | Failure _ -> None

let string_opt s =
  let t = String.trim s in
  if t = "" then None else Some t
