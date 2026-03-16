open Base

let safe_sub s start len =
  let n = String.length s in
  if start >= n || len <= 0 then
    ""
  else
    let len' = if start + len > n then n - start else len in
    String.sub s ~pos:start ~len:len'

let char_opt s idx =
  if idx < 0 || idx >= String.length s then None else Some s.[idx]

let int_opt s =
  Int.of_string_opt (String.strip s)

let float_opt s =
  Float.of_string_opt (String.strip s)

let string_opt s =
  let t = String.strip s in
  if String.is_empty t then None else Some t
