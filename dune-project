(lang dune 3.21)

(name BIOCAML)

(generate_opam_files true)

(source
 (github username/reponame))

(authors "Alex Dukhan")

(maintainers "Alex Dukhan")

(license LICENSE)

(documentation https://url/to/documentation)

(package
 (name BIOCAML)
 (synopsis "Structural bioinformatics and molecular dynamics analysis tool")
 (description
  "A CLI tool and library for biophysics, offering an expanding array of standard tools inspired by the bio3d R package.")
 (depends ocaml base)
 (tags
  (bioinformatics molecular-dynamics pdb-parsing structural-biology ocaml)))

; See the complete stanza docs at https://dune.readthedocs.io/en/stable/reference/dune-project/index.html
