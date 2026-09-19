# TDA's file utilities

TDA was written when a data file could be larger than the machine's
memory, and these commands work a record at a time on an ASCII file
rather than loading it: `esort` sorts the records, `eskip` drops
character positions, `eselect` keeps the records whose key appears in a
second file, `emerge` joins sorted files on a key. They are file-in,
file-out, on raw text records addressed by character position, and they
are that here too: name the input and the output file and the output
path comes back, invisibly, with the run as its `"run"` attribute.
Nothing is read back; for a data frame, use
[`tda_read_table`](https://janmarvin.github.io/TDA/reference/tda_read_table.md)
on the file you asked for.

## Usage

``` r
tda_esort(file, keys, out, options = list(), dir = tempfile("tda"))

tda_eskip(file, drop, out, options = list(), dir = tempfile("tda"))

tda_eselect(
  file,
  with,
  keys,
  with_keys = keys,
  out,
  options = list(),
  dir = tempfile("tda")
)

tda_emerge(
  file,
  with,
  keys,
  with_keys = keys,
  out,
  options = list(),
  dir = tempfile("tda")
)

tda_ejoin(
  data,
  id,
  start,
  end,
  state,
  with = NULL,
  with_id = NULL,
  with_start = NULL,
  with_end = NULL,
  with_state = NULL,
  options = list(),
  dir = tempfile("tda")
)
```

## Arguments

- file:

  an existing ASCII data file.

- keys:

  character positions in the record, as a start and end pair – `c(1, 3)`
  is the first three characters. Required: they are positions in the
  line rather than column numbers. Several keys are given as a list of
  pairs, up to five.

- out:

  the output file to write.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- drop:

  character positions to drop, in the same form.

- with:

  a second file supplying the keys, or the file(s) to merge in.

- with_keys:

  where the key sits in `with`, defaulting to the same positions as in
  `file`. `eselect` and `emerge` need a range in each file.

- data:

  for `tda_ejoin`, a data frame.

- id, start, end, state:

  for `tda_ejoin`, the episode columns in `data`: case identifier, spell
  start and end, and the state during the spell.

- with_id, with_start, with_end, with_state:

  for `tda_ejoin`, the same four columns in `with`, when `with` is given
  and its own column names differ from `data`'s; default to `id`/
  `start`/`end`/`state`'s values.

## Value

The output path, invisibly, for the four file utilities; a data frame
for `tda_ejoin`.

## Details

`tda_emerge` is a merge join and requires both inputs sorted in
ascending order on the key; it says so and stops otherwise, naming the
two records that are out of order.

`tda_ejoin` is different in kind: it joins episode data into levels,
takes data frames, and returns one through the export channel.

## See also

Other episodes:
[`tda_episodes()`](https://janmarvin.github.io/TDA/reference/tda_episodes.md)

## Examples

``` r
f <- tempfile(); writeLines(c("103 30", "101 10", "102 20"), f)
o <- tempfile()
tda_esort(f, keys = c(1, 3), out = o)
readLines(o)
#> [1] "101 10" "102 20" "103 30"
tda_eskip(f, drop = c(1, 4), out = o)
readLines(o)
#> [1] "30" "10" "20"

# eselect/emerge need their inputs already sorted on the key
f1 <- tempfile(); writeLines(c("101 10", "102 20", "103 30"), f1)
f2 <- tempfile(); writeLines(c("101 1", "103 3"), f2)
tda_eselect(f1, with = f2, keys = c(1, 3), out = o); readLines(o)
#> [1] "101 10" "103 30"
tda_emerge(f1, with = f2, keys = c(1, 3), out = o); readLines(o)
#> [1] "101 10 101 1 " "102 20       " "103 30 103 3 "

# ejoin: joins one or two episode datasets, splitting overlapping
# spells into levels -- an activity spanning 0-10, and a state that
# changes mid-activity at t=4, joined so each row is a period where
# neither changes
act <- data.frame(id = 1, start = 0, end = 10, state = 1)
st <- data.frame(id = 1, start = c(0, 4), end = c(4, 10), state = c(1, 2))
tda_ejoin(act, id = "id", start = "start", end = "end", state = "state",
         with = st)
#>   id nspell spellno level start end state1 state2
#> 1  1      2       1     0     0   4      1      1
#> 2  1      2       2     0     4  10      1      2
```
