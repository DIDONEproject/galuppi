# Counts the inversions of the tonic chords of every aria of gp.csv, from the
# Roman-numeral annotations of the MuseScore sources, and writes
# tonic_inversions.csv. The extraction of galuppi.qmd carries chord roots only,
# so the inversions are read here from the .mscx annotations themselves.
# The output is skipped when already present, so delete the CSV to rebuild it.
# Requires the MuseScore corpus, whose path is set in scores_dir below.
# Run from code/discussion/, as the paths are relative to it.

# Required libraries
library(tidyverse)

# Annotation parsing
{

# The MuseScore sources, outside the repository. The tracked
# tonic_inversions.csv makes galuppi.qmd independent of them.
scores_dir <- "~/Documents/Corpus_DIDONE/musescore"

# The Roman-numeral labels of one .mscx file, in order of appearance
read_labels <- function(file) {

  xml <- paste(readLines(con = file, warn = FALSE), collapse = "\n")
  matches <- str_match_all(xml, "(?s)<Harmony>.*?<name>(.*?)</name>")[[1]]
  matches[, 2]

}

# Numeral and inversion figure of one label, or NA for the labels that name no
# chord of the prevailing key, altered qualities (o, %, +) counting as chords
# of their own. Phrase and pedal markers are dropped, a local
# key prefix is resolved (G.I to I), parenthesised figures are added notes
# rather than inversions (I(64) is a root-position tonic), and applied chords
# belong to a secondary key.
parse_label <- function(label) {

  x <- str_remove_all(label, "[{}\\[\\]|]")
  x <- str_remove_all(x, "\\(.*?\\)")
  if (str_detect(x, "/") || x == "") return(c(NA, NA))
  x <- str_split(x, fixed("."))[[1]] |> tail(n = 1)
  numeral <- str_match(x, paste0("^(b|#)?",
                                 "(VII|VI|IV|V|III|II|I|vii|vi|iv|v|iii|ii|i)",
                                 "(o|%|\\+)?"))
  if (is.na(numeral[1])) return(c(NA, NA))
  figure <- str_remove(x, fixed(numeral[1])) |>
    str_remove_all("[^0-9]")
  c(paste0(numeral[2] %|% "", numeral[3], numeral[4] %|% ""), figure)

}

# str_match returns NA for an absent optional group
`%|%` <- function(x, y) if (is.na(x)) y else x

# The inversion figures, by the chord member they put in the bass
root_figures <- c("", "7", "9")
first_figures <- c("6", "65")
second_figures <- c("64", "43")

# Tonic-chord counts of one aria, over the labels of its score
count_tonics <- function(file) {

  parsed <- map(read_labels(file), parse_label)
  numerals <- map_chr(parsed, 1)
  figures <- map_chr(parsed, 2)
  keep <- !is.na(numerals)
  numerals <- numerals[keep]
  figures <- figures[keep]
  tonic <- numerals %in% c("I", "i")
  tibble(n_labels = length(numerals),
         n_tonic = sum(tonic),
         tonic_root = sum(tonic & figures %in% root_figures),
         tonic_first = sum(tonic & figures %in% first_figures),
         tonic_second = sum(tonic & figures %in% second_figures))

}

}

# Corpus pass
{

if (!file.exists("tonic_inversions.csv")) {

# One score per AriaId; the bracketed number of the filename is the AriaId
files <- list.files(path = path.expand(scores_dir), pattern = "\\.mscx$",
                    full.names = TRUE)
ids <- str_match(basename(files), "\\[(\\d+)\\]\\.mscx$")[, 2]
scores <- tibble(AriaId = ids, file = files) |>
  filter(!is.na(AriaId)) |>
  distinct(AriaId, .keep_all = TRUE)

# The arias of the binary corpus, the three doubtful ones included. AriaId is
# read as text, since the identifiers of the filenames are zero-padded.
gp <- read.csv(file = "../gp.csv", header = TRUE,
               colClasses = c("AriaId" = "character")) |>
  select(AriaId, Composer)
missing <- setdiff(gp$AriaId, scores$AriaId)
if (length(missing) > 0) {
  stop("no MuseScore source for AriaId ", paste(missing, collapse = ", "))
}

counts <- gp |>
  left_join(scores, by = "AriaId") |>
  rowwise() |>
  mutate(count_tonics(file)) |>
  ungroup() |>
  select(-file)
write.csv(counts, file = "tonic_inversions.csv", row.names = FALSE)

# Shares by composer, as reported in the paper
counts |>
  filter(Composer != "NA") |>
  group_by(Composer) |>
  summarise(arias = n(), tonics = sum(n_tonic),
            first = sum(tonic_first) / sum(n_tonic),
            second = sum(tonic_second) / sum(n_tonic)) |>
  print()

}

}
