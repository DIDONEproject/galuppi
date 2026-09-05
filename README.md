# galuppi

Repository reproducing the analysis from the paper:

> Llorens, A., García-Portugués, E., and Torrente, Á. (2026). By Baldassare Galuppi? Inferring the authorship of three untraceable *opera seria* arias. *Submitted*.

## Data

Three corpora of eighteenth-century arias, curated at the ERC Didone project and described by their `musif` features, back the three classification problems of the paper:

| Corpus | Arias | Problem |
| --- | --- | --- |
| `code/gp.csv` | 300 | Binary, Baldassare Galuppi vs. Davide Perez; the 3 arias of unknown authorship carry `NA` as `Composer` |
| `code/gp_multi.csv` | 426 | *Demofoonte*, 5 candidate composers |
| `code/gp_ale.csv` | 655 | *Alessandro nell'Indie*, 6 candidate composers |

Each file has an `AriaId` column, a `Composer` column, and the 8,244 `musif` features extracted from the scores. `code/aria_names.csv` maps every `AriaId` to its title and opera. The documentation of the features and of the software that extracts them is available at <https://musif.didone.eu>.

## Reproducing the analysis

From the repository root:

```sh
quarto render code/galuppi.qmd
```

The `knitr` cache of the thirty modelling chunks is shipped in `code/galuppi_cache/`, so this is a warm render: it restores the fitted models instead of refitting them, and reproduces the numbers of the paper exactly.

The render writes:

| Output | Contents |
| --- | --- |
| `code/galuppi.html` | Full analysis report |
| `paper_numbers.tex` | The `\gp...` macros the manuscript inputs |
| `code/paper_numbers.csv` | The same numbers as a table |
| `imgs/` | The three figures of the paper |
| `code/discussion/exemplar_arias.csv`, `code/discussion/ale_exemplar_arias.csv`, `code/discussion/multi_exemplar_arias.csv` | Arias exemplifying each selected stylistic marker |
| `code/marker_correlates.csv`, `code/ale_marker_correlates.csv`, `code/multi_marker_correlates.csv` | Features correlating with each selected marker |
| `code/discussion/marginal_markers.csv`, `code/discussion/doubtful_arias_markers.csv` | Marker summaries for the corpora and for the three doubtful arias |

The folder `code/discussion/` gathers the evidence behind the Discussion section of the paper, namely the exemplar arias of each stylistic marker, the marker summaries, and the two counts taken from the scores themselves along with the scripts that produce them.

`code/coverage_study.R` is the simulation study of the appendix, which measures the coverage of the intervals at the settings of the three problems.

### Requirements

R (>= 4.4), [Quarto](https://quarto.org), a Java runtime for the h2o cluster, and the R packages `tidyverse`, `tidymodels`, `agua`, `h2o`, `glmnet`, `MASS`, `ks`, `doParallel`, `foreach`, `future`, `doFuture`, `progressr`, `withr` and `testthat`; `code/discussion/viola_doubling.py` additionally needs Python with `music21`.