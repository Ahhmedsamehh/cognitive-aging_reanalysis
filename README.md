# BIOT 5206 Term Project — Cognitive Aging Microarray Reanalysis

**Author:** Ahmed Sameh
**Dataset:** GSE9990 (Kadish et al. 2009, J. Neurosci.)
**Language:** R (version 4.3 or higher)

## How to reproduce

1. Open `bioinformatics-project.Rmd` in RStudio.
2. Restart R: Session → Restart R.
3. Knit to PDF: click "Knit" (or run `rmarkdown::render("bioinformatics-project.Rmd")`).

The first knit will install required packages automatically (BiocManager, limma, rae230a.db, AnnotationDbi, clusterProfiler, org.Rn.eg.db, pheatmap, RColorBrewer, ggplot2, igraph, magick, dplyr, tidyr, tibble, WGCNA, styler).
Subsequent knits are faster (~2-3 minutes).

## Required tools

- R 4.3+
- LaTeX with xelatex (for PDF output via the `tinytex` package, or a full
  TeX distribution like TinyTeX, MiKTeX, or MacTeX)

## Data source

The expression matrix is loaded directly from the course's public GitHub
repository (no local files needed):
- https://github.com/ahmedmoustafa/gene-expression-datasets/tree/main/datasets/cognitive_aging

## Expected runtime

About 3-5 minutes on a modern laptop, dominated by the 1000-iteration bootstrap in the network bonus section.