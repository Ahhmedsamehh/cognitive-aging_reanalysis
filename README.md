# Cognitive Aging Microarray Reanalysis

![R Version](https://img.shields.io/badge/R-4.3%2B-blue)
![Dataset](https://img.shields.io/badge/Dataset-GSE9990-green)

**Author:** Ahmed Sameh  

**👉 [Click here to view the full rendered analysis and plots](analysis-files/cognitive_aging_re-analysis.md)**

## 📖 Overview
This repository contains an R Markdown workflow developed for the BIOT 5206 course as part of the Biotechnology Master's program at the American University in Cairo. It performs a comprehensive reanalysis of the **GSE9990** microarray dataset originally published in *Kadish et al. 2009, J. Neurosci.* ## 🗂️ Repository Structure
The project files are organized into the following directories:

* **`analysis-files/`** * `cognitive_aging_re-analysis.Rmd`: The primary R Markdown script containing the complete, reproducible analysis pipeline.
  * `cognitive_aging_re-analysis.md`: The GitHub-optimized markdown report (renders code and plots directly in the browser).
  * `cognitive_aging_re-analysis_files/`: Automatically generated folder containing the plots for the markdown report.
  * `cognitive_aging_re-analysis.pdf`: The finalized, static report.
  * `cognitive_aging_re-analysis.html`: An interactive HTML version of the report.
* **`presentation/`**
  * `cognitive-aging.pdf`: The slide deck summarizing the project's methodology, differential expression findings, and network comparisons.

## ⚙️ Prerequisites
To run this pipeline, ensure your system has the following installed:
- **R** (version 4.3 or higher)
- **RStudio**
- A **LaTeX** distribution supporting `xelatex` (e.g., [TinyTeX](https://yihui.org/tinytex/), MiKTeX, or MacTeX) to generate the PDF output.

## 🚀 Getting Started

### 1. Data Source
The workflow is entirely self-contained. The expression matrix is loaded directly from a public GitHub repository, so **no local data downloads are necessary**. 
🔗 [Cognitive Aging Dataset Repository](https://github.com/ahmedmoustafa/gene-expression-datasets/tree/main/datasets/cognitive_aging)

### 2. How to Reproduce
1. Open `analysis-files/cognitive_aging_re-analysis.Rmd` in RStudio.
2. Ensure a clean environment by restarting R (`Session` → `Restart R`).
3. Knit the document by clicking the **"Knit"** button in RStudio.

**Package Installation:** During the first knit, the script will automatically install all required Bioconductor and CRAN packages (including `limma`, `clusterProfiler`, `WGCNA`, `ggplot2`, and `igraph`). 

## ⏱️ Expected Runtime
The pipeline takes about **3–5 minutes** to execute on a modern laptop. The first knit will take longer due to package installations. Subsequent knits are significantly faster, with the runtime primarily dominated by the 1000-iteration bootstrap in the network analysis bonus section.