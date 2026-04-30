# KRAS sequential treatment

Multi-omic characterisation of the response to KRAS-G12C inhibitors —
**sotorasib (AMG510)** and **adagrasib (MRTX849)** — across two human
KRAS-G12C cell lines (**H23**, **CALU-1**) and  KRasloxKRASMUT lines
(**MEM**, **MEP**) under three conditions per line: NAIVE, sotorasib,
adagrasib.


## Script directory

Scripts live under `Script/` grouped by analysis topic. Every script is a
**stand-alone Rscript**: drop into the project root and run e.g.
`Rscript Script/cnv/cnv_jaccard.R`. They share helpers via two
`source()`-able utility files in `Script/utils/`. Within a topic group
the order is "build before plot" where applicable.

```
Script/
├── utils/
│   ├── utils_RNA.R                    ← get_heatmap3() and small RNA helpers
│   └── utils_ORA.R                    ← ORA() + add_DOSE_measure_to_ORA()
│
├── circos/                            ← per-line circos figures
│   ├── Calu_circos.R                  ← Calu-1 (human) — Figures/circos/CALU_circos.pdf
│   ├── H23_circos.R                   ← H23   (human) — Figures/circos/H23_circos.pdf
│   ├── MEP_circos.R                   ← MEP   (mouse) — Figures/circos/MEP_circos.pdf
│   └── sv_circos_export.R             ← export the SV slice each circos uses to xlsx
│
├── cnv/                               ← copy-number visualisations
│   ├── cnv_upset_heatmaps.R           ← UpSet + ComplexHeatmaps across the 3 lines
│   │                                   (H23 + CALU share a global CNt scale; MEP keeps its own)
│   ├── cnv_jaccard.R                  ← triangular Jaccard heatmaps for Amp / Loss / Amp ∪ Loss
│   └── cnv_drug_only_signature_heatmap.R
│                                       ← signature-restricted "drug-only" CNV heatmaps
│                                       (Amp/Loss in sotorasib or adagrasib, absent in naive)
│
├── mutations/                         ← SNV / small-variant figures
│   ├── mutations_upset.R              ← SNV UpSet across H23/CALU/MEP
│   ├── Calu_mutations_heatmap.R
│   ├── H23_mutations_heatmap.R
│   └── MEP_mutations_heatmap.R
│
├── ora/                               ← over-representation analysis
│   ├── run_ORA_H23_CALU.R             ← runs ORA on the 4 human contrasts → ora_H23.rds, ora_CALU.rds
│   ├── plot_ORA_H23_CALU.R            ← consumes the rds → dot-plots in Figures/ORA/
│   ├── run_ORA_MEM_MEP.R              ← mouse counterpart (msigdbr at runtime, no go_list)
│   └── plot_ORA_MEM_MEP.R             ← Figures/ORA/MEM_MEP/
│
└── rna/                               ← RNA-seq signature & DEG plots
    ├── heatmaps_MAPK_ERK.R       ← per-line signature heatmaps (H23, CALU)
    ├── heatmaps_MAPK_ERK_MEM_MEP.R
    ├── summary_MAPK_ERK.R             ← A4-friendly summary panels (H23, CALU)
    ├── summary_MAPK_ERK_MEM_MEP.R     ← mouse counterpart
    └── alluvial_DEG_signatures.R      ← alluvial flow + matching UpSet across cell-line / drug pairs
```

Run order, end-to-end: `utils/` is sourced as needed; `ora/run_*` produces
the rds files consumed by `ora/plot_*`; everything else is independent and
can run in any order. Typical full-rerun time on a laptop: ~15 minutes
(ORA + clustering on the larger CNV heatmaps dominates).

---

## Conventions

A handful of project-wide conventions, called out so future contributors
don't have to rederive them.

### Sample → treatment mapping (RNA counts)

The `H23_CALU_counts.rds` and `MEM_MEP_counts.rds` files have **no
metadata sidecar** — sample columns are encoded as
`caks-<line>-<lineCode>-<batch>-<XYZ>` and the **second digit (Y) of the
trailing condition code is the treatment**:

| Y | treatment |
|---|---|
| `0` | NAIVE |
| `1` | sotorasib (AMG510) |
| `2` | adagrasib (MRTX) |

Verified empirically (Spearman ≈ +1) by comparing inferred log2FC against
the deg `log2FoldChange` for both contrasts in all four cell lines.
**Don't re-derive** — every sample-handling script applies the rule via
`samples_for_line()` helpers.

### Drug labels

Display as **`sotorasib`** / **`adagrasib`** in plots; `AMG510` / `MRTX`
stay only in raw deg `contrast` strings and intermediate file paths. The
project-wide `relabel_contrast()` helpers strip `_vs_NAIVE` and apply
the substitution.

### DEG-input filter

Project-wide cutoffs for any signature-restricted or summary plot:

```
padj < 0.01  AND  |log2FoldChange| > 0.5
```

Tighter than originally explored; settled after iteration. NA padj rows
are dropped before applying.

### ORA plot-side filters

Two ontology groups, two filters:

| Group | Members | Filter |
|---|---|---|
| `kh` | kegg + hallmarks | `p.adjust ≤ 0.1` |
| `go` | go-bp | `p.adjust < 0.001` AND `FoldEnrichment ≥ 2` |

go-bp at 0.1 returned ~600 hits per panel; the tighter cutoff brings it
down to readable territory while keeping kegg / hallmarks browseable.


---

## Running

R ≥ 4.4 with the package set used by the various scripts:

```r
# Core
install.packages(c("dplyr","tidyr","tibble","stringr","forcats","ggplot2",
                   "patchwork","ggh4x","scico","ggalluvial","writexl",
                   "RColorBrewer","circlize","msigdbr","classInt"))
# Bioconductor
BiocManager::install(c("ComplexHeatmap","clusterProfiler","DOSE"))
```

Then from the project root:

```bash
# Genomic events (independent of RNA)
Rscript Script/circos/Calu_circos.R
Rscript Script/circos/H23_circos.R
Rscript Script/circos/MEP_circos.R
Rscript Script/circos/sv_circos_export.R          # → SV_circos_selected.xlsx

Rscript Script/cnv/cnv_upset_heatmaps.R
Rscript Script/cnv/cnv_jaccard.R
Rscript Script/cnv/cnv_drug_only_signature_heatmap.R

Rscript Script/mutations/mutations_upset.R
Rscript Script/mutations/Calu_mutations_heatmap.R
Rscript Script/mutations/H23_mutations_heatmap.R
Rscript Script/mutations/MEP_mutations_heatmap.R

# RNA pipeline (ORA must run before the ORA plotters)
Rscript Script/ora/run_ORA_H23_CALU.R              # → ora_H23.rds, ora_CALU.rds
Rscript Script/ora/plot_ORA_H23_CALU.R             # → Figures/ORA/
Rscript Script/ora/run_ORA_MEM_MEP.R               # → ora_MEM.rds, ora_MEP.rds
Rscript Script/ora/plot_ORA_MEM_MEP.R              # → Figures/ORA/MEM_MEP/

Rscript Script/rna/heatmap_RAS84_MAPK_ERK.R
Rscript Script/rna/heatmap_RAS84_MAPK_ERK_MEM_MEP.R
Rscript Script/rna/summary_MAPK_ERK.R
Rscript Script/rna/summary_MAPK_ERK_MEM_MEP.R
Rscript Script/rna/volcano_pathways.R
Rscript Script/rna/alluvial_DEG_signatures.R
```

Each script resolves its own location (`commandArgs(--file=)`) and walks
**two levels up** to find the project root — so it works regardless of
the caller's working directory.

---

