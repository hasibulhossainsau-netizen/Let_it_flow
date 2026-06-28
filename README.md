# Snakemake-CRISPR: Automated End-to-End CRISPR/Cas9 gRNA Design & Evaluation Pipeline

An automated, computational biology pipeline implemented in **Snakemake** to design, validate, and rank highly specific single-guide RNAs (sgRNAs) targeting antibiotic resistance genes, specifically focusing on the beta-lactamase gene **_blaCTX-M-15_**.

This pipeline integrates robust thermodynamic modeling, exact mismatch scanning, and cross-validation to minimize off-target risks across the bacterial genome while maximizing cutting efficiency.

---

## 🔬 Scientific Context & Significance

Antimicrobial Resistance (AMR) is one of the greatest global health crises. The **_blaCTX-M-15_** gene encodes an Extended-Spectrum Beta-Lactamase (ESBL) enzyme, rendering enterobacterial pathogens resistant to crucial antibiotics like penicillins and cephalosporins.

This repository provides an automated solution to design **sequence-specific antimicrobials** via CRISPR-Cas9. By targeting and cleaving this specific gene sequence within the bacterial host, we can selectively knock out resistance or eliminate target plasmids, effectively resensitizing superbugs to conventional treatments.

---

## 🛠️ Pipeline Architecture

The workflow is completely modularized and executed deterministically via Snakemake:

1. **Input Validation:** Audits reference genome and target gene FASTA files for structural integrity.
2. **PAM Scanning:** Scans both forward and reverse strands of the target gene to identify all `NGG` protospacer adjacent motifs.
3. **Biophysical Filtering:** Filters candidates based on strict GC content thresholds ($40\% - 65\%$) and assesses seed region (positions 9–20) composition.
4. **Exact Off-Target Mapping (Cas-OFFinder):** Maps guides across the host genome allowing up to 3 mismatches to capture true spacer-PAM interactions.
5. **Cross-Validation Screen (BLASTn):** Executes a high-sensitivity `blastn-short` screen against the host genome to flag partial sequence homologies that bypass standard PAM-restricted matchers.
6. **Thermodynamic Secondary Structure Folding (RNAfold):** Computes the Minimum Free Energy (MFE) of the spacer conjugated with the truncated sgRNA scaffold to optimize target site accessibility.
7. **Composite Scoring & Ranking:** Ranks candidates using a weighted scoring matrix ($25\%$ Spacer GC, $20\%$ Seed GC, $25\%$ MFE Score, $30\%$ Off-Target Penalty) and assigns definitive selection verdicts.

---

## 📈 Results Summary

When executed against a representative host genome harboring the **_blaCTX-M-15_** sequence, the pipeline yielded the following actionable data compiled into `blaCTX-M-15_gRNA_results_final.xlsx`:

- **Total Screened gRNAs:** 92
- **TOP CANDIDATES (Green List):** 74 — Complete zero off-target signature (`ot_0mm_adjusted == 0` and `blast_offtarget_hits == 0`), optimal GC profile, and highly accessible secondary structures. Directly synthesis-ready.
- **CANDIDATES (Yellow List):** 6 — High specificity but marginal biophysical/folding scores. Serves as robust experimental backups.
- **REJECTED (Red List):** 12 — Eliminated automatically due to critical 0-mismatch off-target mapping or severe multi-copy host genome cross-reactivity.

### Top Performing Candidate Profiles

| Rank | gRNA ID | Spacer Sequence (5' → 3') | PAM | Strand | Position (nt) | Efficiency Score | Selection Verdict |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **1** | `DQ915955_-_828` | `CTTCCTAACAACAGCGTGAC` | GGT | `-` | 828 | **1.0000** | 🟢 TOP CANDIDATE |
| **2** | `DQ915955_-_50` | `CGACGCTAATACATCGCGAC` | GGC | `-` | 50 | **1.0000** | 🟢 TOP CANDIDATE |
| **3** | `DQ915955_+_438` | `AGCTGATTTCTCACGTTGGC` | GGC | `+` | 438 | **1.0000** | 🟢 TOP CANDIDATE |

---

## 🚀 Getting Started

### Prerequisites

Ensure you have Miniconda/Anaconda installed on a Linux/WSL environment. The core tools required are:

- Snakemake (>= 7.0)
- Cas-OFFinder
- NCBI BLAST+ (`blastn`, `makeblastdb`)
- ViennaRNA Suite (`RNAfold`)
- Python 3.11+ with packages: `pandas`, `openpyxl`

### Installation & Execution

1. Clone the repository and navigate into the workspace:

   ```bash
   git clone https://github.com/YOUR_USERNAME/grna_pipeline.git
   cd grna_pipeline
   ```

2. Activate your target Conda environment containing your pipeline dependencies:

   ```bash
   conda activate grna_env
   ```

3. Perform a dry-run to inspect the Directed Acyclic Graph (DAG) of jobs:

   ```bash
   snakemake --cores 4 --dry-run
   ```

4. Execute the fully automated workflow:

   ```bash
   snakemake --cores 4
   ```

---

## 📂 Repository Structure

```
├── config/
│   └── config.yaml             # Global configuration parameters and file paths
├── workflow/
│   └── Snakefile               # Core pipeline definition and dependency routing
├── scripts/
│   ├── 01_validate_inputs.py   # FASTA formatting audit
│   ├── 02_pam_finder.py        # Protospacer & PAM extraction
│   ├── 02b_inspect_candidates.py # Size & primary composition check
│   ├── 03_casoffinder_prep.py  # Generation of Cas-OFFinder inputs
│   ├── 04_parse_offtargets.py  # Processing Cas-OFFinder alignment hits
│   ├── 05_parse_blast.py       # Short-read sequence homology filtering
│   ├── 06_rnafold_score.py     # Secondary structure thermodynamic calculations
│   ├── 07_composite_scorer.py  # Multi-criteria scoring and candidate ranking
│   └── 08_export_excel.py      # Conditional openpyxl formatting engine
├── data/
│   ├── raw/                    # Input reference and gene files
│   └── processed/              # Merged intermediate files and databases
└── results/
    ├── blaCTX-M-15_gRNA_results_final.csv   # Flat table containing final scores
    └── blaCTX-M-15_gRNA_results_final.xlsx  # Color-coded, production-ready final spreadsheet
```

---

## 🎓 Academic Application & Reusability

This pipeline is engineered for scalability. To target any other antibiotic resistance gene (e.g., *NDM-1*, *mcr-1*) or a completely different bacterial host pathogen, simply update the reference paths and parameter constraints inside `config/config.yaml`. Running `snakemake` will automatically regenerate a dedicated, production-ready Excel evaluation report for the new targets.
