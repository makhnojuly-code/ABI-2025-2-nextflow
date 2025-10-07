# ABI-2025-2 Nextflow Pipelines

This repository contains training workflows developed during the **ABI 2025 Nextflow course**.
## 🧬 Pipelines
### 1️⃣ `sra_download.nf`
Downloads FASTQ data from the SRA database, processes it, and prepares it for downstream analysis.
### 2️⃣ `sap.nf`
Performs genome assembly using **SPAdes** and evaluates assembly quality with **QUAST**.
## ⚙️ Configuration
The `nextflow.config` file contains default parameters and directory settings.
## ▶️ How to run
Example:
```bash
nextflow run sra_download.nf -profile singularity
nextflow run sap.nf -profile sap
