# Variant Annotation Pipeline for Rare Genetic Diseases

A reproducible **GRCh38-based rare-disease annotation workflow** for:

- small sequence variants: SNVs and small insertions/deletions in VCF format;
- copy-number variants: deletions and duplications in BED format.

The project was developed and tested using synthetic educational datasets for:

1. Tay-Sachs disease
2. Kabuki syndrome
3. Sotos syndrome
4. Noonan syndrome

> **Important:** All variants in this repository are synthetic educational data. They do not represent real patients and must not be used for clinical diagnosis or patient management.

---

## Project objective

The aim of this project is to convert raw, unannotated genomic variants into organized and biologically interpretable results using a reproducible Linux-based workflow.

The pipeline:

- validates and normalizes small variants against GRCh38;
- annotates genes, transcripts, consequences and HGVS descriptions;
- adds ClinVar, ClinGen, gnomAD and SpliceAI evidence;
- annotates and classifies CNVs using multiple independent tools;
- creates disease-specific result folders;
- stores logs, intermediate files and final outputs;
- supports analysis of one disease or all four diseases.

---

## Diseases and synthetic targets

| Disease | Main gene | Inheritance | Synthetic small-variant target | Synthetic CNV target |
|---|---|---|---|---|
| Tay-Sachs disease | `HEXA` | Autosomal recessive | `rs387906309`, `NM_000520.6:c.1274_1277dup`, `p.Tyr427IlefsTer5` | `chr15:72340923-72376014 DEL` |
| Kabuki syndrome | `KMT2D` | Autosomal dominant | `rs398123721`, `NM_003482.4:c.14710C>T`, `p.Arg4904Ter` | `chr12:49018977-49060794 DEL` |
| Sotos syndrome | `NSD1` | Autosomal dominant | `rs587784148`, `NM_022455.5:c.5431C>T`, `p.Arg1811Ter` | `chr5:177131797-177300213 DEL` |
| Noonan syndrome | `PTPN11` | Autosomal dominant | `rs121918459`, `NM_002834.5:c.188A>G`, `p.Tyr63Cys` | `chr12:112418946-112509918 DUP` |

The sequence-variant and CNV files are separate synthetic analytical scenarios. They should not automatically be interpreted as two findings from the same patient.

---

## Workflow overview

### Small-variant branch

```text
Unannotated VCF
    |
    v
bgzip + tabix
    |
    v
bcftools norm
    |
    v
VEP 115
    |
    v
SnpEff
    |
    v
ClinVar annotation
    |
    v
ClinGen annotation
    |
    v
SpliceAI
    |
    v
Final annotated VCF.gz
```

### CNV branch

```text
Four-column CNV BED
        |
        +-------------------+
        |                   |
        v                   v
     AnnotSV           ClassifyCNV
        |                   |
        +---------+---------+
                  |
                  v
               ISV-CNV
                  |
                  v
      Tool-specific outputs and summary
```

The CNV tools independently analyse the same CNV BED input. AnnotSV output is not used as the input for ClassifyCNV or ISV-CNV.

---

## Main tools

| Tool | Purpose |
|---|---|
| Apptainer | Reproducible container execution |
| bcftools | VCF validation, normalization, querying and annotation |
| bgzip | Block compression of genomic files |
| tabix | Indexing compressed VCF and BED files |
| Ensembl VEP 115 | Gene, transcript, consequence, HGVS, MANE and frequency annotation |
| SnpEff | Independent functional consequence annotation |
| SnpSift | Transfer of ClinVar annotations |
| ClinVar | Clinical significance and disease association |
| ClinGen | Haploinsufficiency and triplosensitivity evidence |
| SpliceAI | Splice gain/loss prediction |
| AnnotSV | Detailed structural-variant annotation |
| ClassifyCNV | ACMG/ClinGen-style CNV classification |
| ISV-CNV | Machine-learning CNV pathogenicity prediction with SHAP explanations |
| bedtools | Genomic interval operations |
| Python / pandas | Data processing and CNV result handling |
| Bash / awk / grep / sed | Automation, filtering and quality control |

---

## Databases and resources

The local workflow uses:

- GRCh38 reference genome;
- Ensembl VEP cache release 115;
- SnpEff `GRCh38.mane.1.2.ensembl` database;
- ClinVar VCF;
- ClinGen dosage-sensitivity resources;
- AnnotSV human annotation package;
- Gene2Phenotype supporting data;
- gnomAD frequency fields available through VEP;
- manual checks using ClinVar, gnomAD, SpliceAI, wInterVar, OMIM and ClinGen.

Large external databases are intentionally not included in the GitHub repository.

---

## Project structure

```text
rare_disease_project/
├── containers/
│   ├── core_tools.sif
│   ├── vep_release115.sif
│   ├── snpeff.sif
│   ├── spliceai.sif
│   ├── annotsv.sif
│   └── isv.sif
│
├── resources/
│   ├── reference/
│   │   ├── hg38.fa
│   │   └── hg38.fa.fai
│   ├── vep_cache/
│   ├── snpeff_data/
│   ├── clinvar/
│   ├── clingen/
│   ├── annotsv_annotations/
│   └── g2p/
│
├── tools/
│   └── ClassifyCNV/
│
├── input/
│   ├── snv/
│   └── cnv/
│
├── pipeline/
│   └── run_rare_disease_pipeline.sh
│
└── results/
    ├── tay_sachs/
    ├── kabuki/
    ├── sotos/
    └── noonan/
```

---

## Input files

### Small variants

Stored in:

```text
input/snv/
```

Expected files:

```text
synthetic_tay_sachs_small_variants_GRCh38_unannotated.vcf
synthetic_kabuki_small_variants_GRCh38_unannotated.vcf
synthetic_sotos_small_variants_GRCh38_unannotated.vcf
synthetic_noonan_small_variants_GRCh38_unannotated.vcf
```

Each initial VCF contains one disease-associated variant and 18 synthetic background variants.

Some background records may be removed during normalization when their `REF` allele does not match GRCh38.

### Copy-number variants

Stored in:

```text
input/cnv/
```

Expected files:

```text
synthetic_tay_sachs_cnv_GRCh38_unannotated.bed
synthetic_kabuki_cnv_GRCh38_unannotated.bed
synthetic_sotos_cnv_GRCh38_unannotated.bed
synthetic_noonan_cnv_GRCh38_unannotated.bed
```

Each final CNV BED contains one disease-associated CNV and nine synthetic distractor CNVs.

BED format:

```text
chromosome    start    end    cnv_type
```

Example:

```text
chr15    72340923    72376014    DEL
```

BED starts are zero-based.

---

## Requirements

### Operating system

Tested on:

- Ubuntu through WSL2;
- Linux command line;
- Apptainer-based environments.

### Required software

At minimum:

```text
Apptainer
Bash
Git
```

The remaining tools are executed through `.sif` containers.

### Recommended hardware

- at least 4 CPU threads;
- at least 8 GB RAM;
- substantial free storage for VEP and annotation databases;
- a Linux-compatible filesystem.

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/Wahid-25/Genomics_project.git
cd Genomics_project
```

### 2. Create the local project structure

```bash
mkdir -p ~/rare_disease_project/{containers,resources,input,pipeline,results,tools}
mkdir -p ~/rare_disease_project/resources/{reference,vep_cache,snpeff_data,clinvar,clingen,annotsv_annotations,g2p}
mkdir -p ~/rare_disease_project/input/{snv,cnv}
```

### 3. Copy scripts and synthetic inputs

```bash
cp pipeline/run_rare_disease_pipeline.sh ~/rare_disease_project/pipeline/
cp input/snv/*.vcf ~/rare_disease_project/input/snv/
cp input/cnv/*.bed ~/rare_disease_project/input/cnv/
```

### 4. Add external resources

The following resources must be added locally:

```text
resources/reference/hg38.fa
resources/reference/hg38.fa.fai
resources/vep_cache/homo_sapiens/115_GRCh38/
resources/snpeff_data/data/GRCh38.mane.1.2.ensembl/
resources/clinvar/clinvar.chr.vcf.gz
resources/clinvar/clinvar.chr.vcf.gz.tbi
resources/clingen/clingen_dosage.hg38.bed.gz
resources/clingen/clingen_dosage.hg38.bed.gz.tbi
resources/annotsv_annotations/AnnotSV_annotations/
```

### 5. Add containers

Place these files in `~/rare_disease_project/containers/`:

```text
core_tools.sif
vep_release115.sif
snpeff.sif
spliceai.sif
annotsv.sif
isv.sif
```

### 6. Add ClassifyCNV

Place ClassifyCNV in:

```text
~/rare_disease_project/tools/ClassifyCNV/
```

The main script should be available as:

```text
tools/ClassifyCNV/ClassifyCNV.py
```

---

## Running the pipeline

Move to the project root:

```bash
cd ~/rare_disease_project
```

Run one disease:

```bash
THREADS=4 bash pipeline/run_rare_disease_pipeline.sh tay_sachs
THREADS=4 bash pipeline/run_rare_disease_pipeline.sh kabuki
THREADS=4 bash pipeline/run_rare_disease_pipeline.sh sotos
THREADS=4 bash pipeline/run_rare_disease_pipeline.sh noonan
```

Run all diseases:

```bash
THREADS=4 bash pipeline/run_rare_disease_pipeline.sh all
```

Replace an existing result:

```bash
THREADS=4 bash pipeline/run_rare_disease_pipeline.sh tay_sachs --force
```

The `--force` option removes and recreates only the selected disease result folder.

---

## Small-variant workflow details

### Normalization

```bash
bcftools norm \
  -f resources/reference/hg38.fa \
  -m -any \
  -c x
```

The `-c x` option excludes records with a REF allele that does not match GRCh38.

### VEP

VEP adds:

- gene and transcript;
- consequence and impact;
- exon or intron number;
- coding and protein positions;
- amino-acid and codon change;
- HGVS coding and protein descriptions;
- canonical and MANE status;
- known identifiers;
- available gnomAD frequencies.

Main field:

```text
INFO/CSQ
```

### SnpEff

SnpEff adds an independent consequence annotation in:

```text
INFO/ANN
```

### ClinVar

Important fields include:

```text
CLNSIG
CLNREVSTAT
CLNDN
CLNDISDB
CLNHGVS
CLNVC
CLNVCSO
GENEINFO
```

### ClinGen

ClinGen adds dosage-sensitivity information, including haploinsufficiency and triplosensitivity evidence.

### SpliceAI

SpliceAI reports:

```text
DS_AG
DS_AL
DS_DG
DS_DL
```

These correspond to acceptor gain, acceptor loss, donor gain and donor loss.

---

## CNV workflow details

### AnnotSV

AnnotSV adds:

- overlapping genes and transcripts;
- cytobands;
- dosage-sensitive regions;
- pathogenic and benign structural-variant overlaps;
- OMIM-related information;
- ranking score;
- ACMG class.

Output:

```text
results/<disease>/cnv/<disease>.AnnotSV.tsv
```

### ClassifyCNV

ClassifyCNV applies ACMG/ClinGen-style evidence scoring and reports:

- CNV classification;
- total score;
- dosage-sensitive genes;
- protein-coding genes;
- evidence categories.

Output folder:

```text
results/<disease>/cnv/<disease>.ClassifyCNV/
```

### ISV-CNV

ISV-CNV provides:

- pathogenicity score;
- threshold-based prediction;
- SHAP feature contributions.

Output:

```text
results/<disease>/cnv/<disease>.ISV_with_SHAP.tsv
```

ISV-CNV is supporting computational evidence and should not replace curated evidence.

---

## Result structure

```text
results/<disease>/
├── final/
├── snv/
├── cnv/
├── logs/
├── work/
└── disease_variant_only/
```

Example final outputs:

```text
results/tay_sachs/final/
├── tay_sachs.final.small_variants.annotated.vcf.gz
├── tay_sachs.final.small_variants.annotated.vcf.gz.tbi
├── tay_sachs.final.cnv.summary.tsv
└── tay_sachs.pipeline_outputs.txt
```

The `disease_variant_only/` folder contains compact disease-specific records for screenshots and comparison.

---

## Main findings

### Tay-Sachs disease

```text
Gene: HEXA
Variant: rs387906309
Coding change: NM_000520.6:c.1274_1277dup
Protein change: p.Tyr427IlefsTer5
Consequence: frameshift
ClinVar: pathogenic
CNV: chr15:72340923-72376014 DEL
```

CNV results:

- AnnotSV: pathogenic;
- ClassifyCNV: likely pathogenic, score `0.90`;
- ISV-CNV: approximately `0.0696`;
- autosomal recessive interpretation requires consideration of both alleles and zygosity.

### Kabuki syndrome

```text
Gene: KMT2D
Variant: rs398123721
Coding change: NM_003482.4:c.14710C>T
Protein change: p.Arg4904Ter
Consequence: stop-gained
ClinVar: pathogenic
CNV: chr12:49018977-49060794 DEL
```

CNV results:

- AnnotSV: ACMG class 5;
- ClassifyCNV: likely pathogenic, score `0.90`;
- ISV-CNV: approximately `0.8899`;
- consistent with KMT2D haploinsufficiency.

### Sotos syndrome

```text
Gene: NSD1
Variant: rs587784148
Coding change: NM_022455.5:c.5431C>T
Protein change: p.Arg1811Ter
Consequence: stop-gained
ClinVar: pathogenic
CNV: chr5:177131797-177300213 DEL
```

CNV results:

- AnnotSV: pathogenic;
- ClassifyCNV: pathogenic, score `1.00`;
- ISV-CNV: approximately `0.9479`;
- strongly supports NSD1 haploinsufficiency.

### Noonan syndrome

```text
Gene: PTPN11
Variant: rs121918459
Coding change: NM_002834.5:c.188A>G
Protein change: p.Tyr63Cys
Consequence: missense
ClinVar: pathogenic
CNV: chr12:112418946-112509918 DUP
```

CNV results:

- AnnotSV: ACMG class 3;
- ClassifyCNV: variant of uncertain significance, score `0.00`;
- ISV-CNV: approximately `0.0194`;
- pathogenic missense variant remains the principal finding;
- the duplication is not supported by established PTPN11 triplosensitivity.

---

## GPT-curated vs pipeline comparison

The repository may include an HTML comparison between:

- GPT-curated educational reference VCFs;
- actual automated pipeline outputs.

Main findings:

- all four disease small-variant targets matched;
- all four disease CNV targets matched;
- the Tay-Sachs insertion was represented differently after normalization but remained biologically equivalent;
- pipeline VCFs contained fewer background variants because REF-mismatching synthetic records were removed;
- pipeline outputs contained richer annotations than the curated reference files.

Recommended file:

```text
GPT_Curated_vs_Pipeline_Comparison.html
```

---

## Quality-control checks

The pipeline verifies:

- input VCF readability;
- four-column CNV BED format;
- numeric start and end coordinates;
- `chr` chromosome prefix;
- `DEL` or `DUP` CNV type;
- required containers and resources;
- successful creation of final files;
- presence of the disease-associated gene;
- presence of VEP, SnpEff, ClinVar, ClinGen and SpliceAI fields;
- non-empty AnnotSV, ClassifyCNV and ISV-CNV outputs;
- successful completion messages in log files.

---

## Common issues and solutions

### REF allele mismatch

```text
Reference allele mismatch
```

Solution:

```bash
bcftools norm -c x
```

This excludes incompatible synthetic background records instead of rewriting them.

### Windows metadata files

Remove files such as `:Zone.Identifier`:

```bash
find . -name '*:Zone.Identifier' -delete
```

### View compressed VCF files

```bash
bcftools view file.vcf.gz > readable_file.vcf
```

### BED and VCF coordinate difference

```text
BED start = VCF POS - 1
BED end = VCF INFO/END
```

### ISV-CNV executable not found

ISV-CNV is run through its Python API inside `isv.sif`.

### Large database storage

The complete gnomAD and dbNSFP datasets are not downloaded. Available VEP frequency fields are used, and important variants are checked manually on gnomAD.

---

## Interpretation rules

1. A rare variant is not automatically pathogenic.
2. A ClinVar match must be checked for the correct disease and review status.
3. A low SpliceAI score does not make a nonsense or frameshift variant benign.
4. A CNV overlapping a disease gene is not automatically pathogenic.
5. The CNV type must match the known disease mechanism.
6. ISV-CNV is supporting machine-learning evidence.
7. Inheritance and zygosity must be considered.
8. Automated results require manual review.
9. Synthetic results must not be treated as clinical diagnoses.

---

## Manual verification resources

The main disease-associated variants were manually checked using:

- ClinVar;
- gnomAD;
- SpliceAI Lookup;
- wInterVar;
- OMIM;
- ClinGen;
- Gene2Phenotype.

---

## GitHub upload recommendations

Upload:

```text
README.md
.gitignore
pipeline scripts
Apptainer definition files
synthetic input VCF and BED files
selected compact outputs
comparison HTML
documentation
small screenshots
```

Do not upload:

```text
*.sif
GRCh38 FASTA
VEP cache
SnpEff database
full ClinVar or ClinGen files
AnnotSV annotation package
licensed databases
temporary files
large intermediate files
real patient data
passwords or secrets
```

---

## Suggested `.gitignore`

```gitignore
# Apptainer images
*.sif

# Large reference and annotation resources
resources/reference/*
resources/vep_cache/*
resources/snpeff_data/*
resources/clinvar/*
resources/clingen/*
resources/annotsv_annotations/*
resources/g2p/*

# Preserve empty directories
!resources/reference/.gitkeep
!resources/vep_cache/.gitkeep
!resources/snpeff_data/.gitkeep
!resources/clinvar/.gitkeep
!resources/clingen/.gitkeep
!resources/annotsv_annotations/.gitkeep
!resources/g2p/.gitkeep

# Temporary pipeline files
results/*/work/*
results/*/logs/*
*.tmp
*.temp

# Python cache
__pycache__/
*.pyc

# Windows and editor metadata
*:Zone.Identifier
.vscode/
.idea/
.DS_Store
Thumbs.db

# Secrets
.env
*.key
*.pem
```

---

## Limitations

- synthetic educational data;
- no real patient phenotype;
- no HPO-based prioritization;
- no family segregation;
- no phasing;
- limited genotype interpretation;
- no wet-lab confirmation;
- no clinical diagnosis;
- database versions may change;
- website results may differ from local resources;
- background variants were not clinically validated as benign.

---

## Future improvements

Possible extensions include:

- phenotype-driven prioritization using HPO terms;
- Nextflow or Snakemake conversion;
- automated HTML reporting;
- MultiQC-style summaries;
- FASTQ/BAM processing and variant calling;
- RNA-seq support;
- structural-variant breakpoint validation;
- automated ACMG evidence extraction;
- Docker support;
- continuous integration testing;
- release-based resource version tracking.

---

## Educational disclaimer

This repository is intended for learning, research training, pipeline development and reproducibility demonstrations.

It is not intended for clinical diagnosis, patient reporting, treatment selection or replacement of a qualified clinical geneticist or molecular pathologist.

---

## Author

**Muhammad Abdul Wahid**

Project:

**Development and Evaluation of an Automated GRCh38 Rare-Disease Variant Annotation Pipeline for SNVs and CNVs**

Instructor:

**Sir Tanveer**

---

## Acknowledgements

This project uses open bioinformatics tools and genomic resources maintained by the scientific community, including Ensembl, NCBI, ClinVar, ClinGen, gnomAD, SnpEff, SpliceAI, AnnotSV, ClassifyCNV and ISV-CNV.

---

## License

The original scripts, synthetic educational inputs and project documentation may be released under the MIT License.

External tools, databases and annotation resources remain subject to their own licences and terms of use.

---

## Final summary

This project demonstrates how unannotated GRCh38 VCF and CNV BED files can be converted into structured and biologically interpretable results through:

```text
normalization
+
functional annotation
+
clinical annotation
+
population evidence
+
splice prediction
+
CNV classification
+
machine-learning prediction
+
manual database verification
```

Reliable interpretation requires agreement between technical correctness, gene relevance, molecular consequence, inheritance, population rarity, clinical evidence, dosage sensitivity, disease mechanism and expert review
