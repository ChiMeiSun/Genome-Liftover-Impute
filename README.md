# Genome-Liftover-Impute

## Overview

This workflow facilitates:

1. **Creation of chain file** – for genome assembly liftover
2. **Liftover** – conversion of genetic data from one genome build to another
3. **Phasing** – preparation of haplotype data
4. **Imputation** – filling in missing genotypes using reference data
5. **Cross-validation** – assessment of imputation accuracy

This workflow is implemented with Snakemake 5.3.0
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.3.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)


## Workflow Instructions

### 1. Setup

Navigate to the folder as your working directory:

```sh
git clone https://github.com/ChiMeiSun/Genome-Liftover-Impute.git
cd Genome-Liftover-Impute
```

Create and activate the conda environment:

```sh
mamba env create -n env_name -f environment.yaml
conda activate env_name
```

### 2. Modify `config.yaml` to configure resources, file paths, and names.

### 3. Running Snakemake

#### Dry-run (check shell commands before execution):

```sh
snakemake -np
```

#### Run a specific rule or output file (e.g. you allow 20 cores):

```sh
snakemake --cores 20 rule_name
snakemake --cores 20 output_file_name
```

#### Run the entire workflow (rule all):

```sh
snakemake --cores 20
```

#### Re-run or re-generate an output:

Snakemake will re-run a rule if any of its input files have been updated or modified. If no input file changes have occurred but you still want to force a rule to be re-run, you can use the --force` flag.

Note that many outputs in this workflow are protected (i.e., read-only or not meant to be overwritten). If you need to regenerate these files despite no updates to the inputs, you must remove the existing protected files first. Once deleted, Snakemake will recognize the missing outputs and re-run the corresponding rules to regenerate them.

---

## Step-by-Step Execution using sample Data

### 1. Create Chain File for Liftover

This process may take more than a day for a complete genome. The mapping is not limited to single chromosomes but includes cross-chromosome mapping.

Modify `config.yaml`:

```yaml
Assembly:
  old: "Data/Assembly/Gallus_gallus.Gallus_gallus-5.0.dna.chromosome.162530.fa"
  new: "Data/Assembly/Gallus_gallus.GRCg6a.dna_rm.chromosome.162530.fa"
Output:
  chain: "Output/Chain/gg5Chr162530_to_GRCg6aChr162530.chain"
```

Run either command:

```sh
snakemake --cores 20 make_chain
snakemake --cores 20 "Output/Chain/gg5Chr162530_to_GRCg6aChr162530.chain"
```

### 2. Liftover using BCFtools

Ensure chromosome labels in the VCF match those in the chain file.

Modify `config.yaml`:

```yaml
Data:
  liftgt: "Data/Liftgt_20k_rename.vcf.gz"
Output:
  genolift: "Output/Liftover/Liftgt_20k_lift.vcf.gz"
  genorej: "Output/Liftover/Liftgt_20k_lift_rej.vcf.gz"
```

Run either command:

```sh
snakemake --cores 20 bcftools_liftover
snakemake --cores 20 "Output/Liftover/Liftgt_20k_lift.vcf.gz"
```

### 3. Imputation with Beagle

#### **Important Notes:**

- Reference and target sets must be aligned to the same genome build (e.g., GRCg6a).
- Multiallelic variants are removed.

  If you want to keep them, go to `Snakefile` and remove ``` bcftools-1.20/bcftools view -M2 -m2 | ``` in rule `beagle_phase_ref` and rule `beagle_phase_gt`
- Contigs absent in the reference panel will be removed from the target set.
- To disable imputation of ungenotyped markers (markers in reference panel but not in the target set), set `impute: "false"` in `config.yaml`.

Modify `config.yaml`:

```yaml
Data:
  ref: "Data/Reference.vcf.gz"
  target: "Data/Target.vcf.gz"
Output:
  ref_phase: "Output/Phase/Reference_phase.vcf.gz"
  target_phase: "Output/Phase/Target_phase.vcf.gz"
  target_impute: "Output/Impute/Target_phase_impute.vcf.gz"
Params:
  impute: "true"
```

### 4. Cross-validation for Imputation Accuracy

This step masks variants in folds and compares results with known values.

Modify `config.yaml`:

```yaml
Impacc:
  dir: "ImputationAccuracy"
  rep: [1, 2] # Repetitions
  fold: [5, 10] # Fold values
```

**Outputs:**

- Accuracy reports in both **PDF** and **TXT** format
- Metrics include IQS, Pearson correlation, dosage R², and concordance

### Output Column Descriptions

| Column        | Description                                                              |
| ------------- | ------------------------------------------------------------------------ |
| rep           | Repetition                                                               |
| fold          | Fold number                                                              |
| fsnp          | Number of SNPs in fold                                                   |
| fsnps\_imp    | Number of SNPs remained after imputation                                 |
| fsnps\_idk    | Number of SNPs new in imputation (not found in gt)                       |
| fsnps\_ana    | Number of SNPs overlapped in gt and after imputation (used for analysis) |
| fsnps\_miss   | = fsnp - fsnps\_ana                                                        |
| prop\_RAmatch | Proportion of variants that match both Ref and Alt                       |
| m\_iqs        | Mean Imputation Quality Score (IQS) across variants                      |
| m\_p0         | Mean concordance                                                         |
| m\_cor        | Mean Pearson’s correlation                                               |
| min\_cor      | Minimum Pearson’s correlation                                            |
| m\_DR2        | Mean DR2 (dosage r-square)                                               |

---




### Sample Data Preparation (Provided in folder Data)

A Fast Example of Chicken 20k SNPs on chr 16,25,30 aligned to Gallus_gallus-5.0

data from:

Geibel, J., Reimer, C., Weigend, S., Weigend, A., Pook, T., & Simianer, H. (2020). Data for "How Array Design Creates SNP Ascertainement Bias" [Data set]. Zenodo. https\://doi.org/10.5281/zenodo.4320456                                                                                                                                                                  &#x20;
### 1. Download Reference Genome Assemblies

```sh
mkdir -p Data/Assembly/
wget -P Data/Assembly/ https://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna_rm.chromosome.{16,25,30}.fa.gz
gzip -d Data/Assembly/*.fa.gz
wget -P Data/Assembly/ https://42basepairs.com/download/web/ensembl/release-92/fasta/gallus_gallus/dna/Gallus_gallus.Gallus_gallus-5.0.dna.chromosome.{16,25,30}.fa.gz
gzip -d Data/Assembly/*.fa.gz
cat Data/Assembly/Gallus_gallus.Gallus_gallus-5.0.dna.chromosome.{16,25,30}.fa > Data/Assembly/Gallus_gallus.Gallus_gallus-5.0.dna.chromosome.162530.fa
cat Data/Assembly/Gallus_gallus.GRCg6a.dna_rm.chromosome.{16,25,30}.fa > Data/Assembly/Gallus_gallus.GRCg6a.dna_rm.chromosome.162530.fa
```

### 2. Download and Prepare Data for Liftover

```sh
wget -P Data/ https://zenodo.org/records/4320456/files/IndandPool_all.raw.vcf.gz?download=1
wget -P Data/ https://zenodo.org/records/4320456/files/IndandPool_all.raw.vcf.gz.tbi?download=1

bcftools view -r chr16 Data/IndandPool_all.raw.vcf.gz -Oz -o Data/IndandPool_all_chr16.vcf.gz
bcftools view -r chr25 Data/IndandPool_all.raw.vcf.gz -Oz -o Data/IndandPool_all_chr25.vcf.gz
bcftools view -r chr30 Data/IndandPool_all.raw.vcf.gz -Oz -o Data/IndandPool_all_chr30.vcf.gz

bcftools index -f Data/IndandPool_all_chr16.vcf.gz
bcftools index -f Data/IndandPool_all_chr25.vcf.gz
bcftools index -f Data/IndandPool_all_chr30.vcf.gz

bcftools concat -a \
  Data/IndandPool_all_chr16.vcf.gz \
  Data/IndandPool_all_chr25.vcf.gz \
  Data/IndandPool_all_chr30.vcf.gz \
  -Oz -o Data/Liftgt.vcf.gz
```

Subset 20k SNPs:

```sh
bcftools query -f '%CHROM\t%POS\n' Data/Liftgt.vcf.gz | shuf -n 20000 > Data/random_snps_20k.txt
bcftools view -T Data/random_snps_20k.txt Data/Liftgt.vcf.gz -Oz -o Data/Liftgt_20k.vcf.gz
```

Rename chromosomes to match chain file:

```sh
awk '{print $1, substr($1, 4)}' Data/random_snps_20k.txt > Data/chr.txt
bcftools annotate --rename-chrs Data/chr.txt Data/Liftgt_20k.vcf.gz -Oz -o Data/Liftgt_20k_rename.vcf.gz
rm Data/chr.txt
```

### 3. Prepare Data for Imputation

Reference set (chr 16, 25):

```sh
bcftools view -r 16,25 Output/Liftover/Liftgt_20k_lift.vcf.gz -Oz -o Data/Reference.vcf.gz
bcftools-1.20/bcftools query -f '%CHROM' Data/Reference.vcf.gz | wc -l
```

Target set (random 2k variants from chr 16, 25, 30):

```sh
bcftools query -f '%CHROM\t%POS\n' Output/Liftover/Liftgt_20k_lift.vcf.gz | shuf -n 2000 > Data/random_snps_2k.txt
bcftools view -T Data/random_snps_2k.txt Data/tmp.vcf.gz -Oz -o Data/Target.vcf.gz
bcftools-1.20/bcftools query -f '%CHROM' Data/Target.vcf.gz | uniq -c
```

### 4. Check output
```sh
bcftools-1.20/bcftools query -f '%CHROM' Output/Phase/Reference_phase.vcf.gz | uniq -c
bcftools-1.20/bcftools query -f '%CHROM' Output/Phase/Target_phase.vcf.gz | uniq -c
bcftools-1.20/bcftools query -f '%CHROM' Output/Impute/Target_phase_impute.vcf.gz | uniq -c
```

Chr 30 was removed since it did not present in the reference panel





