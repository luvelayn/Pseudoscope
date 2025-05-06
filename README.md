# Pseudoscope

A comprehensive pseudogene identification pipeline for genomic analysis.

## Overview

Pseudoscope is a command-line tool for automated identification and annotation of pseudogenes in genomic sequences, based on homology search and structural analysis. It identifies potential pseudogenes by comparing protein sequences to a genome, identifying homologous regions outside of known gene locations.

## Features

- Automated gene masking
- Protein homology search using TBLASTN
- Precise re-alignment using TFASTY
- Intron-exon structure refinement
- Pseudogene classification
- Comprehensive results reporting

## Dependencies

Pseudoscope requires the following software to be installed:

- BLAST+ (makeblastdb, tblastn, dustmasker)
- bedtools
- gffread
- FASTA suite (specifically tfasty36)
- Biopython

## Installation

### Setting up a conda environment

The easiest way to install Pseudoscope and its dependencies is through conda:

```bash
# Create a new conda environment
conda create -n pseudoscope python=3.9

# Activate the environment
conda activate pseudoscope

# Install dependencies
conda install -c bioconda blast bedtools gffread fasta biopython

# Clone the repository
git clone https://github.com/yourusername/pseudoscope.git
cd pseudoscope

# Install the package
pip install -e .
```

## Usage

### Basic usage

```bash
python pseudoscope.py -g genome.fasta -a annotation.gff -out output_directory
```

### Command-line arguments

```
pseudoscope.py [-h] -g GENOME -a ANNOTATION [-p PROTEINS] [-out OUTPUT]
               [-mi MAX_INTRON] [-e EVALUE] [-cov COVERAGE] [-id IDENTITY]
               [-t THREADS] [-v]

Arguments:
  -g, --genome       Path to genome FASTA file
  -a, --annotation   Path to genome annotation GFF file
  -p, --proteins     Path to protein FASTA file (optional)
  -out, --output     Output directory (default: ./pseudoscope_output)
  -mi, --max-intron  Maximum allowed intron length (default: 10000)
  -e, --evalue       E-value threshold for filtering (default: 1e-5)
  -cov, --coverage   Coverage threshold for filtering (default: 0.05)
  -id, --identity    Identity threshold for filtering (default: 0.2)
  -t, --threads      Number of threads to use (default: 1)
  -v, --version      Show program's version number and exit
```

## Pipeline Workflow

Pseudoscope operates through the following steps:

1. **Extract proteins**: If no protein file is provided, extract proteins from the genome annotation.
2. **Mask genes**: Mask functional genes in the genome to prevent false positives.
3. **Homology search**: Run TBLASTN to find protein homology regions in the masked genome.
4. **Filter hits**: Filter BLAST hits based on dynamic identity thresholds.
5. **Create exons**: Merge overlapping or closely located hits into exons.
6. **Create clusters**: Group exons by protein, chromosome, and strand to form pseudogene candidates.
7. **Re-alignment**: Run precise re-alignment with TFASTY for each pseudogene candidate.
8. **Refine structure**: Refine intron-exon structure based on canonical splice sites.
9. **Classify**: Classify pseudogenes based on their molecular characteristics.
10. **Generate output**: Create comprehensive reports of identified pseudogenes.

## Output

The results will be saved in the specified output directory with the following structure:

```
output_directory/
├── logs/
│   └── pseudoscope.log
├── results/
│   ├── pseudogenes.bed
│   ├── pseudogenes.gff
│   └── pseudogenes.fasta
```

## Example

```bash
# Activate the conda environment
conda activate pseudoscope

# Run Pseudoscope on your data
python pseudoscope.py \
  -g example/genome.fasta \
  -a example/annotation.gff \
  -out example_results \
  -t 8 \
```
