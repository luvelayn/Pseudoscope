"""
Pseudoscope - A Pseudogene Identification Pipeline

A command-line tool for automated identification and annotation of pseudogenes
in genomic sequences, based on homology search and structural analysis.
"""

import shutil
from utils._setup_logging import _setup_logging
from utils._check_dependencies import _check_dependencies
from core.extract_proteins import extract_proteins
from core.mask_genes import mask_genes
from core.run_tblastn import run_tblastn
from core.filter_hits import filter_hits
from core.create_exons import create_exons
from core.create_clusters import create_clusters
from core.run_tfasty import run_tfasty
from core.refine_exon_structure import refine_exon_structure
from core.classify_pseudogenes import classify_pseudogenes
from core.generate_output import generate_output

import os
import argparse
from Bio import SeqIO

__version__ = "0.1.0"


class Pseudoscope:
    """Main class for the pseudogene identification pipeline"""

    def __init__(
        self,
        genome_file,
        gff_file,
        output_dir,
        max_intron_length,
        evalue,
        coverage,
        identity,
        threads,
        protein_file=None,
    ):
        """
        Initialize the Pseudoscope pipeline

        Args:
            genome_file (str): Path to genome FASTA file
            gff_file (str): Path to genome annotation GFF file
            output_dir (str): Directory for output files
            max_intron_length (int): Maximum allowed intron length
            evalue (float): E-value threshold for filtering
            coverage (float): Coverage threshold for filtering
            identity (float): Identity threshold for filtering
            threads (int): Number of threads to use
            protein_file (str, optional): Path to protein FASTA file
        """
        self.genome_file = os.path.abspath(genome_file)
        self.gff_file = os.path.abspath(gff_file)
        self.protein_file = os.path.abspath(protein_file) if protein_file else None
        self.output_dir = os.path.abspath(output_dir)
        self.max_intron_length = max_intron_length
        self.evalue = evalue
        self.coverage = coverage
        self.identity = identity
        self.threads = threads

        # Set up working directory
        self.temp_dir = os.path.join(self.output_dir, "tmp")
        self.logs_dir = os.path.join(self.output_dir, "logs")
        self.results_dir = os.path.join(self.output_dir, "results")
        self.mask_genes_dir = os.path.join(self.temp_dir, "mask_genes_out")
        self.extract_proteins_dir = os.path.join(self.temp_dir, "extract_proteins_out")
        self.blast_dir = os.path.join(self.temp_dir, "blast_out")
        self.filter_and_merge_dir = os.path.join(
            self.temp_dir, "filtered_and_merged_hits"
        )
        self.tfasty_dir = os.path.join(self.temp_dir, "tfasty_out")
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.temp_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(self.mask_genes_dir, exist_ok=True)
        os.makedirs(self.extract_proteins_dir, exist_ok=True)
        os.makedirs(self.blast_dir, exist_ok=True)
        os.makedirs(self.filter_and_merge_dir, exist_ok=True)
        os.makedirs(self.tfasty_dir, exist_ok=True)

        # Set up logging
        self.logger = _setup_logging(self.logs_dir)

        # Check for required tools
        _check_dependencies(self.logger)

    def run(self):
        """
        Run the complete pseudogene identification pipeline.

        Executes all steps in sequence:
        1. Extract proteins from genome annotation (if not provided)
        2. Mask functional genes in the genome
        3. Run homology search with tblastn
        4. Filter and merge BLAST hits
        5. Run precise re-alignment with tfasty
        6. Refine exon structure
        7. Classify pseudogenes
        8. Generate final reports

        Returns:
        --------
        None
            Results are written to the specified output directory
        """
        self.logger.info(f"Starting Pseudoscope v{__version__}")
        self.logger.info(f"Genome file: {self.genome_file}")
        self.logger.info(f"GFF file: {self.gff_file}")
        self.logger.info(f"Output directory: {self.output_dir}")

        # Step 1: Extract proteins if not provided
        if not self.protein_file:
            self.logger.info("Extracting proteins from genome annotation...")
            self.protein_file = os.path.join(self.extract_proteins_dir, "proteins.fa")
            extract_proteins(
                self.genome_file, self.gff_file, self.extract_proteins_dir, self.logger
            )
        else:
            self.logger.info(f"Using provided protein file: {self.protein_file}")

        # Load protein sequences
        self.protein_seqs = {}
        for record in SeqIO.parse(self.protein_file, "fasta"):
            self.protein_seqs[record.id] = str(record.seq)

        # Step 2: Mask functional genes
        self.logger.info("Masking functional genes in the genome...")
        genes_masked__genome = mask_genes(
            self.gff_file, self.genome_file, self.mask_genes_dir, self.logger
        )

        # Load genome sequences
        self.genome_seqs = {}
        for record in SeqIO.parse(genes_masked__genome, "fasta"):
            self.genome_seqs[record.id] = str(record.seq).upper()

        # Step 3: Homology search
        self.logger.info("Running homology search with tblastn...")
        blast_out = run_tblastn(
            genes_masked__genome,
            self.protein_file,
            self.threads,
            self.max_intron_length,
            self.blast_dir,
            self.logger,
        )

        # Step 4: Filter and merge hits
        self.logger.info("Filtering and merging BLAST hits...")
        filtered_hits_file = filter_hits(
            blast_out, self.filter_and_merge_dir, self.logger
        )
        exons_file = create_exons(
            filtered_hits_file, self.filter_and_merge_dir, self.logger
        )
        exon_clusters = create_clusters(
            exons_file, self.protein_seqs, self.max_intron_length, self.logger
        )

        # Step 5: Precise re-alignment with tfasty
        self.logger.info("Running precise re-alignment with tfasty and filtering...")
        tfasty_pseudogenes = run_tfasty(
            exon_clusters,
            self.protein_seqs,
            self.genome_seqs,
            self.evalue,
            self.coverage,
            self.identity,
            self.tfasty_dir,
            self.logger,
        )

        # Step 6: Refine intron-exon structure
        self.logger.info(
            "Refining intron-exon structure based on canonical splice sites..."
        )
        refined_pseudogenes = refine_exon_structure(
            tfasty_pseudogenes, self.genome_seqs, self.logger
        )

        # Step 7: Classify pseudogenes
        self.logger.info("Classifying pseudogenes...")
        classified_pseudogenes = classify_pseudogenes(
            refined_pseudogenes, self.genome_seqs, self.logger
        )

        # Step 8: Generate final reports and clean tmp directory
        self.logger.info("Generating final pseudogene reports...")
        generate_output(classified_pseudogenes, self.genome_seqs, self.results_dir)
        shutil.rmtree(self.temp_dir)

        self.logger.info("Pseudoscope pipeline completed successfully")
        self.logger.info(f"Results can be found in: {self.results_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Pseudoscope - A Pseudogene Identification Pipeline"
    )
    parser.add_argument(
        "-g", "--genome", required=True, help="Path to genome FASTA file"
    )
    parser.add_argument(
        "-a", "--annotation", required=True, help="Path to genome annotation GFF file"
    )
    parser.add_argument(
        "-p", "--proteins", help="Path to protein FASTA file (optional)"
    )
    parser.add_argument(
        "-out", "--output", default="./pseudoscope_output", help="Output directory"
    )
    parser.add_argument(
        "-mi",
        "--max-intron",
        type=int,
        default=10000,
        help="Maximum allowed intron length (default: 10000)",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=1e-5,
        help="E-value threshold for filtering (default: 1e-5)",
    )
    parser.add_argument(
        "-cov",
        "--coverage",
        type=float,
        default=0.05,
        help="Coverage threshold for filtering (default: 0.05)",
    )
    parser.add_argument(
        "-id",
        "--identity",
        type=float,
        default=0.2,
        help="Identity threshold for filtering (default: 0.2)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=4, help="Number of threads to use (default: 4)"
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    args = parser.parse_args()

    pipeline = Pseudoscope(
        genome_file=args.genome,
        gff_file=args.annotation,
        output_dir=args.output,
        max_intron_length=args.max_intron,
        evalue=args.evalue,
        coverage=args.coverage,
        identity=args.identity,
        threads=args.threads,
        protein_file=args.proteins,
    )

    pipeline.run()


if __name__ == "__main__":
    main()
