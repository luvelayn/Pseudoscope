"""
Pseudoscope - A Pseudogene Identification Pipeline

A command-line tool for automated identification and annotation of pseudogenes
in genomic sequences, based on homology search and structural analysis.
"""
from core._setup_logging import _setup_logging
from core._check_dependencies import _check_dependencies
from core.extract_proteins import extract_proteins
from core.mask_genes import mask_genes
from core.run_tblastn import run_tblastn
from core.filter_hits import filter_hits
from core.merge_hits import merge_hits
from core.run_exonerate import run_exonerate

import os
import argparse

__version__ = '0.1.0'

class Pseudoscope:
    """Main class for the pseudogene identification pipeline"""

    def __init__(self, genome_file, gff_file, protein_file=None, output_dir="./pseudoscope_output",
                 evalue=0.01, max_intron_length=10000, threads=8):
        """
        Initialize the Pseudoscope pipeline
        
        Args:
            genome_file (str): Path to genome FASTA file
            gff_file (str): Path to genome annotation GFF file
            protein_file (str, optional): Path to protein FASTA file
            output_dir (str, optional): Directory for output files
            evalue (float, optional): E-value threshold for BLAST
            max_intron_length (int, optional): Maximum allowed intron length
            threads (int, optional): Number of threads to use
        """
        self.genome_file = os.path.abspath(genome_file)
        self.gff_file = os.path.abspath(gff_file)
        self.protein_file = os.path.abspath(protein_file) if protein_file else None
        self.output_dir = os.path.abspath(output_dir)
        self.evalue = evalue
        self.max_intron_length = max_intron_length
        self.threads = threads
        
        # Set up working directory
        self.temp_dir = os.path.join(self.output_dir, "tmp")
        self.mask_genes_dir = os.path.join(self.temp_dir, "mask_genes_out")
        self.extract_proteins_dir = os.path.join(self.temp_dir, "extract_proteins_out")
        self.blast_dir = os.path.join(self.temp_dir, "blast_out")
        self.filter_and_merge_dir = os.path.join(self.temp_dir, "filtered_and_merged_hits")
        self.exonerate_dir = os.path.join(self.temp_dir, "exonerate_out")
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.temp_dir, exist_ok=True)
        os.makedirs(self.mask_genes_dir, exist_ok=True)
        os.makedirs(self.extract_proteins_dir, exist_ok=True)
        os.makedirs(self.blast_dir, exist_ok=True)
        os.makedirs(self.filter_and_merge_dir, exist_ok=True)
        os.makedirs(self.exonerate_dir, exist_ok=True)
        
        # Set up logging
        self.logger = _setup_logging(self.output_dir)
        
        # Check for required tools
        _check_dependencies(self.logger)
        
        # # Final output paths
        # self.pseudogenes_gff = os.path.join(self.output_dir, "pseudogenes.gff")
        # self.pseudogenes_tsv = os.path.join(self.output_dir, "pseudogenes.tsv")
        # self.pseudogene_fasta = os.path.join(self.output_dir, "pseudogenes.fa")

    def run(self):
        """Run the complete pseudogene identification pipeline"""
        self.logger.info(f"Starting Pseudoscope v{__version__}")
        self.logger.info(f"Genome file: {self.genome_file}")
        self.logger.info(f"GFF file: {self.gff_file}")
        self.logger.info(f"Output directory: {self.output_dir}")

        # Step 1: Extract proteins if not provided
        if not self.protein_file:
            self.logger.info("Extracting proteins from genome annotation")
            self.protein_file = os.path.join(self.extract_proteins_dir, "proteins.fa")
            extract_proteins(self.genome_file, self.gff_file, self.extract_proteins_dir, self.logger)
        else:
            self.logger.info(f"Using provided protein file: {self.protein_file}")

        # Step 2: Mask functional genes
        self.logger.info("Masking functional genes in the genome")
        genes_masked__genome = mask_genes(self.gff_file, self.genome_file, self.mask_genes_dir, self.logger)

        # Step 3: Homology search
        self.logger.info("Running homology search with tblastn")
        blast_out = run_tblastn(genes_masked__genome, self.protein_file, self.evalue,
                                     self.threads, self.max_intron_length, self.blast_dir, self.logger)
        # blast_out = os.path.join(self.blast_dir, "tblastn_out/tblastn_results.tsv")

        # Step 4: Filter and merge hits
        self.logger.info("Filtering and merging BLAST hits")
        hits = filter_hits(blast_out, self.filter_and_merge_dir, self.logger)
        pseudogenes = merge_hits(hits, self.protein_file, self.max_intron_length, self.filter_and_merge_dir, self.logger)

        # Step 5: Precise re-alignment with exonerate
        self.logger.info("Running precise re-alignment with exonerate")
        run_exonerate(pseudogenes, self.protein_file, genes_masked__genome, self.exonerate_dir, self.logger)

        # # Step 6: Classify pseudogenes
        # self.logger.info("Classifying pseudogenes")
        # classified_pseudogenes = self.classify_pseudogenes()

        # # Step 7: Generate final report
        # self.logger.info("Generating final pseudogene annotations")
        # self.generate_output(classified_pseudogenes)

        self.logger.info("Pseudoscope pipeline completed successfully")
        self.logger.info(f"Results can be found in: {self.output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Pseudoscope - A Pseudogene Identification Pipeline')
    parser.add_argument('-g', '--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('-a', '--annotation', required=True, help='Path to genome annotation GFF file')
    parser.add_argument('-p', '--proteins', help='Path to protein FASTA file (optional)')
    parser.add_argument('-o', '--output', default='./pseudoscope_output', help='Output directory')
    parser.add_argument('-e', '--evalue', type=float, default=0.01, help='E-value threshold for BLAST')
    parser.add_argument('-m', '--max-intron', type=int, default=10000, help='Maximum allowed intron length')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads to use')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()
    
    pipeline = Pseudoscope(
        genome_file=args.genome,
        gff_file=args.annotation,
        protein_file=args.proteins,
        output_dir=args.output,
        evalue=args.evalue,
        max_intron_length=args.max_intron,
        threads=args.threads
    )
    
    pipeline.run()

if __name__ == "__main__":
    main()