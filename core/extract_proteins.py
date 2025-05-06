import os
import subprocess
import sys


def extract_proteins(genome_file, gff_file, out_file, logger):
        """
        Extract protein sequences from genome using GFF annotation.
        
        Uses gffread to extract protein-coding sequences from the genome 
        based on the annotation, then formats headers for downstream analysis.
        
        Parameters:
        -----------
        genome_file : str
            Path to the genome FASTA file
        gff_file : str
            Path to the genome annotation GFF file
        out_file : str
            Path for the output protein FASTA file
        logger : logging.Logger
            Logger object for reporting
            
        Returns:
        --------
        None
            Extracted proteins are written to the specified output file
        """
        # Extract proteins
        gffread_cmd = ['gffread', '-S', '--no-pseudo', '-y', out_file, '-g', genome_file, gff_file]

        # Clean up protein headers
        sed_cmd = ['sed', '-E', 's/>rna-([^[:space:]]+).*/>\\1/', out_file, '-i']

        try:
            logger.debug(f"Running command: {' '.join(gffread_cmd)}")
            subprocess.run(gffread_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.debug(f"Running command: {' '.join(sed_cmd)}")
            subprocess.run(sed_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.info(f"Extracted proteins saved to {out_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error extracting proteins: {e}")
            logger.error(e.stderr.decode())
            sys.exit(1)
