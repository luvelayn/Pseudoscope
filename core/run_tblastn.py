import os
import subprocess
import sys


def _format_blast_out(blast_out, logger):
        """
        Format BLAST output for consistency.
        
        Adds descriptive headers and normalizes sequence identifiers
        to ensure consistent downstream processing.
        
        Parameters:
        -----------
        blast_out : str
            Path to the BLAST output file
        logger : logging.Logger
            Logger object for reporting
            
        Returns:
        --------
        None
            Modifications are made in-place to the BLAST output file
        """
        # Add description line
        with open(blast_out, 'r') as f:
            existing_lines = f.readlines()

        desc_line = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\n"
        updated_lines = [desc_line] + existing_lines

        with open(blast_out, 'w') as f:
            f.writelines(updated_lines)

        # Format sequence identifiers
        cmd = ['sed', '-E', 's/ref\|([^|]+)\|/\\1/', blast_out, '-i']
        
        logger.debug(f"Running command: {' '.join(cmd)}")

        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.info(f"BLAST output formatted")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error formatting BLAST output: {e}")
            logger.error(e.stderr.decode())
            sys.exit(1)

def run_tblastn(genes_masked__genome, protein_file, threads, max_intron_length, out_dir, logger):
        """
        Run tblastn to align proteins to genome.
        
        Prepares a BLAST database with soft-masking for low-complexity regions,
        then runs TBLASTN with parameters optimized for pseudogene detection.
        
        Parameters:
        -----------
        genes_masked__genome : str
            Path to the genome file with coding genes masked
        protein_file : str
            Path to the protein FASTA file
        threads : int
            Number of threads to use for BLAST
        max_intron_length : int
            Maximum allowed intron length for TBLASTN
        out_dir : str
            Directory to store output files
        logger : logging.Logger
            Logger object for reporting
            
        Returns:
        --------
        str
            Path to the formatted TBLASTN output file
        """
        # Create masking information for low-compexity regions soft-masking
        mask_dir = os.path.join(out_dir, "mask_info")
        os.makedirs(mask_dir, exist_ok=True)
        mask_info = os.path.join(mask_dir, "genome_mask_info.asnb")

        mask_cmd = ['dustmasker',
                    '-in', genes_masked__genome,
                    '-infmt', 'fasta',
                    '-parse_seqids',
                    '-outfmt', 'maskinfo_asn1_bin',
                    '-out', mask_info]
        
        logger.debug(f"Running command: {' '.join(mask_cmd)}")
        
        try:
            subprocess.run(mask_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.info("Masking info file created successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating masking info file: {e}")
            logger.error(e.stderr.decode())
            sys.exit(1)
    
        # Create BLAST database
        makeblastdb_dir = os.path.join(out_dir, "makeblastdb_out")
        os.makedirs(mask_dir, exist_ok=True)
        blast_db = os.path.join(makeblastdb_dir, "genome_blast_db")

        makeblastdb_cmd = ['makeblastdb', 
                           '-in', genes_masked__genome, 
                           '-dbtype', 'nucl',
                           '-parse_seqids', 
                           '-mask_data', mask_info,
                           '-out', blast_db]
        
        logger.debug(f"Running command: {' '.join(makeblastdb_cmd)}")
        
        try:
            subprocess.run(makeblastdb_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.info("BLAST database created successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating BLAST database: {e}")
            logger.error(e.stderr.decode())
            sys.exit(1)
            
        # Run tblastn
        tblastn_dir = os.path.join(out_dir, "tblastn_out")
        os.makedirs(tblastn_dir, exist_ok=True)
        blast_out = os.path.join(tblastn_dir, "tblastn_results.tsv")

        tblastn_cmd = [
            'tblastn', 
            '-query', protein_file,
            '-db', blast_db,
            '-out', blast_out,
            '-evalue', '0.01',
            '-outfmt', '6',
            '-num_threads', str(threads),
            '-seg', 'yes',
            '-soft_masking', 'true',
            '-db_soft_mask', 'dust',
            '-word_size', '3',
            '-gapextend', '2',
            '-max_intron_length', str(max_intron_length)
        ]
        
        logger.debug(f"Running command: {' '.join(tblastn_cmd)}")
        
        try:
            subprocess.run(tblastn_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logger.info(f"tblastn completed successfully")
            _format_blast_out(blast_out, logger)
            return blast_out
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running tblastn: {e}")
            logger.error(e.stderr.decode())
            sys.exit(1)