import os
import subprocess
import sys


def mask_genes(gff_file, genome_file, out_dir, logger):
    """
    Mask coding gene regions in the genome using bedtools.

    Extracts protein-coding gene annotations and uses bedtools to hard-mask
    these regions in the genome sequence to prevent false positive identifications.

    Parameters:
    -----------
    gff_file : str
        Path to the genome annotation GFF file
    genome_file : str
        Path to the genome FASTA file
    out_dir : str
        Directory to store output files
    logger : logging.Logger
        Logger object for reporting

    Returns:
    --------
    str
        Path to the masked genome file
    """
    # Extract protein coding genes from GFF annotation file
    logger.debug("Extracting protein coding genes from GFF file")

    genes_file = os.path.join(out_dir, "genes.gff")

    with open(gff_file, "r") as infile, open(genes_file, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            if fields[2] == "gene":
                attributes = fields[-1]
                if "protein_coding" in attributes:
                    outfile.write(line)

    # Hard-mask genes using bedtools
    genes_masked__genome = os.path.join(out_dir, "genes_masked__genome.gff")

    cmd = [
        "bedtools",
        "maskfasta",
        "-fi",
        genome_file,
        "-bed",
        genes_file,
        "-fo",
        genes_masked__genome,
    ]

    logger.debug(f"Running command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info(f"Functional genes masked successfully")
        return genes_masked__genome
    except subprocess.CalledProcessError as e:
        logger.error(f"Error masking genome: {e}")
        logger.error(e.stderr.decode())
        sys.exit(1)
