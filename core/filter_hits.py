import os
import sys


def _calculate_identity_threshold(alignment_length):
    """
    Calculate dynamic identity threshold based on alignment length.

    Implements a sliding scale where shorter alignments require higher
    identity percentages to be considered valid.

    Parameters:
    -----------
    alignment_length : int
        Length of the alignment

    Returns:
    --------
    float
        Identity threshold percentage
    """
    threshold = max(75 - ((alignment_length - 20) * 0.6875), 20)
    return threshold


def filter_hits(tblastn_file, out_dir, logger):
    """
    Filter BLAST hits based on identity threshold and add strand information.

    Applies dynamic identity thresholds based on alignment length, adds
    strand information, and writes filtered hits to output file.

    Parameters:
    -----------
    tblastn_file : str
        Path to the TBLASTN output file
    out_dir : str
        Directory to store output files
    logger : logging.Logger
        Logger object for reporting

    Returns:
    --------
    str
        Path to the filtered hits file
    """
    filtered_hits = []

    # Load BLAST results
    try:
        with open(tblastn_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 12:
                    continue

                qseqid = fields[0]  # Query protein ID
                sseqid = fields[1]  # Subject genome sequence ID
                pident = float(fields[2])  # Percent identity
                length = int(fields[3])  # Alignment length
                qstart = int(fields[6])  # Query start position
                qend = int(fields[7])  # Query end position
                sstart = int(fields[8])  # Subject start position
                send = int(fields[9])  # Subject end position
                evalue = float(fields[10])  # E-value

                # Calculate dynamic identity threshold
                threshold = _calculate_identity_threshold(length)

                # Keep hit if it passes identity threshold
                if pident >= threshold:
                    strand = "+"

                    if sstart > send:
                        sstart, send = send, sstart
                        strand = "-"

                    filtered_hits.append(
                        {
                            "qseqid": qseqid,
                            "sseqid": sseqid,
                            "pident": pident,
                            "length": length,
                            "qstart": qstart,
                            "qend": qend,
                            "sstart": sstart,
                            "send": send,
                            "strand": strand,
                            "evalue": evalue,
                        }
                    )

        logger.info(f"BLAST hits filtered: {len(filtered_hits)} hits retained")

        # Save filtered hits to file
        filtered_hits_file = os.path.join(out_dir, "filtered_hits.tsv")

        with open(filtered_hits_file, "w") as f:
            f.write(
                "qseqid\tsseqid\tpident\tlength\tqstart\tqend\tsstart\tsend\tstrand\tevalue\n"
            )
            for hit in filtered_hits:
                f.write(
                    f"{hit['qseqid']}\t{hit['sseqid']}\t{hit['pident']}\t{hit['length']}\t{hit['qstart']}\t"
                    f"{hit['qend']}\t{hit['sstart']}\t{hit['send']}\t{hit['strand']}\t{hit['evalue']}\n"
                )

        return filtered_hits_file

    except Exception as e:
        logger.error(f"Error filtering BLAST hits: {e}")
        sys.exit(1)
