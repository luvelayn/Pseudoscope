def _find_sites(sequence, motif):
    """
    Find all start positions of a motif in a sequence.

    Parameters:
    -----------
    sequence : str
        DNA sequence to search in.
    motif : str
        Motif to search for.

    Returns:
    --------
    list
        # List of starting positions (0-based) where the motif is found.
    """
    seq_lower = sequence.lower()
    motif_lower = motif.lower()

    positions = []
    start = 0

    while True:
        pos = seq_lower.find(motif_lower, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1

    return positions


def _reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.

    Converts A->T, T->A, G->C, C->G and reverses the resulting sequence.
    Preserves case and handles ambiguous bases (N).

    Parameters:
    -----------
    seq : str
        DNA sequence to reverse complement.

    Returns:
    --------
    str
        Reverse complemented DNA sequence.
    """
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "N": "N",
        "n": "n",
    }
    return "".join(complement.get(base, base) for base in reversed(seq))


def refine_exon_structure(pseudogenes, genome_seqs, logger):
    """
    Refine the intron-exon structure of pseudogenes by checking for canonical splice sites.

    Examines the genomic sequence at putative splice junctions (intron-exon boundaries)
    to confirm the presence of canonical splice site motifs (GT-AG). Exons without
    proper splice signals are merged to improve accuracy of the pseudogene structure.

    Parameters:
    -----------
    pseudogenes : list of dict
        List of pseudogene dictionaries, each containing exon information.
    genome_seqs : dict
        Dictionary of chromosome sequences keyed by chromosome name.
    logger : logging.Logger
        Logger object for reporting progress and statistics.

    Returns:
    --------
    list of dict
        Refined pseudogenes with updated exon structure and recalculated statistics.
        Exons without canonical splice sites are merged.
    """
    # Define canonical splice site patterns
    # GT-AG rule: 5' donor site (GT) and 3' acceptor site (AG)
    DONOR_SITE = "GT"  # 5' splice site starts with GT
    ACCEPTOR_SITE = "AG"  # 3' splice site ends with AG

    # Context size to extract around splice sites for better detection
    SPLICE_CONTEXT = 4  # ±SPLICE_CONTEXT bp around the splice site

    refined_pseudogenes = []
    merged_count = 0

    for pg_idx, pg in enumerate(pseudogenes):
        chrom = pg["chrom"]
        strand = pg["strand"]
        exons = pg["exons"]

        # Sort exons by position
        exons.sort(key=lambda x: x["start"])

        # Skip if only one exon
        if len(exons) == 1:
            refined_pseudogenes.append(pg)
            continue

        # Process multi-exon pseudogenes
        refined_exons = []
        current_merged_exon = None

        for i, exon in enumerate(exons):
            # If this is the first exon or a new merged exon needs to be started
            if current_merged_exon is None:
                current_merged_exon = exon.copy()

            # If this is the last exon, add the current merged exon
            if i == len(exons) - 1:
                if current_merged_exon != exon:  # If we're already merging
                    # Update merged exon end position
                    current_merged_exon["end"] = exon["end"]

                    # Update metrics for the merged exon
                    current_merged_exon["frameshifts"] += exon["frameshifts"]
                    current_merged_exon["insertions"] += exon["insertions"]
                    current_merged_exon["deletions"] += exon["deletions"]
                    current_merged_exon["stop_codons"] += exon["stop_codons"]
                    current_merged_exon["score"] += exon["score"]
                    current_merged_exon["identity"] = (
                        current_merged_exon["identity"] + exon["identity"]
                    ) / 2

                    merged_count += 1

                refined_exons.append(current_merged_exon)
                continue

            # Check for canonical splice sites between this exon and the next
            next_exon = exons[i + 1]
            has_canonical_splice = False

            # Extract genomic sequence around potential splice sites
            if strand == "+":
                # For positive strand:
                # - Donor site (GT) should be at the end of current exon
                # - Acceptor site (AG) should be at the start of next exon

                # Extract sequence around donor site (current exon end)
                donor_region_start = max(exon["start"], exon["end"] - SPLICE_CONTEXT)
                donor_region_end = min(
                    next_exon["start"], exon["end"] + SPLICE_CONTEXT + 2
                )
                donor_region = genome_seqs[chrom][donor_region_start:donor_region_end]

                # Extract sequence around acceptor site (next exon start)
                acceptor_region_start = max(
                    exon["end"], next_exon["start"] - SPLICE_CONTEXT - 2
                )
                acceptor_region_end = min(
                    next_exon["end"], next_exon["start"] + SPLICE_CONTEXT
                )
                acceptor_region = genome_seqs[chrom][
                    acceptor_region_start:acceptor_region_end
                ]

                # Find donor and acceptor positions (if they exist)
                donor_positions = _find_sites(donor_region, DONOR_SITE)
                acceptor_positions = _find_sites(acceptor_region, ACCEPTOR_SITE)

                # Determine the start and end positions of the intron
                intron_start = (
                    donor_region_start + donor_positions[-1]
                    if donor_positions
                    else None
                )
                intron_end = (
                    acceptor_region_start + acceptor_positions[0] + 1
                    if acceptor_positions
                    else None
                )

                # Check if there are both splice sites
                has_canonical_splice = (
                    intron_start is not None and intron_end is not None
                )

            else:  # strand == '-'
                # For negative strand:
                # - Donor site (GT becomes AC on reverse strand) should be at the start of next exon
                # - Acceptor site (AG becomes CT on reverse strand) should be at the end of current exon

                # Extract sequence around donor site (next exon start)
                donor_region_start = max(
                    exon["end"], next_exon["start"] - SPLICE_CONTEXT - 2
                )
                donor_region_end = min(
                    next_exon["end"], next_exon["start"] + SPLICE_CONTEXT
                )
                donor_region = genome_seqs[chrom][donor_region_start:donor_region_end]

                # Extract sequence around acceptor site (current exon end)
                acceptor_region_start = max(exon["start"], exon["end"] - SPLICE_CONTEXT)
                acceptor_region_end = min(
                    next_exon["start"], exon["end"] + SPLICE_CONTEXT + 2
                )
                acceptor_region = genome_seqs[chrom][
                    acceptor_region_start:acceptor_region_end
                ]

                # We look for AC and CT in the sequence (reverse complement of GT and AG)
                donor_positions = _find_sites(
                    donor_region, _reverse_complement(DONOR_SITE)
                )
                acceptor_positions = _find_sites(
                    acceptor_region, _reverse_complement(ACCEPTOR_SITE)
                )

                # Determine the start and end positions of the intron
                intron_start = (
                    acceptor_region_start + acceptor_positions[-1]
                    if acceptor_positions
                    else None
                )
                intron_end = (
                    donor_region_start + donor_positions[0] + 1
                    if donor_positions
                    else None
                )

                # Check if there are both splice sites
                has_canonical_splice = (
                    intron_start is not None and intron_end is not None
                )

            # If canonical splice sites found, add current merged exon and start a new one
            if has_canonical_splice:
                current_merged_exon["end"] = intron_start - 1
                next_exon["start"] = intron_end + 1
                refined_exons.append(current_merged_exon)
                current_merged_exon = next_exon.copy()
            else:
                # Merge with next exon
                current_merged_exon["end"] = next_exon["end"]

                # Update metrics for the merged exon
                current_merged_exon["frameshifts"] += next_exon["frameshifts"]
                current_merged_exon["insertions"] += next_exon["insertions"]
                current_merged_exon["deletions"] += next_exon["deletions"]
                current_merged_exon["stop_codons"] += next_exon["stop_codons"]
                current_merged_exon["score"] += next_exon["score"]
                current_merged_exon["identity"] = (
                    current_merged_exon["identity"] + next_exon["identity"]
                ) / 2

                merged_count += 1

        # Create a copy of the pseudogene with refined exons
        refined_pg = pg.copy()
        refined_pg["exons"] = refined_exons
        refined_pg["exon_count"] = len(refined_exons)

        # Update pseudogene statistics based on merged exons
        if refined_exons:
            refined_pg["start"] = min(exon["start"] for exon in refined_exons)
            refined_pg["end"] = max(exon["end"] for exon in refined_exons)

            # Recalculate stats for pseudogene based on new exon structure
            total_frameshifts = sum(exon["frameshifts"] for exon in refined_exons)
            total_insertions = sum(exon["insertions"] for exon in refined_exons)
            total_deletions = sum(exon["deletions"] for exon in refined_exons)
            total_stop_codons = sum(exon["stop_codons"] for exon in refined_exons)
            avg_identity = sum(exon["identity"] for exon in refined_exons) / len(
                refined_exons
            )

            refined_pg["frameshifts"] = total_frameshifts
            refined_pg["insertions"] = total_insertions
            refined_pg["deletions"] = total_deletions
            refined_pg["stop_codons"] = total_stop_codons
            refined_pg["identity"] = avg_identity

        refined_pseudogenes.append(refined_pg)

    logger.info(f"Exon refinement completed: {merged_count} exons merged")

    return refined_pseudogenes
