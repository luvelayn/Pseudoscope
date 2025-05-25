from collections import defaultdict


def classify_pseudogenes(pseudogene_results, genome_seqs, logger):
    """
    Classify pseudogenes into different types based on their features.

    Categorizes pseudogenes as "retrotransposed" (processed pseudogenes derived from mRNA),
    "duplicated" (unprocessed pseudogenes from DNA duplication), or "fragment" (partial pseudogenes)
    based on exon structure, coverage, and the presence of poly-A tails.

    Parameters:
    -----------
    pseudogene_results : list
        List of dictionaries containing pseudogene information from run_tfasty.
    genome_seqs : dict
        Dictionary mapping chromosome names to their sequences.
    logger : logging.Logger
        Logger object for reporting.

    Returns:
    --------
    list
        The same list of pseudogenes with an additional "type" field in each entry.
        Types are: "retrotransposed", "duplicated", or "fragment".
    """
    # Set parameters
    coverage_threshold = 0.7
    polya_window = 20
    polya_min_count = 11
    polya_search_range = 500

    # Sort pseudogenes by chromosome and position for detecting adjacent pseudogenes
    sorted_pseudogenes = sorted(
        pseudogene_results, key=lambda pg: (pg["chrom"], pg["start"])
    )

    # Create a lookup dictionary for pseudogenes by chromosome
    pg_by_chrom = defaultdict(list)
    for pg in sorted_pseudogenes:
        pg_by_chrom[pg["chrom"]].append(pg)

    # Track classification statistics
    stats = {"retrotransposed": 0, "duplicated": 0, "fragment": 0}

    # Process each pseudogene
    for i, pg in enumerate(pseudogene_results):
        # Check if coverage is below threshold - classify as fragment
        if pg["coverage"] < coverage_threshold:
            pg["type"] = "fragment"
            stats["fragment"] += 1

        # Check if single exon
        elif pg["exon_count"] == 1:
            # Check for poly-A tail
            chrom = pg["chrom"]

            # Determine search boundaries based on strand
            if pg["strand"] == "+":
                # For positive strand, search downstream from the end position
                search_start = pg["end"]

                # Find the next pseudogene on this chromosome to limit our search
                search_end = search_start + polya_search_range
                for other_pg in pg_by_chrom.get(chrom, []):
                    if (
                        other_pg["start"] > search_start
                        and other_pg["start"] < search_end
                    ):
                        search_end = other_pg["start"]

                # Don't go beyond chromosome end
                search_end = min(search_end, len(genome_seqs[chrom]))

                # Get the sequence to search in
                search_seq = genome_seqs[chrom][search_start:search_end].lower()

                # Look for a 20-base window with at least 15 'A's
                has_polya = False
                for i in range(len(search_seq) - polya_window + 1):
                    window = search_seq[i : i + polya_window]
                    if window.count("a") >= polya_min_count:
                        has_polya = True
                        break

            else:  # Negative strand
                # For negative strand, search upstream from the start position
                search_end = pg["start"]

                # Find the previous pseudogene on this chromosome to limit our search
                search_start = max(0, search_end - polya_search_range)
                for other_pg in pg_by_chrom.get(chrom, []):
                    if other_pg["end"] < search_end and other_pg["end"] > search_start:
                        search_start = other_pg["end"]

                # Get the sequence to search in
                search_seq = genome_seqs[chrom][search_start:search_end].lower()

                # Look for a 20-base window with at least 15 'T's (reverse complement of poly-A)
                has_polya = False
                for i in range(len(search_seq) - polya_window + 1):
                    window = search_seq[i : i + polya_window]
                    if window.count("t") >= polya_min_count:
                        has_polya = True
                        break

            # Classify based on poly-A presence
            if has_polya:
                pg["type"] = "retrotransposed"
                stats["retrotransposed"] += 1
            else:
                pg["type"] = "duplicated"
                stats["duplicated"] += 1
        else:
            # Multiple exons - classify as duplicated
            pg["type"] = "duplicated"
            stats["duplicated"] += 1

    # Log classification statistics
    logger.info(f"Pseudogene classification complete:")
    logger.info(f"  Retrotransposed pseudogenes: {stats['retrotransposed']}")
    logger.info(f"  Duplicated pseudogenes: {stats['duplicated']}")
    logger.info(f"  Pseudogene fragments: {stats['fragment']}")

    # Return the updated pseudogene results
    return pseudogene_results
