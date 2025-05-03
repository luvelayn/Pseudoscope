from collections import defaultdict
import csv

def _resolve_overlapping(pseudogenes_list):
    """Resolve overlapping pseudogenes from different protein or strand"""
    # Sort by chromosome and position
    pseudogenes_list.sort(key=lambda x: (x['chrom'], x['start']))
    
    # Identify overlapping pseudogenes
    non_overlapping = []
    i = 0
    while i < len(pseudogenes_list):
        current = pseudogenes_list[i]
        overlapping = [current]
        j = i + 1
        
        # Check if overlaps with previously found non-overlapping
        if (non_overlapping and 
        current['chrom'] == non_overlapping[-1]['chrom'] and 
        current['start'] <= non_overlapping[-1]['end']):
        
            if current['length'] > non_overlapping[-1]['length']:
                non_overlapping.pop()
            else:
                i += 1
                continue

        # Find all pseudogenes that overlap with current
        while (j < len(pseudogenes_list) and 
        pseudogenes_list[j]['chrom'] == current['chrom'] and
        pseudogenes_list[j]['start'] <= current['end']):
            overlapping.append(pseudogenes_list[j])
            j += 1
            
        # If none overlaps
        if len(overlapping) == 1:
            non_overlapping.append(current)
            i += 1
        else:
            # Select best pseudogene based on length
            best = max(overlapping, key=lambda x: x['length'])
            non_overlapping.append(best)
        
            # Skip all overlapping pseudogenes
            i = j
    
    return non_overlapping

def create_clusters(input_tsv, max_intron_length, logger):
    """Create pseudogene candidates from exons
    
    Returns:
        list: A list of lists of dictionaries, where each inner list represents a pseudogene
              and contains dictionaries with exon information (protein, chrom, start, end, strand, evalue)
    """
    # Read exons from TSV file
    exons = []
    
    try:
        with open(input_tsv, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Convert numeric fields
                exon = {
                    'protein': row['protein'],
                    'chrom': row['chrom'],
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'strand': row['strand'],
                    'evalue': float(row['evalue']),
                    'score': float(row['score'])
                }
                exons.append(exon)
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        return []
    
    # Group exons by protein, chromosome and strand direction
    grouped_exons = defaultdict(list)
    for exon in exons:
        key = (exon['protein'], exon['chrom'], exon['strand'])  # (protein, chromosome, strand)
        grouped_exons[key].append(exon)
    
    # Sort exons within each group by position
    for key in grouped_exons:
        grouped_exons[key].sort(key=lambda x: x['start'])
    
    pseudogenes = []
    
    # Process each group
    for key, group_exons in grouped_exons.items():
        protein, chrom, strand = key
        
        # Special case: if only one exon
        if len(group_exons) == 1:
            exon = group_exons[0]
            pseudogenes.append({
                'protein': protein,
                'chrom': chrom,
                'start': exon['start'],
                'end': exon['end'],
                'length': exon['end'] - exon['start'] + 1,
                'strand': strand,
                'evalue': exon['evalue'],
                'score': exon['score'], 
                'exons': [exon]  # Include the entire exon dictionary
            })
            continue
            
        # Process multiple exons
        i = 0
        while i < len(group_exons):
            # Start a new potential cluster
            current_cluster = [group_exons[i]]
            start = group_exons[i]['start']
            end = group_exons[i]['end']
            
            # Look for exons that can be added to this cluster
            j = i + 1
            while j < len(group_exons):
                next_exon = group_exons[j]
                
                # Check if exon belongs to the same pseudogene (within max_intron_length)
                if next_exon['start'] <= end + max_intron_length:
                    # Add to current cluster
                    current_cluster.append(next_exon)
                    
                    # Update cluster end boundary
                    end = max(end, next_exon['end'])
                    j += 1
                else:
                    # Next exon is too far away, stop adding to current cluster
                    break
            
            # Create pseudogene from the current group
            best_evalue = min([e['evalue'] for e in current_cluster])
            best_score = max([e['score'] for e in current_cluster])
            
            pseudogenes.append({
                'protein': protein,
                'chrom': chrom,
                'start': start,
                'end': end,
                'length': end - start + 1,
                'strand': strand,
                'evalue': best_evalue,
                'score': best_score, 
                'exons': current_cluster.copy()  # Include all exon dictionaries
            })
            
            # Move to the next unprocessed exon
            i = j

    logger.info(f"Exons merged into clusters, {len(pseudogenes)} clusters created")

    # Filtering out overlapping pseudogenes
    non_overlapping_pgs = _resolve_overlapping(pseudogenes)
    logger.info(f"Overlapping clusters filtered, {len(non_overlapping_pgs)} clusters retained")
    
    exon_clusters = []
    for pg in non_overlapping_pgs: 
        exon_clusters.append(pg['exons'])

    return exon_clusters