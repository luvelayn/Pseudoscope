import csv
import os
from collections import defaultdict

def create_exons(input_tsv, out_dir, logger):
    """
    Merge overlapping hits and closely located hits into exons.
    
    Parameters:
    -----------
    input_tsv : str
        Path to the input TSV file containing filtered BLAST hits
    out_dir : str
        Directory to store output files
    logger : logging.Logger
        Logger object for reporting
        
    Returns:
    --------
    str
        Path to the output exons file
    """
    # Read hits from TSV file
    hits = []
    try:
        with open(input_tsv, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Convert numeric fields with proper error handling
                hit = {
                    'qseqid': row['qseqid'],
                    'sseqid': row['sseqid'],
                    'pident': float(row['pident']),
                    'length': int(row['length']),
                    'qstart': int(row['qstart']),
                    'qend': int(row['qend']),
                    'sstart': int(row['sstart']),
                    'send': int(row['send']),
                    'strand': row['strand'],
                    'evalue': float(row['evalue'])
                }
                hits.append(hit)
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        return None
        
    # Group hits by protein, chromosome and strand direction
    grouped_hits = defaultdict(list)
    for hit in hits:
        key = (hit['qseqid'], hit['sseqid'], hit['strand'])  # (protein, chromosome, strand)
        grouped_hits[key].append(hit)
    
    # Sort hits within each group by position
    for key in grouped_hits:
        grouped_hits[key].sort(key=lambda x: x['sstart'])
    
    exons = []
    
    # Process each group
    for key, group_hits in grouped_hits.items():
        protein, chrom, strand = key
        
        # Special case: if only one hit
        if len(group_hits) == 1:
            hit = group_hits[0]
            exons.append({
                'protein': protein,
                'chrom': chrom,
                'start': hit['sstart'],
                'end': hit['send'],
                'qstart': hit['qstart'],
                'qend': hit['qend'],
                'strand': hit['strand'],
                'evalue': hit['evalue'],
                'score': hit['length'] * hit['pident'],  # Add score metric for better ranking
            })
            continue
            
        # Process multiple hits: merge overlapping or very close hits
        i = 0
        while i < len(group_hits):
            # Start a new potential exon
            current_merge = [group_hits[i]]
            best_evalue = group_hits[i]['evalue']
            best_score = group_hits[i]['length'] * group_hits[i]['pident']
            start = group_hits[i]['sstart']
            end = group_hits[i]['send']
            qstart = group_hits[i]['qstart']
            qend = group_hits[i]['qend']
            
            # Look for hits that can be merged into this exon
            j = i + 1
            while j < len(group_hits):
                next_hit = group_hits[j]
                
                # Check if hit overlaps or is very close (within 30bp)
                if next_hit['sstart'] <= end + 30:
                    # Merge this hit
                    current_merge.append(next_hit)
                    
                    # Update best metrics
                    if next_hit['evalue'] < best_evalue:
                        best_evalue = next_hit['evalue']
                    
                    hit_score = next_hit['length'] * next_hit['pident']
                    if hit_score > best_score:
                        best_score = hit_score
                    
                    # Expand exon boundaries
                    end = max(end, next_hit['send'])
                    qstart = min(qstart, next_hit['qstart'])
                    qend = max(qend, next_hit['qend'])
                    j += 1
                else:
                    # Next hit is too far away, stop merging
                    break
            
            # Create an exon from the current merged hits
            exons.append({
                'protein': protein,
                'chrom': chrom,
                'start': start,
                'end': end,
                'qstart': qstart,
                'qend': qend,
                'strand': strand,
                'evalue': best_evalue,
                'score': best_score
            })
            
            # Move to the next unprocessed hit
            i = j

    logger.info(f"Hits merged into exons: {len(exons)} exons created")

    # Save exons to file
    try:
        exons_file = os.path.join(out_dir, "exons.tsv")
        
        with open(exons_file, 'w') as f:
            # Add more fields to the output
            f.write("protein\tchrom\tstart\tend\tqstart\tqend\tstrand\tevalue\tscore\n")
            for exon in exons:
                f.write(f"{exon['protein']}\t{exon['chrom']}\t{exon['start']}\t{exon['end']}\t"
                        f"{exon['qstart']}\t{exon['qend']}\t{exon['strand']}\t"
                        f"{exon['evalue']}\t{exon['score']}\n")
        
        return exons_file
    except Exception as e:
        logger.error(f"Error writing exons file: {e}")
        return None