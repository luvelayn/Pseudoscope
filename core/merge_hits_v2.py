import csv
import os
from Bio import SeqIO
from Bio import pairwise2

def calculate_protein_similarity(seq1, seq2):
    """Calculate sequence identity between two protein sequences"""
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    matches = sum(a == b for a, b in zip(best_alignment.seqA, best_alignment.seqB))
    identity = (matches / len(best_alignment.seqA)) * 100
    return identity

def merge_hits(input_tsv, protein_file, max_intron_length, out_dir, logger):
        """Merge overlapping hits and closely located hits into pseudogene candidates"""
        # Load protein sequences and lengths
        protein_sequences = {}
        protein_lengths = {}
        for record in SeqIO.parse(protein_file, "fasta"):
            protein_sequences[record.id] = str(record.seq)
            protein_lengths[record.id] = len(record.seq)

        # Read hits from TSV file
        hits = []
        try:
            with open(input_tsv, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    # Convert numeric fields
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
            logger.error(f"Error reading input file {input_tsv}: {e}")
            return []
        
        # Group hits by protein, chromosome and strand direction
        grouped_hits = {}
        for hit in hits:
            key = (hit['qseqid'], hit['sseqid'], hit['strand'])  # (protein, chromosome, strand)
            if key not in grouped_hits:
                grouped_hits[key] = []
            grouped_hits[key].append(hit)

        # Sort hits within each group by position
        for key in grouped_hits:
            grouped_hits[key].sort(key=lambda x: x['sstart'])

        pseudogene_candidates = []

        # Function to check if two hits should be in the same cluster 
        # Suggesting prev_hit start is before curr_hit start
        def can_be_merged(prev_hit, curr_hit):
                # Check genomic proximity
                if curr_hit['sstart'] > prev_hit['send'] + max_intron_length:
                    return False
                    
                # Check query positions based on strand
                if prev_hit['strand'] == '+':
                    # For forward strand, curr_hit's qstart should be reasonably close to prev_hit's qend
                    # Allow for some overlap or gap, but maintain overall order
                    return curr_hit['qstart'] > prev_hit['qstart'] and abs(curr_hit['qstart'] - prev_hit['qend']) < 50
                else:  # '-' strand
                    # For reverse strand, curr_hit's qend should be reasonably close to prev_hit's qstart
                    # Remember query coordinates are from the protein perspective
                    return curr_hit['qend'] < prev_hit['qend'] and abs(prev_hit['qstart'] - curr_hit['qend']) < 50
                
        # Process each group
        for key, group_hits in grouped_hits.items():
            protein, chrom, strand = key
            
            # Processing if only one hit
            if len(group_hits) == 1:
                hit = group_hits[0]
                pseudogene_candidates.append({
                    'protein': protein,
                    'chrom': chrom,
                    'start': hit['sstart'],
                    'end': hit['send'],
                    'length': hit['send'] - hit['sstart'] + 1,
                    'strand': hit['strand'],
                    'evalue': hit['evalue'],
                    'coverage': (hit['qend'] - hit['qstart']) / protein_lengths[protein]
                })
                continue
            
            # Create all possible clusters by considering each hit as a starting point
            all_clusters = []
            
            for start_idx in range(len(group_hits)):
                current_cluster = [group_hits[start_idx]]
                
                # Try to extend cluster with subsequent hits
                for i in range(start_idx + 1, len(group_hits)):
                    curr_hit = group_hits[i]
                    prev_hit = current_cluster[-1]

                    if can_be_merged(prev_hit, curr_hit):
                        current_cluster.append(curr_hit)
                
                all_clusters.append(current_cluster)
            
            # Convert clusters to pseudogene candidates
            for cluster in all_clusters:
                best_evalue = min([h['evalue'] for h in cluster])
                start = min([h['sstart'] for h in cluster])
                end = max([h['send'] for h in cluster])
                qstart = min([h['qstart'] for h in cluster])
                qend = max([h['qend'] for h in cluster])
                    
                coverage = (qend - qstart + 1) / protein_lengths[protein]
                
                pseudogene_candidates.append({
                    'protein': protein,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'length': end - start + 1,
                    'strand': strand,
                    'evalue': best_evalue,
                    'coverage': coverage
                })

        # Function to check if two pseudogenes overlap
        def overlaps(pg1, pg2):
            # Only consider overlap if they're on the same chromosome and protein
            if (pg1['protein'] != pg2['protein']) or (pg1['chrom'] != pg2['chrom']) or (pg1['strand'] != pg2['strand']):
                return False
            
            # Check for genomic overlap
            return max(pg1['start'], pg2['start']) <= min(pg1['end'], pg2['end'])

        # Filter overlapping pseudogenes to keep the longest
        # Sort candidates by length (descending)
        pseudogene_candidates.sort(key=lambda x: x['length'], reverse=True)

        merged_pseudogenes = []
        for candidate in pseudogene_candidates:
            # Check if this candidate overlaps with any existing accepted pseudogene
            if not any(overlaps(candidate, pg) for pg in merged_pseudogenes):
                merged_pseudogenes.append(candidate)

        logger.info(f"Created {len(merged_pseudogenes)} pseudogene candidates after merging same protein hits")

        # Resolve overlapping pseudogenes from different proteins
        merged_pseudogenes.sort(key=lambda x: (x['chrom'], x['strand'], x['start']))
        
        # Identify overlapping pseudogenes
        non_overlapping = []
        i = 0
        while i < len(merged_pseudogenes):
            current = merged_pseudogenes[i]
            overlapping = [current]
            j = i + 1
            
            # Check if overlaps with previously found non-overlapping
            if (non_overlapping and 
            current['chrom'] == non_overlapping[-1]['chrom'] and
            current['strand'] == non_overlapping[-1]['strand'] and 
            current['start'] <= non_overlapping[-1]['end']):
            
                if current['length'] > non_overlapping[-1]['length']:
                    non_overlapping.pop()
                else:
                    i += 1
                    continue

            # Find all pseudogenes that overlap with current
            while (j < len(merged_pseudogenes) and 
            merged_pseudogenes[j]['chrom'] == current['chrom'] and 
            merged_pseudogenes[j]['strand'] == current['strand'] and
            merged_pseudogenes[j]['start'] <= current['end']):
                overlapping.append(merged_pseudogenes[j])
                j += 1
                
            # If none overlaps
            if  len(overlapping) == 1:
                non_overlapping.append(current)
                i += 1
            else:
                # Select best pseudogene based on lentgh
                best = max(overlapping, key=lambda x: x['length'])
                non_overlapping.append(best)
            
                # Skip all overlapping pseudogenes
                i = j
        
        logger.info(f"Created {len(non_overlapping)} pseudogene candidates after merging overlapping different protein hits")
        
        # Merge close alignments (<2500) from different proteins
        sorted_alignments = sorted(non_overlapping, key=lambda x: (x['chrom'], x['strand'], x['start']))
        merged_alignments = []

        # Cached similarity calculations to avoid repeating
        similarity_cache = {}
        
        i = 0
        while i < len(sorted_alignments):
            current = sorted_alignments[i]
            merged_group = [current]
            j = i + 1
            
            # Find close alignments (within 2500 bases)
            while (j < len(sorted_alignments) and 
                   sorted_alignments[j]['chrom'] == current['chrom'] and
                   sorted_alignments[j]['strand'] == current['strand'] and
                   sorted_alignments[j]['start'] <= current['end'] + 2500):
                
                next_alignment = sorted_alignments[j]

                # Compute protein similarity if not already cached
                pair_key = tuple(sorted([current['protein'], next_alignment['protein']]))
                if pair_key not in similarity_cache:
                    seq1 = protein_sequences[current['protein']]
                    seq2 = protein_sequences[next_alignment['protein']]
                    similarity_cache[pair_key] = calculate_protein_similarity(seq1, seq2)
                
                # Check if proteins are similar enough
                similarity_threshold = 30
                if similarity_cache[pair_key] >= similarity_threshold:
                    merged_group.append(next_alignment)
                    
                j += 1
            
            if len(merged_group) == 1:
                # No merging needed
                merged_alignments.append(current)
                i += 1
            else:
                # Create a merged alignment
                start = min([a['start'] for a in merged_group])
                end = max([a['end'] for a in merged_group])
                # Use the protein with the best coverage
                best_protein = max(merged_group, key=lambda x: x['coverage'])['protein']
                # Sum coverage, but cap at 1.0
                total_coverage = min(1.0, sum([a['coverage'] for a in merged_group]))
                best_evalue = min([a['evalue'] for a in merged_group])
                
                merged_alignments.append({
                    'protein': best_protein,
                    'chrom': current['chrom'],
                    'start': start,
                    'end': end,
                    'length': end - start + 1,
                    'strand': current['strand'],
                    'evalue': best_evalue,
                    'coverage': total_coverage
                })
                
                # Skip all merged alignments
                i = j
        
        # # Filter out alignments with coverage < 0.2
        # filtered_alignments = [a for a in merged_alignments if a['coverage'] >= 0.2]
        
        logger.info(f"Final result: {len(merged_alignments)} pseudogene candidates after merging close different protein hits")
        
        # Save merged hits to file
        merged_hits_file = os.path.join(out_dir, "merged_hits.tsv")

        with open(merged_hits_file, 'w') as f:
            f.write("protein\tchrom\tstart\tend\tstrand\tevalue\tcoverage\n")
            for pg in merged_alignments:
                f.write(f"{pg['protein']}\t{pg['chrom']}\t{pg['start']}\t{pg['end']}\t"
                        f"{pg['strand']}\t{pg['evalue']}\t{pg['coverage']}\n")
        
        return merged_hits_file