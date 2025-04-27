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
        
        # Resolve overlapping and close hits from same proteins    
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
        
        merged_pseudogenes = []
        
        # Process each group
        for key, group_hits in grouped_hits.items():
            protein, chrom, strand = key
            
            # Processing if only one hit
            if len(group_hits) == 1:
                hit = group_hits[0]
                merged_pseudogenes.append({
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
                
            # Check for overlapping or nearby hits
            current_group = [group_hits[0]]
            for i in range(1, len(group_hits)):
                current_hit = group_hits[i]
                last_hit = current_group[-1]
                
                # Check if hits are the part of the same pseudogene
                if (current_hit['sstart'] <= last_hit['send'] + max_intron_length):
                    current_group.append(current_hit)
                else:
                    # Create a pseudogene from the current group
                    if current_group:
                        best_evalue = min([h['evalue'] for h in current_group])
                        start = min([h['sstart'] for h in current_group])
                        end = max([h['send'] for h in current_group])
                        coverage = (max([h['qend'] for h in current_group]) - min([h['qstart'] for h in current_group])) / protein_lengths[protein] 
                        
                        merged_pseudogenes.append({
                            'protein': protein,
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'length': end - start + 1,
                            'strand': strand,
                            'evalue': best_evalue,
                            'coverage': coverage
                        })
                    
                    # Start a new group with the current hit
                    current_group = [current_hit]
            
            # Don't forget the last group
            if current_group:
                best_evalue = min([h['evalue'] for h in current_group])
                start = min([h['sstart'] for h in current_group])
                end = max([h['send'] for h in current_group])
                coverage = (max([h['qend'] for h in current_group]) - min([h['qstart'] for h in current_group])) / protein_lengths[protein] 
                
                merged_pseudogenes.append({
                    'protein': protein,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'length': end - start + 1,
                    'strand': strand,
                    'evalue': best_evalue,
                    'coverage': coverage
                })

        logger.info(f"Created {len(merged_pseudogenes)} pseudogene candidates after merging same protein hits")

        # Merge overlapping pseudogenes from different but similar proteins
        # Sort by chromosome, strand and position
        merged_pseudogenes.sort(key=lambda x: (x['chrom'], x['strand'], x['start']))
        
        # Create a cache for protein similarity calculations
        similarity_cache = {}
        
        # Define similarity threshold
        similarity_threshold = 50  # Adjust this value as needed

        # Define maximum overlap percentage for comparing pseudogenes
        max_overlap_percentage = 20  # Only compare pseudogenes with less than 20% overlap
        
        # Initialize the list for pseudogenes after similar protein merging
        similar_protein_merged = []
        i = 0
        
        while i < len(merged_pseudogenes):
            current = merged_pseudogenes[i]
            merged_group = [current]
            j = i + 1
            
            # Find pseudogenes that might qualify for merging (they overlap but not too much)
            while (j < len(merged_pseudogenes) and 
                   merged_pseudogenes[j]['chrom'] == current['chrom'] and
                   merged_pseudogenes[j]['strand'] == current['strand'] and
                   # Check if they overlap
                   merged_pseudogenes[j]['start'] <= current['end']):
                
                next_pg = merged_pseudogenes[j]

                # Calculate overlap
                overlap_start = max(current['start'], next_pg['start'])
                overlap_end = min(current['end'], next_pg['end'])
                overlap_length = max(0, overlap_end - overlap_start + 1)
                
                # Calculate overlap percentage relative to both pseudogenes
                current_overlap_percent = (overlap_length / current['length']) * 100
                next_overlap_percent = (overlap_length / next_pg['length']) * 100
                
                # Skip if overlap is too large
                if current_overlap_percent > max_overlap_percentage or next_overlap_percent > max_overlap_percentage:
                    j += 1
                    continue
                
                # Calculate protein similarity if not already in cache
                protein_pair = tuple(sorted([current['protein'], next_pg['protein']]))
                if protein_pair not in similarity_cache:
                    seq1 = protein_sequences[current['protein']]
                    seq2 = protein_sequences[next_pg['protein']]
                    similarity_cache[protein_pair] = calculate_protein_similarity(seq1, seq2)
                
                # If proteins are similar enough, merge the pseudogenes
                if similarity_cache[protein_pair] >= similarity_threshold:
                    merged_group.append(next_pg)
                
                j += 1
            
            if len(merged_group) == 1:
                # No merging needed
                similar_protein_merged.append(current)
                i += 1
            else:
                # Create a merged pseudogene
                start = min([pg['start'] for pg in merged_group])
                end = max([pg['end'] for pg in merged_group])
                
                # Choose the pseudogene with best E-value as the representative
                best_pg = min(merged_group, key=lambda x: x['evalue'])
                
                # Combine coverage, but cap at 1.0
                total_coverage = min(1.0, sum([pg['coverage'] for pg in merged_group]))
                
                similar_protein_merged.append({
                    'protein': best_pg['protein'],  # Use protein from best e-value pseudogene
                    'chrom': current['chrom'],
                    'start': start,
                    'end': end,
                    'length': end - start + 1,
                    'strand': current['strand'],
                    'evalue': best_pg['evalue'],
                    'coverage': total_coverage
                })
                
                # Skip all merged pseudogenes
                i = j
        
        logger.info(f"Created {len(similar_protein_merged)} pseudogene candidates after merging pseudogenes from similar proteins")

        # Resolve remaining overlapping pseudogenes from different protein or strand
        similar_protein_merged.sort(key=lambda x: (x['chrom'], x['start']))
        
        # Identify overlapping pseudogenes
        non_overlapping = []
        i = 0
        while i < len(similar_protein_merged):
            current = similar_protein_merged[i]
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
            while (j < len(similar_protein_merged) and 
            similar_protein_merged[j]['chrom'] == current['chrom'] and
            similar_protein_merged[j]['start'] <= current['end']):
                overlapping.append(similar_protein_merged[j])
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
        
        logger.info(f"Created {len(non_overlapping)} pseudogene candidates after filtering different protein/strand hits")
        
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
                similarity_threshold = 40
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
        
        logger.info(f"Final result: {len(merged_alignments)} pseudogene candidates after merging close alignments from different proteins")
        
        # Save merged hits to file
        merged_hits_file = os.path.join(out_dir, "merged_hits.tsv")

        with open(merged_hits_file, 'w') as f:
            f.write("protein\tchrom\tstart\tend\tstrand\tevalue\tcoverage\n")
            for pg in merged_alignments:
                f.write(f"{pg['protein']}\t{pg['chrom']}\t{pg['start']}\t{pg['end']}\t"
                        f"{pg['strand']}\t{pg['evalue']}\t{pg['coverage']}\n")
        
        return merged_hits_file