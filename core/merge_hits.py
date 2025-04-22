import csv
import os


def merge_hits(input_tsv, max_intron_length, out_dir, logger):
        """Merge overlapping hits and closely located hits into pseudogene candidates"""
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
                    'strand': hit['strand'],
                    'evalue': hit['evalue']
                })
                continue
                
            # Check for overlapping or nearby hits
            current_group = [group_hits[0]]
            for i in range(1, len(group_hits)):
                current_hit = group_hits[i]
                last_hit = current_group[-1]
                
                # Check if hits are close enough to be part of the same pseudogene
                # We consider hits within max_intron_length distance
                # TODO подумать как можно мержить более точно
                if current_hit['sstart'] <= last_hit['send'] + max_intron_length:
                    current_group.append(current_hit)
                else:
                    # Create a pseudogene from the current group
                    if current_group:
                        best_evalue = min([h['evalue'] for h in current_group])
                        
                        merged_pseudogenes.append({
                            'protein': protein,
                            'chrom': chrom,
                            'start': min([h['sstart'] for h in current_group]),
                            'end': max([h['send'] for h in current_group]),
                            'strand': strand,
                            'evalue': best_evalue
                        })
                    
                    # Start a new group with the current hit
                    current_group = [current_hit]
            
            # Don't forget the last group
            if current_group:
                best_evalue = min([h['evalue'] for h in current_group])
                
                merged_pseudogenes.append({
                    'protein': protein,
                    'chrom': chrom,
                    'start': min([h['sstart'] for h in current_group]),
                    'end': max([h['send'] for h in current_group]),
                    'strand': strand,
                    'evalue': best_evalue
                })
        
        # Resolve overlapping pseudogenes from different proteins
        merged_pseudogenes.sort(key=lambda x: (x['chrom'], x['start']))
        
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
            current['start'] <= non_overlapping[-1]['end']):
            
                if current['evalue'] < non_overlapping[-1]['evalue']:
                    non_overlapping.pop()
                else:
                    i += 1
                    continue

            # Find all pseudogenes that overlap with current
            while (j < len(merged_pseudogenes) and 
            merged_pseudogenes[j]['chrom'] == current['chrom'] and 
            merged_pseudogenes[j]['start'] <= current['end']):
                overlapping.append(merged_pseudogenes[j])
                j += 1
                
            # If none overlaps
            if  len(overlapping) == 1:
                non_overlapping.append(current)
                i += 1
            else:
                # Select best pseudogene based on e-value
                # TODO можно сделать поточнее
                best = min(overlapping, key=lambda x: x['evalue'])
                non_overlapping.append(best)
            
                # Skip all overlapping pseudogenes
                i = j
        
        logger.info(f"Created {len(non_overlapping)} pseudogene candidates after merging")
        
        # Save merged hits to file
        merged_hits_file = os.path.join(out_dir, "merged_hits.tsv")

        with open(merged_hits_file, 'w') as f:
            f.write("protein\tchrom\tstart\tend\tstrand\tevalue\n")
            for pg in non_overlapping:
                f.write(f"{pg['protein']}\t{pg['chrom']}\t{pg['start']}\t{pg['end']}\t"
                        f"{pg['strand']}\t{pg['evalue']}\n")
        
        return non_overlapping
