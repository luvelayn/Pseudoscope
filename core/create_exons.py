import csv
import os

def create_exons(input_tsv, out_dir, logger):
        """Merge overlapping hits and closely located hits into exons"""
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
        
        exons = []
        
        # Process each group
        for key, group_hits in grouped_hits.items():
            protein, chrom, strand = key
            
            # Processing if only one hit
            if len(group_hits) == 1:
                hit = group_hits[0]
                exons.append({
                    'protein': protein,
                    'chrom': chrom,
                    'start': hit['sstart'],
                    'end': hit['send'],
                    'strand': hit['strand'],
                    'evalue': hit['evalue']
                })
                continue
                
            # Check for overlapping or very close hits
            current_group = [group_hits[0]]
            for i in range(1, len(group_hits)):
                current_hit = group_hits[i]
                last_hit = current_group[-1]
                
                # Check if hits are the part of the same exon
                if (current_hit['sstart'] <= last_hit['send'] + 30):
                    current_group.append(current_hit)
                else:
                    # Create exon from the current group
                    if current_group:
                        best_evalue = min([h['evalue'] for h in current_group])
                        start = min([h['sstart'] for h in current_group])
                        end = max([h['send'] for h in current_group])
                        
                        exons.append({
                            'protein': protein,
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'evalue': best_evalue
                        })
                    
                    # Start a new group with the current hit
                    current_group = [current_hit]
            
            # Don't forget the last group
            if current_group:
                best_evalue = min([h['evalue'] for h in current_group])
                start = min([h['sstart'] for h in current_group])
                end = max([h['send'] for h in current_group])
                
                exons.append({
                    'protein': protein,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'evalue': best_evalue
                })

        logger.info(f"Hits merged into exons successfully. {len(exons)} exons created.")

        # Save exons to file
        exons_file = os.path.join(out_dir, "exons.tsv")

        with open(exons_file, 'w') as f:
            f.write("protein\tchrom\tstart\tend\tstrand\tevalue\n")
            for exon in exons:
                f.write(f"{exon['protein']}\t{exon['chrom']}\t{exon['start']}\t{exon['end']}\t"
                        f"{exon['strand']}\t{exon['evalue']}\n")
        
        return exons_file