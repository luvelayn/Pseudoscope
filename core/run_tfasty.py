import os
import subprocess
import tempfile

def run_tfasty(exon_clusters, protein_seqs, genome_seqs, evalue, coverage, identity, out_dir, logger):
    """Run precise re-alignment with tfasty for each exon separately, then merge into pseudogenes
    
    Parameters:
    -----------
    exon_clusters : list of lists
        Each inner list contains dictionaries representing exons of a pseudogene
        [
            [
                {'protein': 'protA', 'chrom': 'chr1', 'start': 1000, 'end': 1500, 'strand': '-', 'evalue': 9.26e-17},
                {'protein': 'protA', 'chrom': 'chr1', 'start': 1700, 'end': 2000', 'strand': '-', 'evalue': 9.26e-17}
            ],
            [
                {'protein': 'protB', 'chrom': 'chr2', 'start': 5000, 'end': 5200', 'strand': '+', 'evalue': 9.26e-17}
            ],
            ...
        ]
    protein_file : str
        Path to the protein fasta file
    genome_file : str
        Path to the genome fasta file
    out_dir : str
        Directory to store output files
    logger : logging.Logger
        Logger object for reporting
    """
    if not exon_clusters:
        logger.warning("No pseudogene candidates for re-alignment")
        return
    
    # Create directories for output
    tfasty_output_dir = os.path.join(out_dir, "tfasty_alignments")
    os.makedirs(tfasty_output_dir, exist_ok=True)
    
    # List to store all pseudogene results
    pseudogene_results = []
    
    # Process each cluster separately
    for cluster_idx, exon_cluster in enumerate(exon_clusters):
        if not exon_cluster:
            continue
        
        first_exon = exon_cluster[0]
        chrom = first_exon['chrom']
        protein = first_exon['protein']
        strand = first_exon['strand']
        
        # Check if protein and chromosome exist
        if protein not in protein_seqs:
            logger.warning(f"Protein {protein} not found in protein file. Skipping cluster {cluster_idx+1}.")
            continue
            
        if chrom not in genome_seqs:
            logger.warning(f"Chromosome {chrom} not found in genome file. Skipping cluster {cluster_idx+1}.")
            continue
        
        # Generate unique ID for this pseudogene
        pseudogene_id = f"PG{cluster_idx+1}"
        
        # Process each exon in the cluster
        exon_results = []
        
        for exon_idx, exon in enumerate(exon_cluster):
            try:
                # Add flanking sequence (30 nt on each side for tfasty)
                flank = 30
                exon_start = max(0, exon['start'] - flank)
                exon_end = min(len(genome_seqs[chrom]), exon['end'] + flank)
                
                # Extract sequence with flanking regions
                exon_seq = genome_seqs[chrom][exon_start:exon_end]
                
                # Create temporary files for exon and protein sequences
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as exon_file, \
                     tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as prot_file:
                    
                    exon_id = f"{pseudogene_id}_exon{exon_idx+1}"
                    exon_file.write(f">{exon_id}\n{exon_seq}\n")
                    prot_file.write(f">{protein}\n{protein_seqs[protein]}\n")
                    
                    exon_file_path = exon_file.name
                    prot_file_path = prot_file.name
                
                # Output file for tfasty results
                alignment_out = os.path.join(tfasty_output_dir, f"{exon_id}_alignment.txt")
                
                # Run tfasty
                tfasty_cmd = f"tfasty36 -m 10 {prot_file_path} {exon_file_path} > {alignment_out}"
                
                logger.debug(f"Running command: {tfasty_cmd}")
                subprocess.run(tfasty_cmd, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                
                # Parse tfasty output
                exon_info = {
                    'id': exon_id,
                    'chrom': chrom,
                    'start': exon_start,
                    'end': exon_end,
                    'protein': protein,
                    'protein_seq': protein_seqs[protein],
                    'strand': strand
                }
                
                parsed_alignment = _parse_tfasty_output(alignment_out, exon_info, logger)
                
                if parsed_alignment:
                    # Add the exon alignment to our results
                    exon_results.append({
                        'id': exon_id,
                        'alignment': parsed_alignment,
                        'offset': exon_start  # Store for coordinate conversion
                    })
                
                # Clean up temporary files
                try:
                    os.unlink(exon_file_path)
                    os.unlink(prot_file_path)
                except:
                    pass
                    
            except Exception as e:
                logger.error(f"Error processing exon {exon_idx+1} in cluster {cluster_idx+1}: {e}")
                import traceback
                logger.error(traceback.format_exc())
        
        # Skip if no exons were successfully aligned
        if not exon_results:
            logger.warning(f"No successful alignments for cluster {cluster_idx+1}. Skipping.")
            continue
        
        # Now merge exon alignments into a pseudogene
        merged_pseudogene = _merge_exons(pseudogene_id, exon_results, chrom, protein, strand, protein_seqs)
        
        # Filter pseudogenes based on evalue, coverage and identity
        if (merged_pseudogene['evalue'] <= evalue and 
            merged_pseudogene['coverage'] >= coverage and
            merged_pseudogene['identity'] >= identity):
            # Add to our results
            pseudogene_results.append({
                'id': pseudogene_id,
                'chrom': chrom,
                'protein': protein,
                'strand': strand,
                'start': merged_pseudogene['start'],
                'end': merged_pseudogene['end'],
                'frameshifts': merged_pseudogene['frameshifts'],
                'insertions': merged_pseudogene['insertions'],
                'deletions': merged_pseudogene['deletions'],
                'stop_codons': merged_pseudogene['stop_codons'],
                'exons': merged_pseudogene['exons'],
                'exon_count': len(merged_pseudogene['exons']),
                'score': merged_pseudogene['score'],
                'evalue': merged_pseudogene['evalue'],
                'identity': merged_pseudogene['identity'],
                'coverage': merged_pseudogene['coverage']
            })
    
    logger.info(f"Completed tfasty re-alignment: {len(pseudogene_results)} pseudogenes remained")
    
    return pseudogene_results

def _parse_tfasty_output(tfasty_output, exon_info, logger):
    """Parse tfasty output (format 10) and extract alignment information"""
    alignments = []
    
    try:
        # Read the tfasty output file
        with open(tfasty_output, 'r') as f:
            content = f.readlines()
        
        # Extract ids
        exon_id = exon_info['id']
        protein_id = exon_info['protein']
        
        # Parse the file
        i = 0
        while i < len(content) and not content[i].startswith(">>><<<"):
            line = content[i].strip()
            
            # Look for the start of an alignment block
            if line.startswith(f">>{exon_id}") or line.startswith(">--"):
                # Initialize alignment data structure
                alignment = {
                    "name": exon_id,
                    "strand": exon_info['strand'],  # Default to exon strand
                    "score": 0,
                    "evalue": 0,
                    "identity": 0,
                    "overlap": 0,
                    "protein_range": {},
                    "dna_range": {},
                    "protein_seq": "",
                    "dna_seq": "",
                    "frameshifts": 0,
                    "insertions": 0,
                    "deletions": 0,
                    "stop_codons": 0
                }
                
                # Parse alignment details
                i += 1
                while i < len(content) and not content[i].startswith(f">>{exon_id}") and not content[i].startswith(">--"):
                    line = content[i].strip()
                    
                    # Strand information (override if specified)
                    if line.startswith("; tfy_frame:"):
                        frame = line.split(":", 1)[1].strip()
                        alignment["strand"] = "+" if frame == "f" else "-"
                        i += 1

                    # Score and evalue information
                    elif line.startswith("; tfy_opt:"):
                        alignment["score"] = float(line.split(":", 1)[1].strip())
                        i += 1
                    elif line.startswith("; tfy_expect:"):
                        alignment["evalue"] = float(line.split(":", 1)[1].strip())
                        i += 1
                    
                    # Identity and overlap information
                    elif line.startswith("; sw_ident:"):
                        alignment["identity"] = float(line.split(":", 1)[1].strip())
                        i += 1
                    elif line.startswith("; sw_overlap:"):
                        alignment["overlap"] = int(line.split(":", 1)[1].strip())
                        i += 1
                    
                    # Get to the protein section
                    elif line.startswith(f">{protein_id}"):
                        # Skip to the protein coordinates
                        i += 1
                        while i < len(content) and content[i].startswith(";"):
                            line = content[i].strip()
                            if line.startswith("; al_start:"):
                                alignment["protein_range"]["start"] = int(line.split(":", 1)[1].strip())
                            elif line.startswith("; al_stop:"):
                                alignment["protein_range"]["end"] = int(line.split(":", 1)[1].strip())
                            i += 1
                        
                        # Collect the protein sequence (may span multiple lines)
                        protein_seq = ""
                        while i < len(content) and not content[i].startswith(">"):
                            protein_seq += content[i].strip()
                            i += 1
                        alignment["protein_seq"] = protein_seq
                        
                        # Get to the dna section (skip > line)
                        if i < len(content):
                            i += 1

                        # Get dna coordinates
                        while i < len(content) and content[i].startswith(";"):
                            line = content[i].strip()
                            if line.startswith("; al_start:"):
                                alignment["dna_range"]["start"] = int(line.split(":", 1)[1].strip())
                            elif line.startswith("; al_stop:"):
                                alignment["dna_range"]["end"] = int(line.split(":", 1)[1].strip())
                            i += 1
                        
                        # Collect the DNA sequence (may span multiple lines)
                        dna_seq = ""
                        while i < len(content) and not content[i].startswith(";") and not content[i].startswith(">"):
                            dna_seq += content[i].strip()
                            i += 1
                        alignment["dna_seq"] = dna_seq

                        # Skip to the next alignment block
                        while i < len(content) and not content[i].startswith(">"):
                            i += 1
                        
                        if i >= len(content) or content[i].startswith(">>>"):
                            break  # End of the file
                    
                    else:
                        # Skip to the next line
                        i += 1
                
                # Analyze the alignment for indels and stop codons
                if alignment["protein_seq"] and alignment["dna_seq"]:
                    # Count frameshifts (both forward '\' and backward '/')
                    forward_frameshifts = alignment["dna_seq"].count('\\')
                    backward_frameshifts = alignment["dna_seq"].count('/')
                    alignment["frameshifts"] = forward_frameshifts + backward_frameshifts

                    # Count stop codons in the pseudogene sequence
                    alignment["stop_codons"] = alignment["dna_seq"].count('*')
                    
                    # Count insertions (gaps in protein seq)
                    alignment["insertions"] = alignment["protein_seq"].count('-')
                    
                    # Count deletions (gaps in dna seq)
                    alignment["deletions"] = alignment["dna_seq"].count('-')

                # Add the alignment to our collection
                alignments.append(alignment)
            else:
                i += 1
        
        # Sort alignments by score (higher scores first)
        if alignments:
            alignments.sort(key=lambda x: x["score"], reverse=True)

            # Filter alignments by strand - keep only those matching the original strand
            alignments = [a for a in alignments if a["strand"] == exon_info['strand']]
            
            # Return only the best alignment
            return alignments[0]
        
    except Exception as e:
        logger.error(f"Error parsing tfasty output for {exon_info['id']}: {e}")
        import traceback
        logger.error(traceback.format_exc())
    
    return []

def _merge_exons(pseudogene_id, exon_results, chrom, protein, strand, protein_seqs):
    """Merge aligned exons from the same cluster into a pseudogene"""
    if not exon_results:
        return None
    
    # Initialize pseudogene data
    merged_pseudogene = {
        'id': pseudogene_id,
        'chrom': chrom,
        'protein': protein,
        'strand': strand,
        'start': float('inf'),
        'end': 0,
        'score': 0,
        'evalue': float('inf'),
        'identity': 0,
        'coverage': 0,
        'frameshifts': 0,
        'insertions': 0,
        'deletions': 0,
        'stop_codons': 0,
        'exons': []
    }
    
    # Track protein identity and overlap
    total_identity = 0
    pg_prot_start = float('inf')
    pg_prot_end = 0
    
    # Process each exon
    for exon_idx, exon in enumerate(exon_results):
        offset = exon["offset"]  # Local to genomic coordinate conversion
        alignment = exon["alignment"]
            
        # Get DNA coordinates
        dna_start = min(alignment["dna_range"].get("start", 0), alignment["dna_range"].get("end", 0))
        dna_end = max(alignment["dna_range"].get("start", 0), alignment["dna_range"].get("end", 0))
        
        # Convert to genomic coordinates
        genomic_start = offset + dna_start
        genomic_end = offset + dna_end
        
        # Update pseudogene bounds
        merged_pseudogene['start'] = min(merged_pseudogene['start'], genomic_start)
        merged_pseudogene['end'] = max(merged_pseudogene['end'], genomic_end)
        
        # Update statistics
        merged_pseudogene['score'] += alignment['score']
        merged_pseudogene['evalue'] = min(merged_pseudogene['evalue'], alignment['evalue'])
        merged_pseudogene['frameshifts'] += alignment['frameshifts']
        merged_pseudogene['insertions'] += alignment['insertions']
        merged_pseudogene['deletions'] += alignment['deletions']
        merged_pseudogene['stop_codons'] += alignment['stop_codons']
        
        # Track protein coverage
        exon_prot_start = min(alignment["protein_range"].get("start", 0), alignment["protein_range"].get("end", 0))
        exon_prot_end = max(alignment["protein_range"].get("start", 0), alignment["protein_range"].get("end", 0))
        pg_prot_start = min(pg_prot_start, exon_prot_start)
        pg_prot_end = max(pg_prot_end, exon_prot_end)

        # Update for average identity and coverage calculation
        total_identity += alignment["identity"]
        
        # Add this exon to the merged pseudogene
        merged_pseudogene['exons'].append({
            'id': f"{pseudogene_id}_exon{len(merged_pseudogene['exons'])+1}",
            'start': genomic_start,
            'end': genomic_end,
            'score': alignment['score'],
            'strand': strand,
            'protein_range': f"{exon_prot_start}-{exon_prot_end}",
            'identity': alignment['identity'],
            'frameshifts': alignment['frameshifts'],
            'insertions': alignment['insertions'],
            'deletions': alignment['deletions'],
            'stop_codons': alignment['stop_codons']
        })
    
    overlap = pg_prot_end - pg_prot_start

    # Calculate average identity 
    merged_pseudogene['identity'] = total_identity / len(exon_results)
    
    # Calculate avarage coverage
    protein_length = len(protein_seqs[protein])
    merged_pseudogene['coverage'] = (overlap / protein_length)
    
    # Sort exons by position
    merged_pseudogene['exons'].sort(key=lambda x: x['start'])
    
    return merged_pseudogene