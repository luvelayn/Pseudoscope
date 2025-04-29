import os
import subprocess
import tempfile
import csv
import re
from Bio import SeqIO

def run_tfasty(pseudogenes_file, protein_file, genome_file, out_dir, logger):
        """Run precise re-alignment with tfasty for each pseudogene candidate"""
        pseudogenes = []

        with open(pseudogenes_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for i, row in enumerate(reader, 1):
                try:
                    # Convert start and end to integers
                    start = int(row['start'])
                    end = int(row['end'])
                    
                    # Convert evalue and coverage to float
                    evalue = float(row['evalue'])
                    coverage = float(row['coverage'])
                    strand = row['strand']
                    
                    pseudogenes.append({
                        'protein': row['protein'],
                        'chrom': row['chrom'],
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'evalue': evalue,
                        'coverage': coverage
                    })
                    
                except ValueError as e:
                    logger.warning(f"Error parsing line {i}: {e}")
                except Exception as e:
                    logger.error(f"Unexpected error processing line {i}: {e}")

        if not pseudogenes:
            logger.warning("No pseudogene candidates for re-alignment")
            return

        # Extract protein sequences
        protein_seqs = {}
        for record in SeqIO.parse(protein_file, "fasta"):
            protein_seqs[record.id] = str(record.seq)
            
        # Extract genome sequences
        genome_seqs = {}
        for record in SeqIO.parse(genome_file, "fasta"):
            genome_seqs[record.id] = str(record.seq)
            
        # Extract pseudogene candidates with flanking sequence
        pseudogene_regions = []
        for i, pg in enumerate(pseudogenes):
            try:
                chrom = pg['chrom']
                if chrom not in genome_seqs:
                    logger.warning(f"Chromosome {chrom} not found in genome, skipping pseudogene")
                    continue
                    
                # Add flanking sequence (200 nt on each side for tfasty)
                flank = 200
                start = max(0, pg['start'] - flank)
                end = min(len(genome_seqs[chrom]), pg['end'] + flank)
                
                seq = genome_seqs[chrom][start:end]
                protein = pg['protein']
                
                if protein not in protein_seqs:
                    logger.warning(f"Protein {protein} not found, skipping pseudogene")
                    continue
                
                pseudogene_regions.append({
                    'id': f"PG{i+1}",
                    'chrom': chrom,
                    'start': start,
                    'original_start': pg['start'],  # Keep original coordinates
                    'end': end,
                    'original_end': pg['end'],      # Keep original coordinates 
                    'seq': seq,
                    'protein': protein,
                    'protein_seq': protein_seqs[protein],
                    'strand': pg['strand']
                })
            except Exception as e:
                logger.error(f"Error processing pseudogene {i}: {e}")
                
        logger.info(f"Prepared {len(pseudogene_regions)} pseudogene regions for re-alignment")
                
        # Create temporary file for all pseudogene sequences
        pseudogene_seqs_file = os.path.join(out_dir, "pseudogene_candidates_seqs.fa")
        with open(pseudogene_seqs_file, 'w') as f:
            for pg in pseudogene_regions:
                f.write(f">{pg['id']} chrom={pg['chrom']} start={pg['start']} end={pg['end']} protein={pg['protein']}\n")
                f.write(f"{pg['seq']}\n")
        
        # Run tfasty for each pseudogene-protein pair
        results = []
        tfasty_output_dir = os.path.join(out_dir, "tfasty_alignments")
        os.makedirs(tfasty_output_dir, exist_ok=True)
        
        for pg in pseudogene_regions:
            # Create temporary files for protein and dna sequences
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as pg_file, \
                 tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as prot_file:
                
                pg_file.write(f">{pg['id']}\n{pg['seq']}\n")
                prot_file.write(f">{pg['protein']}\n{pg['protein_seq']}\n")
                
                pg_file_path = pg_file.name
                prot_file_path = prot_file.name
            
            # Output for tfasty results
            alignment_out = os.path.join(tfasty_output_dir, f"{pg['id']}_alignment.txt")
            gff_out = os.path.join(tfasty_output_dir, f"{pg['id']}_alignment.gff")
            
            # Set up tfasty commad
            tfasty_cmd = f"tfasty36 -m 10 {prot_file_path} {pg_file_path} > {alignment_out}"
            
            try:
                logger.debug(f"Running command: {tfasty_cmd}")
                # Use shell=True to handle the redirection properly
                subprocess.run(tfasty_cmd, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                
                # Parse tfasty output and convert to GFF format
                gff_lines = _parse_tfasty_to_gff(alignment_out, pg, logger)
                
                # Save GFF results
                with open(gff_out, 'w') as out:
                    for line in gff_lines:
                        out.write(line + '\n')
                
                results.append({
                    'id': pg['id'],
                    'chrom': pg['chrom'],
                    'protein': pg['protein'],
                    'gff_file': gff_out,
                    'gff_lines': gff_lines
                })
                
            except subprocess.CalledProcessError as e:
                logger.error(f"Error running tfasty for {pg['id']}: {e}")
                if hasattr(e, 'stderr'):
                    logger.error(e.stderr)
            
            # Clean up temporary files
            try:
                os.unlink(pg_file_path)
                os.unlink(prot_file_path)
            except:
                pass
                
        # Combine all GFF files
        tfasty_out = os.path.join(out_dir, "tfasty_results.gff")

        with open(tfasty_out, 'w') as f:
            f.write("##gff-version 3\n")
            for result in results:
                for line in result['gff_lines']:
                    f.write(line + '\n')
                    
        logger.info(f"Completed tfasty re-alignment for {len(results)} pseudogenes")
        
        # Generate a summary file
        summary_file = os.path.join(out_dir, "pseudogene_summary.tsv")
        with open(summary_file, 'w') as f:
            f.write("pseudogene_id\tchrom\tstart\tend\tstrand\tparent_protein\t"
                    "frameshifts\tinsertions\tdeletions\tstop_codons\t"
                    "exon_count\tcoverage\tidentity\tevalue\n")
            
            for result in results:
                if not result['gff_lines']:
                    continue
                
                # Extract pseudogene features
                pseudogene_line = None
                exon_count = 0
                
                for line in result['gff_lines']:
                    if "\tpseudogene\t" in line:
                        pseudogene_line = line
                    elif "\texon\t" in line:
                        exon_count += 1
                
                if pseudogene_line:
                    parts = pseudogene_line.strip().split('\t')
                    if len(parts) >= 9:
                        chrom = parts[0]
                        start = parts[3]
                        end = parts[4]
                        strand = parts[6]
                        attr = parts[8]
                        
                        # Extract attributes
                        id_match = re.search(r'ID=([^;]+)', attr)
                        parent_match = re.search(r'Parent=([^;]+)', attr)
                        evalue_match = re.search(r'evalue=([^;]+)', attr)
                        identity_match = re.search(r'identity=([^;]+)', attr)
                        coverage_match = re.search(r'coverage=([^;]+)', attr)
                        frameshifts_match = re.search(r'frameshifts=([^;]+)', attr)
                        insertions_match = re.search(r'insertions=([^;]+)', attr)
                        deletions_match = re.search(r'deletions=([^;]+)', attr)
                        stop_codons_match = re.search(r'stop_codons=([^;]+)', attr)
                        
                        pg_id = id_match.group(1) if id_match else result['id']
                        parent = parent_match.group(1) if parent_match else result['protein']
                        evalue = evalue_match.group(1) if evalue_match else "NA"
                        identity = identity_match.group(1) if identity_match else "NA"
                        coverage = coverage_match.group(1) if coverage_match else "NA"
                        frameshifts = frameshifts_match.group(1) if frameshifts_match else "0"
                        insertions = insertions_match.group(1) if insertions_match else "0"
                        deletions = deletions_match.group(1) if deletions_match else "0"
                        stop_codons = stop_codons_match.group(1) if stop_codons_match else "0"
                        
                        f.write(f"{pg_id}\t{chrom}\t{start}\t{end}\t{strand}\t{parent}\t"
                                f"{frameshifts}\t{insertions}\t{deletions}\t{stop_codons}\t"
                                f"{exon_count}\t{coverage}\t{identity}\t{evalue}\n")
        
        return tfasty_out

def _parse_tfasty_to_gff(tfasty_output, pseudogene_info, logger):
    """Parse tfasty output (format 10) and convert to GFF3 format"""
    gff_lines = []
    
    try:
        # Read the tfasty output file
        with open(tfasty_output, 'r') as f:
            content = f.readlines()
        
        # Storage for alignments
        alignments = []
        
        # Extract ids
        pseudogene_id = pseudogene_info['id']
        protein_id = pseudogene_info['protein']
        
        # Parse the file
        i = 0
        while not content[i].startswith(">>><<<"):
            line = content[i].strip()
            
            # Look for the start of a pseudogene alignment block
            if line.startswith(f">>{pseudogene_id}") or line.startswith(">--"):
                # Initialize alignment data structure
                alignment = {
                    "name": pseudogene_id,
                    "strand": "",
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
                while not content[i].startswith(f">>{pseudogene_id}") and not content[i].startswith(">--"):
                    line = content[i].strip()
                    
                    # Strand information
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
                        while content[i].startswith(";"):
                            line = content[i].strip()
                            if line.startswith("; al_start:"):
                                alignment["protein_range"]["start"] = int(line.split(":", 1)[1].strip())
                            elif line.startswith("; al_stop:"):
                                alignment["protein_range"]["end"] = int(line.split(":", 1)[1].strip())
                            i += 1
                        
                        # Collect the protein sequence (may span multiple lines)
                        protein_seq = ""
                        while not content[i].startswith(">"):
                            protein_seq += content[i].strip()
                            i += 1
                        alignment["protein_seq"] = protein_seq
                        
                        # Get to the dna section (skip > line)
                        i += 1

                        # Get dna coordinates
                        while content[i].startswith(";"):
                            line = content[i].strip()
                            if line.startswith("; al_start:"):
                                alignment["dna_range"]["start"] = int(line.split(":", 1)[1].strip())
                            elif line.startswith("; al_stop:"):
                                alignment["dna_range"]["end"] = int(line.split(":", 1)[1].strip())
                            i += 1
                        
                        # Collect the DNA sequence (may span multiple lines)
                        dna_seq = ""
                        while not content[i].startswith(";"):
                            dna_seq += content[i].strip()
                            i += 1
                        alignment["dna_seq"] = dna_seq

                        # Skip to the next alignment block
                        while not content[i].startswith(">"):
                            i += 1
                        
                        if content[i].startswith(">>>"):
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
        alignments.sort(key=lambda x: x["score"], reverse=True)

        # Filter alignments by strand - keep only those matching the highest-scoring alignment's strand
        strand = alignments[0]["strand"]  # Get strand of highest-scoring alignment
        alignments = [a for a in alignments if a["strand"] == strand]
        
        # Filter out overlapping alignments (keep the one with higher score)
        non_overlapping = []
        for alignment in alignments:
            dna_start = min(alignment["dna_range"]["start"], alignment["dna_range"]["end"])
            dna_end = max(alignment["dna_range"]["start"], alignment["dna_range"]["end"])

            prot_start = min(alignment["protein_range"]["start"], alignment["protein_range"]["end"])
            prot_end = max(alignment["protein_range"]["start"], alignment["protein_range"]["end"])
            
            # Check for overlap with existing alignments
            overlap = False
            for existing in non_overlapping:
                existing_dna_start = min(existing["dna_range"]["start"], existing["dna_range"]["end"])
                existing_dna_end = max(existing["dna_range"]["start"], existing["dna_range"]["end"])

                existing_prot_start = min(existing["protein_range"]["start"], existing["protein_range"]["end"])
                existing_prot_end = max(existing["protein_range"]["start"], existing["protein_range"]["end"])
                
                # Check for overlap in either DNA or protein ranges
                if (dna_start <= existing_dna_end and dna_end >= existing_dna_start
                    or prot_start <= existing_prot_end and prot_end >= existing_prot_start):
                    overlap = True
                    break
            
            if not overlap:
                non_overlapping.append(alignment)
        
        # Re-sort non-overlapping alignments by position
        non_overlapping.sort(key=lambda x: min(x["dna_range"].get("start", 0), x["dna_range"].get("end", 0)))
        
        # If no valid alignments remain, return empty results
        if not non_overlapping:
            return gff_lines
        
        # Calculate gene boundaries and total statistics
        gene_start = float('inf')
        gene_end = 0
        total_score = 0
        min_evalue = float('inf')
        total_identity = 0
        total_overlap = 0
        total_frameshifts = 0 
        total_insertions = 0
        total_deletions = 0
        total_stop_codons = 0
        
        for alignment in non_overlapping:
            dna_start = min(alignment["dna_range"]["start"], alignment["dna_range"]["end"])
            dna_end = max(alignment["dna_range"]["start"], alignment["dna_range"]["end"])
            
            gene_start = min(gene_start, dna_start)
            gene_end = max(gene_end, dna_end)
            
            total_score += alignment["score"]
            min_evalue = min(min_evalue, alignment["evalue"])
            total_identity += alignment["identity"] * alignment["overlap"]  # Weight by overlap
            total_overlap += alignment["overlap"]
            total_frameshifts += alignment["frameshifts"] 
            total_insertions += alignment["insertions"]
            total_deletions += alignment["deletions"] 
            total_stop_codons += alignment["stop_codons"]
        
        # Calculate average identity across all blocks
        avg_identity = total_identity / total_overlap if total_overlap > 0 else 0
        
        # Calculate coverage based on protein length
        protein_length = len(pseudogene_info['protein_seq'])
        coverage = (total_overlap / protein_length) if protein_length > 0 else 0
        
        # Generate GFF for the gene and its exons
        chrom = pseudogene_info['chrom']
        local_offset = pseudogene_info['start']
        
        # Calculate genome coordinates
        pseudogene_start = local_offset + gene_start
        pseudogene_end = local_offset + gene_end

        gff_lines.append(
            f"{chrom}\tpseudoscope\tpseudogene\t{pseudogene_start}\t{pseudogene_end}\t{total_score}\t{strand}\t.\t"
            f"ID={pseudogene_id};Parent={protein_id};evalue={min_evalue};"
            f"identity={avg_identity:.2f};coverage={coverage:.2f};"
            f"frameshifts={total_frameshifts};insertions={total_insertions};"
            f"deletions={total_deletions};stop_codons={total_stop_codons}"
        )
        
        # Create exon features for each alignment block
        for i, alignment in enumerate(non_overlapping):
            dna_start = min(alignment["dna_range"]["start"], alignment["dna_range"]["end"])
            dna_end = max(alignment["dna_range"]["start"], alignment["dna_range"]["end"])
            
            exon_start = local_offset + dna_start
            exon_end = local_offset + dna_end
            
            prot_start = alignment["protein_range"].get("start", 0)
            prot_end = alignment["protein_range"].get("end", 0)
            
            gff_lines.append(
                f"{chrom}\tpseudoscope\texon\t{exon_start}\t{exon_end}\t{alignment['score']}\t{strand}\t.\t"
                f"ID={pseudogene_id}_exon{i+1};Parent={pseudogene_id};"
                f"protein_range={prot_start}-{prot_end};identity={alignment['identity']:.2f};"
                f"frameshifts={alignment['frameshifts']};insertions={alignment['insertions']};"
                f"deletions={alignment['deletions']};stop_codons={alignment['stop_codons']}"
            )
    
    except Exception as e:
        logger.error(f"Error parsing tfasty output for {pseudogene_info['id']}: {e}")
        import traceback
        logger.error(traceback.format_exc())
    
    return gff_lines