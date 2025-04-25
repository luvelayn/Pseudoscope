import os
import subprocess
import tempfile
import csv
from Bio import SeqIO

def run_exonerate(pseudogenes_file, protein_file, genome_file, out_dir, logger):
        """Run precise re-alignment with exonerate for each pseudogene candidate"""
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
                    
                # Add flanking sequence (100 nt on each side)
                flank = 100
                start = max(0, pg['start'] - flank)
                end = min(len(genome_seqs[chrom]), pg['end'] + flank)
                
                seq = genome_seqs[chrom][start:end]
                protein = pg['protein']
                
                if protein not in protein_seqs:
                    logger.warning(f"Protein {protein} not found, skipping pseudogene")
                    continue
                
                pseudogene_regions.append({
                    'id': f"pseudogene_{i+1}",
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'seq': seq,
                    'protein': protein,
                    'protein_seq': protein_seqs[protein],
                    'strand': pg['strand']
                })
            except Exception as e:
                logger.error(f"Error processing pseudogene {i}: {e}")
                
        logger.info(f"Prepared {len(pseudogene_regions)} pseudogene regions for re-alignment")
                
        # Create temporary files for sequences
        pseudogene_seqs_file = os.path.join(out_dir, "pseudogene_candidates_seqs.fa")
        with open(pseudogene_seqs_file, 'w') as f:
            for pg in pseudogene_regions:
                f.write(f">{pg['id']} chrom={pg['chrom']} start={pg['start']} end={pg['end']} protein={pg['protein']}\n")
                f.write(f"{pg['seq']}\n")
        
        # Run exonerate for each pseudogene-protein pair
        results = []
        for pg in pseudogene_regions:
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as pg_file, \
                 tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as prot_file, \
                 tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as out_file:
                
                pg_file.write(f">{pg['id']}\n{pg['seq']}\n")
                prot_file.write(f">{pg['protein']}\n{pg['protein_seq']}\n")
                
                pg_file_path = pg_file.name
                prot_file_path = prot_file.name
                out_file_path = out_file.name
            
            # Run exonerate
            exonerate_cmd = [
                'exonerate',
                '--model', 'protein2genome',
                '--showvulgar', 'no',
                '--showalignment', 'no',
                '--showtargetgff', 'yes',
                '--bestn', '1',
                '--refine', 'full',
                prot_file_path,
                pg_file_path
            ]
            
            try:
                result = subprocess.run(exonerate_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                
                # Parse results
                gff_lines = []
                for line in result.stdout.split('\n'):
                    if not line.startswith('#') and line.strip():
                        fields = line.strip().split('\t')
                        if len(fields) >= 8:
                            # Adjust coordinates to genome coordinates
                            if fields[0] == pg['id']:
                                fields[0] = pg['chrom']
                                fields[3] = str(int(fields[3]) + pg['start'])
                                fields[4] = str(int(fields[4]) + pg['start'])
                                gff_lines.append('\t'.join(fields))
                
                # Save results
                with open(out_file_path, 'w') as out:
                    for line in gff_lines:
                        out.write(line + '\n')
                
                results.append({
                    'id': pg['id'],
                    'chrom': pg['chrom'],
                    'protein': pg['protein'],
                    'gff_file': out_file_path,
                    'gff_lines': gff_lines
                })
                
            except subprocess.CalledProcessError as e:
                logger.error(f"Error running exonerate for {pg['id']}: {e}")
                logger.error(e.stderr)
            
            # Clean up temporary files
            try:
                os.unlink(pg_file_path)
                os.unlink(prot_file_path)
                # Keep out_file_path for now
            except:
                pass
                
        # Combine all GFF files
        exonerate_out = os.path.join(out_dir, "exonerate_results.gff")

        with open(exonerate_out, 'w') as f:
            f.write("##gff-version 3\n")
            for result in results:
                for line in result['gff_lines']:
                    f.write(line + '\n')
                    
        logger.info(f"Completed exonerate re-alignment for {len(results)} pseudogenes")
        
        # Clean up individual result files
        for result in results:
            try:
                os.unlink(result['gff_file'])
            except:
                pass