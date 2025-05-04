import os

def generate_output(pseudogene_results, genome_seqs, output_dir):
    # Final output paths
    pseudogenes_gff = os.path.join(output_dir, "pseudogenes.gff")
    pseudogenes_tsv = os.path.join(output_dir, "pseudogenes.tsv")
    pseudogene_fasta = os.path.join(output_dir, "pseudogenes.fa")

    # Write gff file
    with open(pseudogenes_gff, 'w') as f:
        f.write("##gff-version 3\n")
        # Write pseudogene line
        for pg in pseudogene_results:
            f.write(f"{pg['chrom']}\tpseudoscope\tpseudogene\t{pg['start']}\t"
                    f"{pg['end']}\t{pg['score']}\t{pg['strand']}\t.\t"
                    f"ID={pg['id']};Parent={pg['protein']};"
                    f"frameshifts={pg['frameshifts']};insertions={pg['insertions']};"
                    f"deletions={pg['deletions']};stop_codons={pg['stop_codons']};"
                    f"exons={pg['exon_count']};type={pg['type']}\n")
            # Write exon lines
            for exon in pg['exons']:
                f.write(f"{pg['chrom']}\tpseudoscope\texon\t{exon['start']}\t"
                        f"{exon['end']}\t{exon['score']}\t{exon['strand']}\t.\t"
                        f"ID={exon['id']};Parent={pg['id']};"
                        f"protein_range={exon['protein_range']};"
                        f"frameshifts={exon['frameshifts']};insertions={exon['insertions']};"
                        f"deletions={exon['deletions']};stop_codons={exon['stop_codons']}\n")

    # Write tsv file
    with open(pseudogenes_tsv, 'w') as f:
        f.write("# id\tchrom\tstart\tend\tstrand\tparent_protein\t"
                "frameshifts\tinsertions\tdeletions\tstop_codons\t"
                "exons\tcoverage\tidentity\tevalue\ttype\n")
        
        for pg in pseudogene_results:
            f.write(f"{pg['id']}\t{pg['chrom']}\t{pg['start']}\t{pg['end']}\t{pg['strand']}\t{pg['protein']}\t"
                    f"{pg['frameshifts']}\t{pg['insertions']}\t{pg['deletions']}\t{pg['stop_codons']}\t{pg['exon_count']}\t"
                    f"{pg['coverage']:.2f}\t{pg['identity']:.2f}\t{pg['evalue']}\t{pg['type']}\n")
    
    # Write fasta file
    with open(pseudogene_fasta, 'w') as f:
        for pg in pseudogene_results:
            f.write(f">{pg['id']}\n")
            f.write(f"{genome_seqs[pg['chrom']][pg['start']:pg['end']]}\n")