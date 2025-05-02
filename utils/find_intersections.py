def parse_pseudopipe_genes(filename):
    genes = []
    
    with open(filename, 'r') as file:
        for line in file:
            # Skip comment lines
            if line.startswith('#'):
                continue
                
            columns = line.strip().split('\t')
            
            try:
                crom = columns[0]
                start = int(columns[1])
                end = int(columns[2])
                coverage = float(columns[5])
                identity = float(columns[11])
                    
                gene_info = {
                    'chrom': crom,
                    'start': start,
                    'end': end,
                    'length': end - start + 1,
                    'coverage': coverage,
                    'identity': identity,
                }
                genes.append(gene_info)
            except (ValueError, IndexError):
                # Skip lines with non-integer coordinates or insufficient columns
                continue
    
    return genes

def parse_tfasty_genes(filename):
    genes = []
    
    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
        
            try:
                crom = columns[1]
                start = int(columns[2])
                end = int(columns[3])
                coverage = float(columns[11])
                identity = float(columns[12])
                    
                gene_info = {
                    'chrom': crom,
                    'start': start,
                    'end': end,
                    'length': end - start + 1,
                    'coverage': coverage,
                    'identity': identity,
                }
                genes.append(gene_info)
            except (ValueError, IndexError):
                # Skip lines with non-integer coordinates or insufficient columns
                continue
    
    return genes

def calculate_overlap(pseudopipe_gene, pseudoscope_gene):
    """
    Calculate the overlap percentage between two genes.
    
    Args:
        pseudopipe_gene (dict): First gene information
        pseudoscope_gene (dict): Second gene information
        
    Returns:
        float: Overlap percentage relative to pseudopipe_gene
    """
    # Find the overlapping region
    overlap_start = max(pseudopipe_gene['start'], pseudoscope_gene['start'])
    overlap_end = min(pseudopipe_gene['end'], pseudoscope_gene['end'])
    
    # If there's no overlap
    if overlap_start > overlap_end:
        return 0.0
    
    overlap_length = overlap_end - overlap_start + 1
    overlap_percentage = (overlap_length / pseudopipe_gene['length']) * 100

    return overlap_percentage

def find_intersecting_genes(pseudopipe_genes, pseudoscope_genes, threshold=60.0):
    """
    Find genes from pseudopipe_genes that intersect with any gene in pseudoscope_genes.
    Returns the best match (highest overlap percentage) for each gene in pseudopipe_genes.
    
    Args:
        pseudopipe_genes (list): List of genes from file 1
        pseudoscope_genes (list): List of genes from file 2
        threshold (float): Minimum overlap percentage to consider a match
        
    Returns:
        list: List of dictionaries with matching gene information
    """
    matches = []
    
    for pseudopipe_gene in pseudopipe_genes:
        best_match = None
        best_overlap = 0.0
        
        for pseudoscope_gene in pseudoscope_genes:
            # Check if genes are on the same chromosome
            if pseudopipe_gene['chrom'] != pseudoscope_gene['chrom']:
                continue
            
            overlap_percentage = calculate_overlap(pseudopipe_gene, pseudoscope_gene)
            
            if overlap_percentage >= threshold and overlap_percentage > best_overlap:
                best_overlap = overlap_percentage
                best_match = {
                    'ppipe_start': pseudopipe_gene['start'],
                    'ppipe_end': pseudopipe_gene['end'],
                    'pscope_start': pseudoscope_gene['start'],
                    'pscope_end': pseudoscope_gene['end'],
                    'overlap': overlap_percentage,
                    'chrom': pseudoscope_gene['chrom'],
                    'ppipe_cov': pseudopipe_gene['coverage'],
                    'pscope_cov': pseudoscope_gene['coverage'],
                    'ppipe_ident': pseudopipe_gene['identity'],
                    'pscope_ident': pseudoscope_gene['identity'],
                }
        
        # Add the best match if one was found
        if best_match is not None:
            matches.append(best_match)
    
    return matches

def calculate_averages(matches):
    """
    Calculate average values for ppipe_cov, pscope_cov, ppipe_ident, and pscope_ident.
    
    Args:
        matches (list): List of match dictionaries
        
    Returns:
        dict: Dictionary containing average values
    """
    if not matches:
        return {
            'avg_ppipe_cov': 0.0,
            'avg_pscope_cov': 0.0,
            'avg_ppipe_ident': 0.0,
            'avg_pscope_ident': 0.0,
            'avg_overlap': 0.0
        }
    
    total_ppipe_cov = sum(match['ppipe_cov'] for match in matches)
    total_pscope_cov = sum(match['pscope_cov'] for match in matches)
    total_ppipe_ident = sum(match['ppipe_ident'] for match in matches)
    total_pscope_ident = sum(match['pscope_ident'] for match in matches)
    total_overlap = sum(match['overlap'] for match in matches)
    
    return {
        'avg_ppipe_cov': total_ppipe_cov / len(matches),
        'avg_pscope_cov': total_pscope_cov / len(matches),
        'avg_ppipe_ident': total_ppipe_ident / len(matches),
        'avg_pscope_ident': total_pscope_ident / len(matches),
        'avg_overlap': total_overlap / len(matches)
    }

def main():
    import sys
    
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <file1.gff> <file2.gff> <output_file>")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3]
    
    # Parse genes from both files
    pseudopipe_genes = parse_pseudopipe_genes(file1)
    pseudoscope_genes = parse_tfasty_genes(file2)
    
    print(f"Parsed {len(pseudopipe_genes)} genes from {file1}")
    print(f"Parsed {len(pseudoscope_genes)} genes from {file2}")
    
    # Find intersecting genes
    matches = find_intersecting_genes(pseudopipe_genes, pseudoscope_genes, threshold=60.0)
    
    # Calculate average values
    averages = calculate_averages(matches)
    
    # Calculate percentage of genes from file1 that have matches
    match_percentage = (len(matches) / len(pseudopipe_genes)) * 100 if pseudopipe_genes else 0
    
    # Write results to output file
    with open(output_file, 'w') as outfile:
        outfile.write(f"# Intersecting genes (threshold: 60%)\n")
        outfile.write(f"# Total genes in file1: {len(pseudopipe_genes)}\n")
        outfile.write(f"# Total genes in file2: {len(pseudoscope_genes)}\n")
        outfile.write(f"# Total matches found: {len(matches)}\n")
        outfile.write(f"# Match percentage: {match_percentage:.2f}%\n")
        outfile.write(f"# Average values:\n")
        outfile.write(f"#   pseudopipe coverage: {averages['avg_ppipe_cov']:.2f}\n")
        outfile.write(f"#   pseudoscope coverage: {averages['avg_pscope_cov']:.2f}\n")
        outfile.write(f"#   pseudopipe identity: {averages['avg_ppipe_ident']:.2f}\n")
        outfile.write(f"#   pseudoscope identity: {averages['avg_pscope_ident']:.2f}\n")
        outfile.write(f"#   overlap: {averages['avg_overlap']:.2f}%\n\n")
        
        outfile.write("# chrom\tppipe_start\tppipe_end\tpscope_start\tpscope_end\tppipe_cov\tpscope_cov\tppipe_ident\tpscope_ident\toverlap\n")
        
        for match in matches:
            outfile.write(f"{match['chrom']}\t{match['ppipe_start']}\t{match['ppipe_end']}\t{match['pscope_start']}\t{match['pscope_end']}\t{match['ppipe_cov']:.2f}\t{match['pscope_cov']:.2f}\t{match['ppipe_ident']:.2f}\t{match['pscope_ident']:.2f}\t{match['overlap']:.2f}%\n")
    
    print(f"Results written to {output_file}")
    print(f"Match percentage: {match_percentage:.2f}%")
    print(f"Average pseudopipe coverage: {averages['avg_ppipe_cov']:.2f}")
    print(f"Average pseudoscope coverage: {averages['avg_pscope_cov']:.2f}")
    print(f"Average pseudopipe identity: {averages['avg_ppipe_ident']:.2f}")
    print(f"Average pseudoscope identity: {averages['avg_pscope_ident']:.2f}")
    print(f"Average overlap: {averages['avg_overlap']:.2f}%")

if __name__ == "__main__":
    main()