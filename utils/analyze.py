import pandas as pd
import os
import subprocess
from Bio import SeqIO

class PseudogeneAnalyzer:
    def __init__(self, coordinate_file, genome_fasta, protein_fasta, output_dir="results"):
        """
        Initialize the pseudogene analyzer with input files.
        
        Args:
            coordinate_file: Path to file with gene coordinates
            genome_fasta: Path to genome FASTA file
            protein_fasta: Path to protein sequences FASTA file
            output_dir: Directory to store results
        """
        self.coordinate_file = coordinate_file
        self.genome_fasta = genome_fasta
        self.protein_fasta = protein_fasta
        self.output_dir = output_dir
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Load coordinates
        self.coordinates = self.load_coordinates()
        
        # Load genome
        self.genome_dict = self.load_genome()
        
        # Load protein sequences
        self.protein_dict = self.load_proteins()
        
        # Initialize result containers
        self.annotation_coverage = []
        self.pseudoscope_coverage = []
        self.blast_results = None
        
    def load_coordinates(self):
        """Load coordinate file, skipping header lines."""
        with open(self.coordinate_file, 'r') as f:
            lines = f.readlines()
        
        data_lines = [line for line in lines if not line.startswith('#')]
        
        coordinates = []
        for line in data_lines:
            parts = line.strip().split('\t')
            if len(parts) == 7:  # Ensure we have all expected columns
                coordinates.append({
                    'annotation_start': int(parts[0]),
                    'annotation_end': int(parts[1]),
                    'pseudoscope_start': int(parts[2]),
                    'pseudoscope_end': int(parts[3]),
                    'overlap_percentage': parts[4].strip('%'),
                    'chromosome': parts[5],
                    'pseudoscope_coverage': float(parts[6])
                })
        
        return pd.DataFrame(coordinates)
    
    def load_genome(self):
        """Load genome FASTA file into a dictionary."""
        genome_dict = {}
        for record in SeqIO.parse(self.genome_fasta, "fasta"):
            genome_dict[record.id] = str(record.seq)
        return genome_dict
    
    def load_proteins(self):
        """Load protein FASTA file into a dictionary and capture lengths."""
        protein_dict = {}
        protein_lengths = {}
        
        for record in SeqIO.parse(self.protein_fasta, "fasta"):
            protein_dict[record.id] = str(record.seq)
            protein_lengths[record.id] = len(record.seq)
            
        self.protein_lengths = protein_lengths
        return protein_dict
    
    def extract_sequence(self, chromosome, start, end):
        """Extract sequence from genome based on coordinates."""
        if chromosome in self.genome_dict:
            # Convert to 0-based indexing and extract
            return self.genome_dict[chromosome][start-1:end]
        else:
            print(f"Warning: Chromosome {chromosome} not found in genome FASTA")
            return ""
    
    def create_gene_fasta_files(self):
        """Extract sequences for annotation genes and pseudoscope genes."""
        # Create FASTA files for both gene types
        annotation_fasta_path = os.path.join(self.output_dir, "annotation_genes.fasta")
        pseudoscope_fasta_path = os.path.join(self.output_dir, "pseudoscope_genes.fasta")
        
        with open(annotation_fasta_path, 'w') as annot_file, \
             open(pseudoscope_fasta_path, 'w') as pseudo_file:
            
            for idx, row in self.coordinates.iterrows():
                chrom = row['chromosome']
                
                # Extract annotation gene sequence
                annot_seq = self.extract_sequence(chrom, row['annotation_start'], row['annotation_end'])
                annot_file.write(f">annot_{idx}|{chrom}:{row['annotation_start']}-{row['annotation_end']}\n")
                annot_file.write(f"{annot_seq}\n")
                
                # Extract pseudoscope gene sequence
                pseudo_seq = self.extract_sequence(chrom, row['pseudoscope_start'], row['pseudoscope_end'])
                pseudo_file.write(f">pseudo_{idx}|{chrom}:{row['pseudoscope_start']}-{row['pseudoscope_end']}\n")
                pseudo_file.write(f"{pseudo_seq}\n")
        
        print(f"Created FASTA files for annotation genes ({annotation_fasta_path}) and pseudoscope genes ({pseudoscope_fasta_path})")
        
        return annotation_fasta_path, pseudoscope_fasta_path
    
    def create_protein_database(self):
        """Create a BLAST database from protein sequences."""
        db_path = os.path.join(self.output_dir, "proteins_db")
        cmd = f"makeblastdb -in {self.protein_fasta} -dbtype prot -out {db_path}"
        
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"Created protein BLAST database at {db_path}")
            return db_path
        except subprocess.CalledProcessError as e:
            print(f"Failed to create BLAST database: {e}")
            return None
    
    def run_blastx(self, query_fasta, db_path, output_name):
        """Run BLASTX alignment of gene sequences to protein database."""
        output_path = os.path.join(self.output_dir, f"{output_name}_blast_results.txt")
        
        cmd = f"blastx -query {query_fasta} -db {db_path} -out {output_path} " \
              f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' " \
              f"-evalue 1e-5 -max_target_seqs 5"
        
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"BLASTX completed successfully. Results saved to {output_path}")
            return output_path
        except subprocess.CalledProcessError as e:
            print(f"BLASTX failed: {e}")
            return None
    
    def parse_blast_results(self, blast_file):
        """Parse BLAST output file into a pandas DataFrame."""
        cols = ["query", "subject", "pident", "length", "mismatch", "gapopen", 
                "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        
        try:
            results = pd.read_csv(blast_file, sep='\t', names=cols)
            return results
        except Exception as e:
            print(f"Failed to parse BLAST results: {e}")
            return pd.DataFrame(columns=cols)
    
    def calculate_coverage(self, blast_results, gene_type):
        """
        Calculate coverage of proteins by gene sequences, accounting for multiple HSPs.
        
        This improved version calculates:
        1. Individual hit coverage - each HSP's coverage of the protein
        2. Combined coverage - merged non-overlapping regions' coverage
        3. Total span coverage - from first hit to last hit (including gaps)
        """
        # Group results by query-subject pairs 
        grouped = blast_results.groupby(['query', 'subject'])
        
        coverage_data = []
        
        for (query, protein), group_df in grouped:
            # Extract gene ID from FASTA header
            gene_id = query.split('|')[0]
            
            # Get protein length
            protein_len = self.protein_lengths.get(protein, 0)
            
            if protein_len > 0:
                # Sort hits by position on protein
                group_df = group_df.sort_values('sstart')
                
                # Calculate individual hit coverages
                hit_coverages = []
                for _, hit in group_df.iterrows():
                    # Convert to protein coordinates for coverage calculation
                    prot_start = min(hit['sstart'], hit['send'])
                    prot_end = max(hit['sstart'], hit['send'])
                    hit_len = prot_end - prot_start + 1
                    hit_coverage = hit_len / protein_len
                    hit_coverages.append(hit_coverage)
                
                # Calculate average hit coverage
                avg_hit_coverage = sum(hit_coverages) / len(hit_coverages)
                
                # Calculate total span coverage (from first to last hit)
                first_start = min(group_df['sstart'].min(), group_df['send'].min())
                last_end = max(group_df['sstart'].max(), group_df['send'].max())
                total_span = last_end - first_start + 1
                total_span_coverage = min(1.0, total_span / protein_len)  # Cap at 100%
                
                # Calculate combined coverage (merge overlapping regions)
                covered_regions = []
                for _, hit in group_df.iterrows():
                    prot_start = min(hit['sstart'], hit['send'])
                    prot_end = max(hit['sstart'], hit['send'])
                    covered_regions.append((prot_start, prot_end))
                
                # Merge overlapping regions
                if covered_regions:
                    covered_regions.sort()
                    merged_regions = [covered_regions[0]]
                    
                    for current in covered_regions[1:]:
                        prev = merged_regions[-1]
                        
                        # If current overlaps with previous
                        if current[0] <= prev[1]:
                            # Merge them
                            merged_regions[-1] = (prev[0], max(prev[1], current[1]))
                        else:
                            # Add as separate region
                            merged_regions.append(current)
                    
                    # Calculate total coverage of merged regions
                    total_covered = sum(end - start + 1 for start, end in merged_regions)
                    combined_coverage = total_covered / protein_len
                else:
                    combined_coverage = 0.0
                
                # Calculate average identity (weighted by aligned segment length)
                total_length = group_df['length'].sum()
                weighted_identity = sum(hit['pident'] * hit['length'] for _, hit in group_df.iterrows()) / total_length
                
                # Get overall best e-value
                best_evalue = group_df['evalue'].min()
                
                coverage_data.append({
                    "gene_id": gene_id,
                    "protein_id": protein,
                    "identity": weighted_identity,
                    "avg_hit_coverage": avg_hit_coverage,  # Average of individual hit coverages
                    "combined_coverage": combined_coverage,  # Coverage after merging overlapping regions
                    "total_span_coverage": total_span_coverage,  # Coverage from first hit to last hit
                    "evalue": best_evalue,
                    "segments": len(group_df),  # Number of HSPs
                    "gene_type": gene_type,
                    "score": combined_coverage * weighted_identity  # Score for ranking matches
                })
        
        # For each query, find the best matching protein by combined_coverage * identity
        best_matches = []
        for query in blast_results['query'].unique():
            query_matches = [m for m in coverage_data if m['gene_id'] == query.split('|')[0]]
            if query_matches:
                # Sort by score (combined_coverage * identity)
                query_matches.sort(key=lambda x: x['score'], reverse=True)
                best_matches.append(query_matches[0])
        
        coverage_df = pd.DataFrame(best_matches)
        
        if gene_type == "annotation":
            self.annotation_coverage = coverage_df
        else:
            self.pseudoscope_coverage = coverage_df
            
        return coverage_df
    
    def compare_coverage(self):
        """Compare coverage statistics between annotation and pseudoscope genes."""
        if self.annotation_coverage.empty or self.pseudoscope_coverage.empty:
            print("Coverage data not available. Run alignment analysis first.")
            return None
        
        # Combined coverage data
        all_coverage = pd.concat([self.annotation_coverage, self.pseudoscope_coverage])
        
        # Summary statistics for all coverage metrics
        coverage_metrics = ['avg_hit_coverage', 'combined_coverage', 'total_span_coverage', 'identity', 'segments']
        
        # Create summary statistics
        summary = {}
        for metric in coverage_metrics:
            summary[metric] = all_coverage.groupby('gene_type')[metric].agg(['mean', 'median', 'std', 'min', 'max'])
        
        # Convert to DataFrame
        summary_df = pd.DataFrame()
        for metric, stats in summary.items():
            stats.columns = [f"{metric}_{col}" for col in stats.columns]
            if summary_df.empty:
                summary_df = stats
            else:
                summary_df = pd.concat([summary_df, stats], axis=1)
        
        # Save summary to file
        summary_path = os.path.join(self.output_dir, "coverage_summary.csv")
        summary_df.to_csv(summary_path)
        
        # Save detailed coverage data
        detailed_path = os.path.join(self.output_dir, "detailed_coverage.csv")
        all_coverage.to_csv(detailed_path, index=False)
        
        print(f"Coverage summary saved to {summary_path}")
        print(f"Detailed coverage data saved to {detailed_path}")
        
        
        return summary_df
    
    def find_best_protein_matches(self):
        """Find the best protein match for each gene."""
        if self.annotation_coverage.empty or self.pseudoscope_coverage.empty:
            print("Coverage data not available. Run alignment analysis first.")
            return None
        
        # Extract gene IDs from query strings
        def extract_id(gene_id):
            return int(gene_id.split('_')[1])
        
        self.annotation_coverage['original_idx'] = self.annotation_coverage['gene_id'].apply(extract_id)
        self.pseudoscope_coverage['original_idx'] = self.pseudoscope_coverage['gene_id'].apply(extract_id)
        
        # Merge results to compare protein matches
        merged = pd.merge(
            self.annotation_coverage, 
            self.pseudoscope_coverage,
            on='original_idx',
            suffixes=('_annotation', '_pseudoscope')
        )
        
        # Check if genes match to the same protein
        merged['same_protein'] = merged['protein_id_annotation'] == merged['protein_id_pseudoscope']
        
        # Calculate match statistics
        match_stats = {
            'total_pairs': len(merged),
            'same_protein_matches': merged['same_protein'].sum(),
            'match_percentage': merged['same_protein'].mean() * 100,
            'avg_annotation_hit_coverage': merged['avg_hit_coverage_annotation'].mean(),
            'avg_pseudoscope_hit_coverage': merged['avg_hit_coverage_pseudoscope'].mean(),
            'avg_annotation_combined_coverage': merged['combined_coverage_annotation'].mean(),
            'avg_pseudoscope_combined_coverage': merged['combined_coverage_pseudoscope'].mean(),
            'avg_annotation_span_coverage': merged['total_span_coverage_annotation'].mean(),
            'avg_pseudoscope_span_coverage': merged['total_span_coverage_pseudoscope'].mean(),
            'avg_annotation_identity': merged['identity_annotation'].mean(),
            'avg_pseudoscope_identity': merged['identity_pseudoscope'].mean(),
            'avg_annotation_segments': merged['segments_annotation'].mean(),
            'avg_pseudoscope_segments': merged['segments_pseudoscope'].mean()
        }
        
        # Save detailed match data
        match_path = os.path.join(self.output_dir, "protein_matches.csv")
        merged.to_csv(match_path, index=False)
        
        # Save match statistics
        stats_path = os.path.join(self.output_dir, "match_statistics.txt")
        with open(stats_path, 'w') as f:
            f.write("Pseudogene-Annotation Protein Match Statistics\n")
            f.write("=============================================\n")
            f.write(f"Total gene pairs analyzed: {match_stats['total_pairs']}\n")
            f.write(f"Pairs matching same protein: {match_stats['same_protein_matches']}\n")
            f.write(f"Match percentage: {match_stats['match_percentage']:.2f}%\n\n")
            
            f.write("Coverage Statistics:\n")
            f.write("------------------\n")
            f.write(f"Average annotation gene hit coverage: {match_stats['avg_annotation_hit_coverage']:.4f}\n")
            f.write(f"Average pseudoscope gene hit coverage: {match_stats['avg_pseudoscope_hit_coverage']:.4f}\n")
            f.write(f"Ratio (pseudoscope/annotation): {match_stats['avg_pseudoscope_hit_coverage']/match_stats['avg_annotation_hit_coverage']:.4f}\n\n")
            
            f.write(f"Average annotation gene combined coverage: {match_stats['avg_annotation_combined_coverage']:.4f}\n")
            f.write(f"Average pseudoscope gene combined coverage: {match_stats['avg_pseudoscope_combined_coverage']:.4f}\n")
            f.write(f"Ratio (pseudoscope/annotation): {match_stats['avg_pseudoscope_combined_coverage']/match_stats['avg_annotation_combined_coverage']:.4f}\n\n")
            
            f.write(f"Average annotation gene span coverage: {match_stats['avg_annotation_span_coverage']:.4f}\n")
            f.write(f"Average pseudoscope gene span coverage: {match_stats['avg_pseudoscope_span_coverage']:.4f}\n")
            f.write(f"Ratio (pseudoscope/annotation): {match_stats['avg_pseudoscope_span_coverage']/match_stats['avg_annotation_span_coverage']:.4f}\n\n")
            
            f.write("Other Statistics:\n")
            f.write("---------------\n")
            f.write(f"Average annotation gene identity: {match_stats['avg_annotation_identity']:.2f}%\n")
            f.write(f"Average pseudoscope gene identity: {match_stats['avg_pseudoscope_identity']:.2f}%\n")
            f.write(f"Ratio (pseudoscope/annotation): {match_stats['avg_pseudoscope_identity']/match_stats['avg_annotation_identity']:.4f}\n\n")
            
            f.write(f"Average annotation gene segments: {match_stats['avg_annotation_segments']:.2f}\n")
            f.write(f"Average pseudoscope gene segments: {match_stats['avg_pseudoscope_segments']:.2f}\n")
            f.write(f"Ratio (pseudoscope/annotation): {match_stats['avg_pseudoscope_segments']/match_stats['avg_annotation_segments']:.4f}\n")
        
        print(f"Match statistics saved to {stats_path}")
        print(f"Detailed match data saved to {match_path}")
        
        return match_stats
    
    def run_analysis(self):
        """Run the complete analysis pipeline."""
        print("Starting pseudogene analysis pipeline...")
        
        # Step 1: Extract sequences
        annot_fasta, pseudo_fasta = self.create_gene_fasta_files()
        
        # Step 2: Create protein database
        db_path = self.create_protein_database()
        if not db_path:
            print("Failed to create protein database. Aborting analysis.")
            return
        
        # Step 3: Run BLASTX for annotation genes
        annot_blast = self.run_blastx(annot_fasta, db_path, "annotation")
        if not annot_blast:
            print("Failed to run BLASTX for annotation genes. Aborting analysis.")
            return
        
        # Step 4: Run BLASTX for pseudoscope genes
        pseudo_blast = self.run_blastx(pseudo_fasta, db_path, "pseudoscope")
        if not pseudo_blast:
            print("Failed to run BLASTX for pseudoscope genes. Aborting analysis.")
            return
        
        # Step 5: Parse BLAST results
        annot_results = self.parse_blast_results(annot_blast)
        pseudo_results = self.parse_blast_results(pseudo_blast)
        
        # Step 6: Calculate coverage
        self.calculate_coverage(annot_results, "annotation")
        self.calculate_coverage(pseudo_results, "pseudoscope")
        
        # Step 7: Compare coverage
        coverage_summary = self.compare_coverage()
        
        # Step 8: Find best protein matches
        match_stats = self.find_best_protein_matches()
        
        print("Analysis complete!")
        return {
            "coverage_summary": coverage_summary,
            "match_stats": match_stats
        }


# Example usage
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze pseudogenes and their protein matches")
    parser.add_argument("--coord", required=True, help="File with gene coordinates")
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--proteins", required=True, help="Protein FASTA file")
    parser.add_argument("--outdir", default="results", help="Output directory")
    
    args = parser.parse_args()
    
    analyzer = PseudogeneAnalyzer(
        coordinate_file=args.coord,
        genome_fasta=args.genome,
        protein_fasta=args.proteins,
        output_dir=args.outdir
    )
    
    # Run the complete analysis
    results = analyzer.run_analysis()