import pandas as pd
import argparse

def load_bed12(file_path):
    """Load BED12 file into a DataFrame and parse relevant columns."""
    df = pd.read_csv(
        file_path, sep='\t', header=None,
        names=["chrom", "start", "end", "name", "score", "strand", "thick_start", 
               "thick_end", "item_rgb", "block_count", "block_sizes", "block_starts"]
    )
    return df

def parse_junctions(row):
    """Parse internal junctions and additional criteria for comparison."""
    chrom = row["chrom"]
    block_count = row["block_count"]
    block_sizes = [int(size) for size in row["block_sizes"].split(",") if size]
    block_starts = [int(start) for start in row["block_starts"].split(",") if start]
    
    # Calculate end of first exon and start of last exon
    first_exon_end = row["start"] + block_sizes[0]
    last_exon_start = row["start"] + block_starts[-1]

    if block_count == 1:
        # Skip single-exon isoforms
        return None
    
    elif block_count == 2:
        # Only return minimal information for two-exon isoforms
        return chrom, block_count, None, row["start"], row["end"], first_exon_end, last_exon_start

    elif block_count > 2:
        # Calculate internal junctions (start and end of each block except first and last)
        junctions = []
        for i in range(1, block_count - 1):  # Exclude first and last exon
            exon_start = row["start"] + block_starts[i]
            exon_end = exon_start + block_sizes[i]
            junctions.append((exon_start, exon_end))
        
        # Output chrom, exon count, list of junctions, start, end, first exon end, last exon start
        return chrom, block_count, tuple(junctions), row["start"], row["end"], first_exon_end, last_exon_start

def find_matching_isoforms(file1, file2, outfile, max_diff=500):
    """Find isoforms with matching splice junctions, start/end within 500bp, and matching first/last exons."""
    with open(outfile, "w") as out:
        out.write(file1 + "\t" + file2 + "\n")
        
        # Load and parse files
        df1 = load_bed12(file1)
        df2 = load_bed12(file2)
        
        # Create dictionaries to store internal junctions by isoform name
        junctions1 = {row["name"]: parse_junctions(row) for _, row in df1.iterrows() if parse_junctions(row)}
        junctions2 = {row["name"]: parse_junctions(row) for _, row in df2.iterrows() if parse_junctions(row)}
        
        # Find matches between the two files
        matches = []
        for iso1, juncs1 in junctions1.items():
            for iso2, juncs2 in junctions2.items():
                # Check if block count is 2 (two-exon isoforms)
                if juncs1[1] == 2 and juncs2[1] == 2:
                    # Check chromosome, block count, first exon end, and last exon start
                    if juncs1[0] == juncs2[0] and juncs1[5] == juncs2[5] and juncs1[6] == juncs2[6]:
                        start_diff = abs(juncs1[3] - juncs2[3])
                        end_diff = abs(juncs1[4] - juncs2[4])
                        if start_diff <= max_diff and end_diff <= max_diff:
                            matches.append((iso1, iso2))
                            out.write(f"{iso1}\t{iso2}\n")
                else:
                    # For isoforms with more than 2 exons, check additional criteria
                    if juncs1[:3] == juncs2[:3]:  # Compare chrom, exon count, and junctions
                        start_diff = abs(juncs1[3] - juncs2[3])
                        end_diff = abs(juncs1[4] - juncs2[4])
                        
                        # Check first exon end and last exon start
                        first_exon_match = juncs1[5] == juncs2[5]
                        last_exon_match = juncs1[6] == juncs2[6]
                        
                        # Only consider a match if all conditions are met
                        if start_diff <= max_diff and end_diff <= max_diff and first_exon_match and last_exon_match:
                            matches.append((iso1, iso2))
                            out.write(f"{iso1}\t{iso2}\n")
        return matches

# Example usage:
# file1 = "ont.bed12"
# file2 = "pacbio.bed12"
# outfile = "matching_isoforms_new.txt"
argparser = argparse.ArgumentParser()
argparser.add_argument("file1", help="First BED12 file")
argparser.add_argument("file2", help="Second BED12 file")
argparser.add_argument("outfile", help="Output file")
args = argparser.parse_args()
file1 = args.file1
file2 = args.file2
outfile = args.outfile
matching_isoforms = find_matching_isoforms(file1, file2, outfile, max_diff=500)
# print("Matching isoforms with identical internal splice junctions and matching first/last exon ends:")
# for iso1, iso2 in matching_isoforms:
#     print(f"{iso1} in file1 matches {iso2} in file2")
