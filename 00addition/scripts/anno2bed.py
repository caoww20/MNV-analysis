"""
Purpose: batch-convert hg38_*.txt annotation files (columns: type chr start end ...)
to BED4 (chr start-1 end type) and write to data/bed_output.
Usage: run python anno2bed.py in a directory containing hg38_*.txt;
default output dir: /home/caow/03mnv/analyse3/00addition/data/bed_output.
"""

import os
import glob

# 1. Get all txt files starting with hg38_
files = glob.glob("hg38_*.txt")

# Create output directory
output_dir = "/home/caow/03mnv/analyse3/00addition/data/bed_output"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for file in files:
    # --- File name processing ---
    # 1. Strip path, keep filename only
    basename = os.path.basename(file)
    
    # 2. Remove .txt suffix
    name_core = basename.replace(".txt", "")
    
    # 3. Remove hg38_ prefix
    if name_core.startswith("hg38_"):
        name_core = name_core.replace("hg38_", "", 1) # Replace first only
        
    # 4. Remove _aj suffix if present
    if name_core.endswith("_aj"):
        name_core = name_core[:-3] # Drop last 3 chars (_aj)
        
    # 5. Build new filename
    new_filename = name_core + ".bed"
    out_path = os.path.join(output_dir, new_filename)
    
    print(f"Convert: {file} -> {out_path}")
    
    # --- Content conversion (BED4) ---
    with open(file, 'r') as f_in, open(out_path, 'w') as f_out:
        for line in f_in:
            line = line.strip()
            if not line: continue
            
            parts = line.split()
            
            # Ensure at least 4 columns (Type, Chr, Start, End)
            if len(parts) >= 4:
                type_name = parts[0]
                chrom = parts[1]
                start = parts[2]
                end = parts[3]
                
                try:
                    # Coordinate conversion: Start - 1
                    start_0based = int(start) - 1
                    
                    # Write: Chr  Start-1  End  Type
                    f_out.write(f"{chrom}\t{start_0based}\t{end}\t{type_name}\n")
                    
                except ValueError:
                    print(f"  [Skip non-numeric line] {line}")

print(f"\nAll done! Output folder: {output_dir}")