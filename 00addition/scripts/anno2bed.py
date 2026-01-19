"""
功能：批量把当前目录下的 hg38_*.txt 注释文件（列顺序：type chr start end ...）转为 BED4 (chr start-1 end type)，输出到 data/bed_output。
用法：在包含 hg38_*.txt 的目录执行 python anno2bed.py；输出目录默认为 /home/caow/03mnv/analyse3/00addition/data/bed_output。
"""

import os
import glob

# 1. 获取所有以 hg38_ 开头的 txt 文件
files = glob.glob("hg38_*.txt")

# 创建输出目录
output_dir = "/home/caow/03mnv/analyse3/00addition/data/bed_output"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for file in files:
    # --- 文件名处理逻辑 ---
    # 1. 去掉路径，只留文件名 (以防万一)
    basename = os.path.basename(file)
    
    # 2. 去掉 .txt 后缀
    name_core = basename.replace(".txt", "")
    
    # 3. 去掉 hg38_ 前缀
    if name_core.startswith("hg38_"):
        name_core = name_core.replace("hg38_", "", 1) # 只替换第一个
        
    # 4. 去掉 _aj 后缀 (如果存在)
    if name_core.endswith("_aj"):
        name_core = name_core[:-3] # 去掉最后3个字符 (_aj)
        
    # 5. 拼接新文件名
    new_filename = name_core + ".bed"
    out_path = os.path.join(output_dir, new_filename)
    
    print(f"转换: {file} -> {out_path}")
    
    # --- 内容转换逻辑 (BED4 格式) ---
    with open(file, 'r') as f_in, open(out_path, 'w') as f_out:
        for line in f_in:
            line = line.strip()
            if not line: continue
            
            parts = line.split()
            
            # 确保至少有4列 (Col1=Type, Col2=Chr, Col3=Start, Col4=End)
            if len(parts) >= 4:
                type_name = parts[0]
                chrom = parts[1]
                start = parts[2]
                end = parts[3]
                
                try:
                    # 坐标转换: Start - 1
                    start_0based = int(start) - 1
                    
                    # 写入: Chr  Start-1  End  Type
                    f_out.write(f"{chrom}\t{start_0based}\t{end}\t{type_name}\n")
                    
                except ValueError:
                    print(f"  [跳过非数字行] {line}")

print(f"\n全部完成！请查看文件夹: {output_dir}")