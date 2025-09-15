# 打开输入文件a.txt和输出文件b.txt
with open('1000G_gf_phase_fix.vcf', 'r') as input_file, open('1000G_gf_phase_fix_filter.vcf', 'w') as output_file:
    for line in input_file:
        if '#' in line:  # 如果行包含 #
             output_file.write(line)  
        else:
            columns = set(line.strip().split('\t')[9:])  # 按空格分割行
            if len(columns) > 1:  # 如果行数大于1
                output_file.write(line)  # 写入输出文件
            elif len(columns) == 1 and '0|0' not in columns:
                output_file.write(line)
