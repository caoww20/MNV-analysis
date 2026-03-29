import sys
input_file = sys.argv[1] # 
output_file = sys.argv[2]

with open(input_file) as f, open(output_file,'w') as w:
    for line in f:
        i = line.split('\t')[3].split(',')
        if len(i) >5:
            continue
        w.write(line)