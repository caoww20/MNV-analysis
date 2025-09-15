import sys
input1=sys.argv[1]
input2=sys.argv[2]
output=sys.argv[3]


with open(input1, 'r') as file_a, open(input2, 'r') as file_b, open(output, 'w') as output_file:
    a_lines = set(file_a.readlines())
    b_lines = set(file_b.readlines())

    difference = a_lines.difference(b_lines)

    output_file.writelines(difference)