import os
import sys

files = sys.argv[1].split(',')
output = open(sys.argv[2], 'w')

sources = [os.path.basename(x.split('-')[0]) for x in files]

header_written = False
for i in range(0, len(sources)):
    source = sources[i]
    input_file = files[i]

    f = open(input_file)
    if not header_written:
        header = '\t'.join(f.readline().strip().split('\t') + ['source'])
        output.write(header + '\n')
        header_written = True
    else:
        f.readline()

    for line in f.readlines():
        output.write('\t'.join(line.strip().split('\t') + [source]) + '\n')

output.close()
