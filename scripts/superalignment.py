#!/usr/bin/python
import sys, os
folder = sys.argv[1]
data = {}
key = None
for f in os.listdir(folder):
    if f.endswith('trimmed'):
        for line in open('%s%s' % (folder, f), 'r'):
            if line.startswith('>'):
                key=line.strip()[1:]
                if key in data:
                    pass
                else:
                    data.update({key: ''})	
            else:
                data[key] += line.strip()
# output
out = open('%s/supermatrix.aln.faa' % folder, 'w')
for key in data:
    out.write('>%s\n' % key)
    out.write('%s\n' % data[key]) 
out.close()

