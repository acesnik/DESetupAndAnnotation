import sys
import numpy as np

USAGE = "python make_rsem_dataframe.py <column_index> <outFile> file1 file2 ... fileN"
if len(sys.argv) < 4: print USAGE; exit()
column = int(sys.argv[1])
outf, files = sys.argv[2], sys.argv[3:]
line_ct = sum(1 for line in open(files[0]))

transcript_ids = np.array([line.split('\t')[0] for line in open(files[0])])
columns = [transcript_ids]
for file in files:
    columns.append(np.array([line.split('\t')[column] for line in open(file)]))
dataframe = np.column_stack(columns)
dataframe[0, 1:] = ["expt_" + str(i) for i in range(1, len(files) + 1)]

np.savetxt(outf, dataframe, delimiter="\t", fmt="%s")
