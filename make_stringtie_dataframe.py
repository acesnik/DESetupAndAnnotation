import sys
import numpy as np

USAGE = "python make_stringtie_dataframe.py <outFile> file1 file2 ... fileN"
if len(sys.argv) < 4: print USAGE; exit()
outf, files = sys.argv[1], sys.argv[2:]
line_ct = sum(1 for line in open(files[0]))

transcript_ids = np.array([[val for val in line.split('\t')[8].split("; ") if val.startswith("transcript_id")][0][15:].replace(';', '').replace('"','') for line in open(files[0]) if not line.startswith("#") and line.split('\t')[2] == "transcript"])
columns = [transcript_ids]
for file in files:
    columns.append(np.array([[val for val in line.split('\t')[8].split("; ") if val.startswith("cov")][0][5:].replace(';', '').replace('"','') for line in open(file) if not line.startswith("#") and line.split('\t')[2] == "transcript"]))
dataframe = np.column_stack(columns)
dataframe[0, 1:] = ["expt_" + str(i) for i in range(1, len(files) + 1)]

np.savetxt(outf, dataframe, delimiter="\t", fmt="%s")
