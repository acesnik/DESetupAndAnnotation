import sys
import numpy as np

USAGE = "python make_rsem_edgeR_dataframe.py <outFile> file1 file2 ... fileN"
if len(sys.argv) < 3: print USAGE; exit()

outf, files = sys.argv[1], sys.argv[2:]
line_ct = sum(1 for line in open(files[0]))

all_transcripts = set()
cufflinks_results = {} # transcript_id : [expt_coverage1, expt_coverage2, ..., expt_coverageN]
for i, f in enumerate(files):
    print "processing file " + f
    transcript_ids = np.array([line.split('\t')[0] for j, line in enumerate(open(f)) if j != 0])
    coverage = np.array([line.split('\t')[8] for j, line in enumerate(open(f)) if j != 0])
    for j, x in enumerate(transcript_ids):
        if x in cufflinks_results: cufflinks_results[x][i] = coverage[j]
        else: cufflinks_results[x] = [0 if k != i else coverage[j] for k, f in enumerate(files)]

outf = open(outf, 'w')
outf.write("transfrag_id" + '\t' + '\t'.join("expt_" + str(i) for i, f in enumerate(files)) + '\n')
for transfrag_id in cufflinks_results:
    outf.write(transfrag_id + '\t' + '\t'.join(str(x) for x in cufflinks_results[transfrag_id]) + '\n')
outf.close()