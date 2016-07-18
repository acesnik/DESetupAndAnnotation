# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "anthonycesnik"
__date__ = "$Dec 7, 2015 10:30:36 PM$"

import sys
import os.path
from subprocess import check_output
import optparse

#def make_transcript_gtf(gtf_path):
#    transcript_gtf_txt = check_output(['grep', '\ttranscript\t', gtf_path])
#    trascript_gtf_path = os.path.abspath(gtf_path[:-3] + 'transcript.gtf')
#    outF = open(trascript_gtf_path, 'w')
#    outF.write(transcript_gtf_txt)
#    outF.close()
#    return trascript_gtf_path
    
def get_transcript_annotation(transcript_line):
    line = transcript_line.split('\t')
    chrom_start, chrom_end, chromosome, strand = line[3], line[4], line[0], line[6]
    attribute_list = line[8].split('; ')
    attributes = {} #This is what's consistent across Ensembl gene models...
    for item in attribute_list:
        item = item.split(' ')
        attributes[item[0]] = item[1][1:-1]
        attributes['chrom_start'] = chrom_start
        attributes['chrom_end'] = chrom_end
        attributes['chromosome'] = chromosome
        attributes['strand'] = strand
    return attributes

#Parse Command Line
parser = optparse.OptionParser()
parser.add_option( '-i', '--rsem_isoform_results', dest='rsem_isoform_results', help='RSEM version 1.2.25 isoform abundance results file.' )
parser.add_option( '-g', '--gene_model', dest='gene_model', default=None, help='GTF gene model file used to run RSEM. If not already only transcripts, a transcript.gtf file will be generated. Used to annotate Ensembl accessions.')
(options, args) = parser.parse_args()
    
try:
    #Prepare transcript-level gene mode
    geneModelPath = os.path.abspath(options.gene_model)
    geneModelFile = open(geneModelPath, 'r')
    geneModel = {}
    for line in geneModelFile:
        if not line.startswith('#') and line.split('\t')[2] == 'transcript': 
            attributes = get_transcript_annotation(line)
            geneModel[attributes['transcript_id']] = attributes
    geneModelFile.close()
except Exception, e:
    print >> sys.stderr, "Parsing gene model and creating transcript_gtf failed: %s" % e
    exit(2)

#Annotate results
try:
    rsem_results, rsem_annotated_results = os.path.abspath(options.rsem_isoform_results), os.path.abspath(options.rsem_isoform_results + '.annotated')
    linect = sum(1 for line in open(rsem_results))
    rsem_results, rsem_annotated_results = open(rsem_results, 'r'), open(rsem_annotated_results, 'w')
    print "Using transcripts listed at " + geneModelPath + " to annotate RSEM results."
    for i, line in enumerate(rsem_results):
        if i % 10000 == 0: print "RSEM results line " + str(i) + " of " + str(linect)
        if line.startswith('transcript_id'): 
            rsem_annotated_results.write(line.strip() + "\tchromosome\tstrand\tchrom_start\tchrom_end\tgene_name\tgene_biotype\ttranscript_name\ttranscript_biotype\n")
            continue
        transcript_id, gene_id, length, effective_length, expected_count, tpm, fpkm, isoform_percentage = line.strip().split('\t')
        annotations = geneModel[transcript_id]
        rsem_annotated_results.write(line.strip() + '\t' + annotations['chromosome'] + '\t' + annotations['strand'] + '\t' + annotations['chrom_start'] + '\t' + annotations['chrom_end'] + '\t' + annotations['gene_name'] + '\t' + annotations['gene_biotype'] + '\t' + annotations['transcript_name'] + '\t' + annotations['transcript_biotype'] + '\n')
    rsem_results.close()
    rsem_annotated_results.close()
        
except Exception, e:
    print >> sys.stderr, "Parsing RSEM-1.2.25 isoform results file failed: %s" % e