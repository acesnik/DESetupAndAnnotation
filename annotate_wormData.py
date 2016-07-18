import sys
from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'


# Reads the headers and sequences of a fasta file into RAM
def read_protein_fasta(protein_fasta):
    proteinFasta = ([],[]) #headers, #sequences
    line = protein_fasta.readline()
    while line.startswith ('#'):
        line = protein_fasta.readline()
    sequence = ""
    while line != "":
        if line.startswith(">"):
            accession = line.strip().split(' ')[0][1:]
            proteinFasta[0].append(accession) #modified for this specific case
            line = protein_fasta.readline()
            while not line.startswith(">"):
                if line == "": break
                sequence += line
                line = protein_fasta.readline()
            proteinFasta[1].append(sequence.replace('\n','').replace('\r',''))
            sequence = ""
    protein_fasta.close()
    return proteinFasta


# Structure for the WormBase GFF entry
class GffTranscript:
    def __init__(self, line, worm_pep_fa):
        line = line.strip()
        fields = line.split('\t')
        info_fields = fields[8].split(';')
        self.gff_id = info_fields[0].split('=')[1]
        self.protein_seq = worm_pep_fa[1][worm_pep_fa[0].index(self.gff_id[11:])]
        self.interpro_id = []
        self.interpro_desc = []
        if len(info_fields) > 3:
            for interpro in info_fields[3][5:].split(' %0A'):
                ipr_acce_idx = interpro.find('IPR')
                ipr_desc_idx = interpro.find('description:') + 12
                self.interpro_id.append(interpro[ipr_acce_idx: ipr_acce_idx + 9])
                self.interpro_desc.append(interpro[ipr_desc_idx:])

    def __str__(self):
        return '\t'.join([';'.join(self.interpro_id), ';'.join(self.interpro_desc)])


USAGE = "python make_rsem_edgeR_dataframe.py <outFile> <expression_data.txt> <worm.gff / ensembl.gtf> <uniprot.xml> <worm.pep.fa / ensembl.pep.fa>"
if len(sys.argv) != 6: print USAGE; exit()
outfn, expression_data, worm_gff, worm_xml, worm_pep_fa = sys.argv[1:]
out = open(outfn, 'w')

# worm_pep_fa
worm_fasta = read_protein_fasta(open(worm_pep_fa))
worm_seqs = set(worm_pep_fa[1])

# worm_gff transcripts only
worm_gff_mRNA = filter(None, [line if line.find('mRNA') >= 0 else None for line in open(worm_gff)])
worm_transcripts = { GffTranscript(t, worm_fasta).gff_id : GffTranscript(t, worm_fasta) for t in worm_gff_mRNA }

# uniprot_xml
db = et.parse(worm_xml)
r = db.getroot()
r.remove(r.find(UP+'copyright'))
dbRef_interpro = filter(None, [e if e.get('type') == 'InterPro' else None for e in r.findall('.//'+UP+'dbReference')])
xml_seqs = [e for e in r.findall('.//'+UP+'sequence') if e.text is not None]
print len(xml_seqs)

# iterate through all lines of the expression data and annotate with the information in the references
line_ct = sum(1 for line in open(expression_data))
header = ''
for i, line in enumerate(open(expression_data)):
    if i % 1000 == 0: print "processed " + str(i) + " lines out of " + str(line_ct)
    line = line.strip()

    # handle the header of the expression data
    if i == 0:
        header = line
        out.write('transcript_id\t' + header + '\t' + 'InterPro_accession\tInterPro_description\tUniprot_accessions\tUniProt_names\tUniprot_relatedInterProAccessions\tGO_id1\tGO_term1\tGO_evidence1\GO_project1\tGO_id2 ... etc\n')
        continue
    transcript_id = line.split('\t')[0]
    gff_transcript = worm_transcripts[transcript_id]

    # get the uniprot entries with the same sequence, if any
    uniprot_entries = filter(None, [seq.getparent() if seq.text.replace('\n','').replace('\r','') == gff_transcript.protein_seq else None for seq in xml_seqs])
    # print "num of uniprot entries:" + str(len(uniprot_entries)) + " out of " + str(len(r))
    if len(uniprot_entries) == 0: out.write('\t'.join([line, str(gff_transcript), '', '', '', '']) + '\n'); continue
    uniprot_accessions = ';'.join(entry.find(UP+'accession').text for entry in uniprot_entries)
    uniprot_names = ';'.join(entry.find('.//'+UP+'fullName').text for entry in uniprot_entries)

    related_dbRefs = []
    for entry in uniprot_entries:
        for dbRef in entry.findall('.//'+UP+'dbReference'):
            related_dbRefs.append(dbRef)

    # get the related interpro terms
    related_interpro_dbRefs = [e for e in related_dbRefs if e.get('type') == 'InterPro']
    # print "num of interpro refs:" + str(len(related_interpro_dbRefs))
    related_interpro_ids = ';'.join([','.join(e.get('id') for e in related_interpro_dbRefs)])

    # get the related go terms
    go_terms = [e for e in related_dbRefs if e.get('type') == 'GO']
    go_properties = { term.get('id') : term.findall('.//'+UP+'property') for term in go_terms }
    go_str = '\t'.join(term + '\t' + '\t'.join(p.get('value') for p in go_properties[term]) for term in go_properties)

    out.write('\t'.join([line, str(gff_transcript), uniprot_accessions, uniprot_names, related_interpro_ids, go_str]) + '\n')

out.close()