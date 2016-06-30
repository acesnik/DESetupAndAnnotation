import sys
from lxml import etree as et

HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

class InterProTranscript():
    def __init__(self, line):
        fields = line.strip().split('\t')[8].split(';')
        self.gff_id = fields[0].split('=')[1]
        self.interpro_id = []
        self.interpro_desc = []
        if len(fields) > 3:
            for interpro in fields[3][5:].split(' %0A'):
                ipr_acce_idx = interpro.find('IPR')
                ipr_desc_idx = interpro.find('description:') + 12
                self.interpro_id.append(interpro[ipr_acce_idx: ipr_acce_idx + 9])
                self.interpro_desc.append(interpro[ipr_desc_idx:])

    def __str__(self):
        return '\t'.join([';'.join(self.interpro_id), ';'.join(self.interpro_desc)])


USAGE = "python make_rsem_edgeR_dataframe.py <outFile> <expression_data.txt> <worm.gff> <worm_uniprot.xml>"
if len(sys.argv) != 5: print USAGE; exit()

outfn, expression_data, worm_gff, worm_xml = sys.argv[1:]
out = open(outfn, 'w')
worm_transcripts = filter(None, [line if line.find('mRNA') >= 0 else None for line in open(worm_gff)])
interpro_transcripts = { InterProTranscript(t).gff_id : InterProTranscript(t) if t.find('info') >= 0 else None for t in worm_transcripts }
db = et.parse(worm_xml)
r = db.getroot()
dbRef_interpro = filter(None, [e if e.get('type') == 'InterPro' else None for e in r.findall('.//'+UP+'dbReference')])

line_ct = sum(1 for line in open(expression_data))
header = ''
for i, line in enumerate(open(expression_data)):
    if i % 1000 == 0: print "processed " + str(i) + " lines out of " + str(line_ct)
    line = line.strip()
    if i == 0:
        header = line
        out.write('transcript_id\t' + header + '\t' + 'InterPro_accession\tInterPro_description\tUniprot_accessions\tUniProt_names\tUniprot_relatedInterProAccessions\tGO_id1\tGO_term1\tGO_evidence1\GO_project1\tGO_id2 ... etc\n')
        continue
    transcript_id = line.split('\t')[0]
    interpro_proteins = interpro_transcripts[transcript_id]
    if interpro_proteins == None: out.write(line + '\n'); continue
    interpro_proteins_dbRef_elements = filter(None, [e if e.get('id') in interpro_proteins.interpro_id else None for e in dbRef_interpro])
    print "num of interpro refs:" + str(len(interpro_proteins_dbRef_elements))
    uniprot_entries = set([dbRef.getparent() for dbRef in interpro_proteins_dbRef_elements])
    print "num of uniprot entries:" + str(len(uniprot_entries)) + " out of " + str(len(r))
    uniprot_accessions = ';'.join(entry.find(UP+'accession').text for entry in uniprot_entries)
    uniprot_names = ';'.join(entry.find('.//'+UP+'fullName').text for entry in uniprot_entries)
    related_dbRefs = []
    for entry in uniprot_entries:
        for dbRef in entry.findall('.//'+UP+'dbReference'):
            related_dbRefs.append(dbRef)

    related_interpro_dbRefs = filter(None, [e if e.get('type') == 'InterPro' else None for e in related_dbRefs])
    related_interpro_ids = ';'.join([','.join(e.get('id') for e in related_interpro_dbRefs)])

    go_terms = filter(None, [e if e.get('type') == 'GO' else None for e in related_dbRefs])
    go_properties = { term.get('id') : term.findall('.//'+UP+'property') for term in go_terms }
    go_str = '\t'.join(term + '\t' + '\t'.join(p.get('value') for p in go_properties[term]) for term in go_properties)

    out.write('\t'.join([line, str(interpro_proteins), uniprot_accessions, uniprot_names, related_interpro_ids, go_str]) + '\n')

out.close()