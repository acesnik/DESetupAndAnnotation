# Program for annotating RSEM-1.2.25 isoform results on a cufflinks assembly

__author__ = "anthonycesnik"
__date__ = "$July 18, 2016$"

import os.path
import optparse
from lxml import etree as et

class Exon:
    def __init__(self, chrom, start, end, strand, transcript_id, exon_no, attributes):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.transcript_id = transcript_id
        self.exon_no = exon_no
        self.attributes = attributes

    def __str__(self):
        return str(self.transcript_id) + "." + str(self.exon_no) + \
            "|".join(str(x) for x in [self.chrom, self.start, self.end, self.strand, self.info])


class Transcript:
    def __init__(self, attributes):
        self.chrom, self.start, self.end, self.strand = "", -1, -1, ""
        self.transcript_id = ""
        self.protein_seq = ""
        self.exons = []
        self.attributes = {}
        self.add_attributes(attributes)

    # get the chrom segment info, assuming exons has been initialized
    def get_chrom_segment_info(self):
        id_set = set([exon.transcript_id for exon in self.exons])
        chrom_set = set([exon.chrom for exon in self.exons])
        if not self.exons:
            print "Error processing empty exons for a transcript."
        elif len(id_set) == 1 and len(chrom_set) == 1:
            self.start = min([exon.start for exon in self.exons])
            self.end = max([exon.end for exon in self.exons])
            self.chrom = next(iter(chrom_set))
            self.transcript_id = next(iter(id_set))
        else:
            print "Error processing exons " + ",".join(str(exon) for exon in self.exons) + \
                  ". Exons must all belong to the same transcript to construct a Transcript object."
        
    # Remove attributes that don't match; add the rest
    def add_attributes(self, attributes):
        attributes = { key : item for key, item in attributes.items() if key not in self.attributes or item == self.attributes[key] }
        for key, item in attributes.items():
            self.attributes[key] = item
        return self

    def __str__(self):
        return ";".join(str(x) for x in [self.transcript_id, self.chrom, self.start, self.end, self.strand] + \
                        [str(exon) for exon in self.exons]) + '\t' + \
                ';'.join(str(key) + ":" + str(value) for key, value in self.attributes)


class GeneModelHandler:
    def __init__(self, gene_model, fasta, is_gtf):
        self.transcripts = {}
        curr_transcript_id = ""

        # Handle attributes. Add new transcript with these attributes if the transcript ID doesn't match the last.
        line_ct = sum(1 for line in open(gene_model))
        print "Processing " + gene_model + " (" + str(line_ct) + " lines)."
        for i, line in enumerate(open(gene_model)):
            if i % 1000 == 0: print "Processed " + str(i) + " lines out of " + str(line_ct) + " from " + gene_model + "."
            if not line.startswith('#'):
                line = line.strip().rstrip(';').split('\t')
                attribute_list = line[8].replace('; ', ';').split(';')
                new_attributes = self.get_attributes(attribute_list, is_gtf)
                transcript_id_attrib = 'transcript_id' if is_gtf else 'transcript'
                if transcript_id_attrib not in new_attributes:
                    continue
                elif new_attributes[transcript_id_attrib] == curr_transcript_id:
                    self.transcripts[curr_transcript_id].add_attributes(new_attributes)
                else:
                    if curr_transcript_id in self.transcripts: self.wrap_up_transcript(curr_transcript_id)
                    curr_transcript_id = new_attributes[transcript_id_attrib]
                    self.transcripts[curr_transcript_id] = Transcript(new_attributes)

                # Handle exon definitions
                if line[2] == 'exon':
                    chrom, start, end, strand = line[0], line[3], line[4], line[6]
                    transcript_id = new_attributes[transcript_id_attrib]
                    exon_number = new_attributes['exon_number'] if is_gtf else new_attributes['exon'].split('.')[1]
                    exon = Exon(chrom, start, end, strand, transcript_id, exon_number, new_attributes)
                    self.transcripts[curr_transcript_id].exons.append(exon)

        self.wrap_up_transcript(curr_transcript_id)

    def wrap_up_transcript(self, curr_transcript_id):
        self.transcripts[curr_transcript_id].get_chrom_segment_info()
        self.transcripts[curr_transcript_id].protein_seq = fasta.get_seq_by_id(curr_transcript_id)

    def get_attributes(self, attribute_list, is_gtf):
        attributes = {}
        for item in attribute_list:
                #GTF attributes
                if is_gtf:
                    item = item.split(' ')
                    attributes[item[0]] = item[1][1:-1]

                #GFF attributes
                else:
                    item = item.split('=')
                    if item[0] != 'info' and item[1].contains(':'):
                        gff_field = item[1].split(':')
                        attributes[gff_field[0]] = gff_field[1]
                    elif item[0] == 'info' and item[1].contains('InterPro'):
                        interpro_id = []
                        interpro_desc = []
                        for interpro in item[1][5:].split(' %0A'):
                            ipr_acce_idx = interpro.find('IPR')
                            ipr_desc_idx = interpro.find('description:') + 12
                            interpro_id.append(interpro[ipr_acce_idx: ipr_acce_idx + 9])
                            interpro_desc.append(interpro[ipr_desc_idx:])
                        attributes['interpro_id'] = interpro_id
                        attributes['interpro_desc'] = interpro_desc
        return attributes


# Reads the headers and sequences of a fasta file into RAM
class FastaManager:
    def __init__(self, protein_fasta):
        self.headers = []
        self.sequences = []
        line = protein_fasta.readline()
        while line.startswith ('#'):
            line = protein_fasta.readline()
        sequence = ""
        while line != "":
            if line.startswith(">"):
                accession = line.strip().split(' ')[0][1:]
                self.headers.append(accession) #modified for this specific case
                line = protein_fasta.readline()
                while not line.startswith(">"):
                    if line == "": break
                    sequence += line
                    line = protein_fasta.readline()
                self.sequences.append(sequence.replace('\n','').replace('\r',''))
                sequence = ""
        protein_fasta.close()
        
    def get_seq_by_id(self, id):
        for i, header in enumerate(fasta.headers):
                if header.find(id) >= 0:
                    return fasta.sequences[i]
        return ""


HTML_NS = "http://uniprot.org/uniprot"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
NAMESPACE_MAP = {None:HTML_NS, "xsi":XSI_NS}
UP = '{'+HTML_NS+'}'

def condense_xml_entry(entry):
    for element in entry:
        if element.tag not in [UP+'protein',UP+'accession',UP+'name',UP+'gene',UP+'organism',UP+'proteinExistence', \
                               UP+'depth',UP+'sequence',UP+'feature',UP+'dbReference']:
            entry.remove(element)
        elif element.tag == UP+'organism':
            for field in element:
                if field.tag != UP+'name': element.remove(field)
        else: continue
    return entry


# Parse Command Line
parser = optparse.OptionParser()
parser.add_option('-i', '--rsem_isoform_results', dest='rsem_isoform_results', help='RSEM version 1.2.25 isoform abundance results file.')
parser.add_option('-o', '--output_file', dest='output_file', help='Annotated datframe output.')
parser.add_option('-a', '--pep_all_fasta', dest='pep_all_fasta', help='Protein fasta, corresponding to the genome and gene model used to run cufflinks.')
parser.add_option('-c', '--gene_model_cuff', dest='gene_model_cuff', help='Cuffmerge GTF file used to run RSEM. Should removed zero-abundance transcripts beforehand and used the reference gene model to annotate.')
parser.add_option('-f', '--gene_model_gff', dest='gene_model_gff', default=None, help='GFF gene model file used to run cufflinks.')
parser.add_option('-g', '--gene_model_gtf', dest='gene_model_gtf', default=None, help='GTF gene model file used to run cufflinks.')
parser.add_option('-x', '--uniprot_xml', dest='uniprot_xml', help='UniProt-XML for the organism.')
(options, args) = parser.parse_args()

### FILE SETUP
expression_data = os.path.abspath(options.rsem_isoform_results)
out = os.path.abspath(options.output_file)
pep_fa = os.path.abspath(options.pep_all_fasta)
uniprot_xml = os.path.abspath(options.uniprot_xml)
cuff = os.path.abspath(options.gene_model_cuff)
if options.gene_model_gff and options.gene_model_gtf:
    print "Both GFF and GTF were specified. Please choose one or the other."
    exit(2)
elif options.gene_model_gff: gene_model_ref = os.path.abspath(options.gene_model_gff)
elif options.gene_model_gtf: gene_model_ref = os.path.abspath(options.gene_model_gtf)
else: print "Something's wrong..."

# pep_all_fasta from reference
fasta = FastaManager(open(pep_fa))
seqs = set(fasta.sequences)
print "Processed " + str(len(seqs)) + " protein sequences from the protein fasta."

# cuffmerge assembly
cufflinks_assembly = GeneModelHandler(cuff, fasta, not options.gene_model_gff)
print "Processed " + str(len(cufflinks_assembly.transcripts)) + " transcripts from the cufflinks assembly."

# gene model transcripts
ref_transcripts = GeneModelHandler(gene_model_ref, fasta, not options.gene_model_gff)
aggregate_attributes = sorted(list(set([transcript.attributes for transcript in ref_transcripts.transcripts] + \
                           [transcript.attributes for transcript in cufflinks_assembly.transcripts])))
print "Processed " + str(len(ref_transcripts.transcripts)) + " transcripts from the " + \
      ("GTF" if options.gene_model_gtf else "GFF") + \
      " gene model."

# uniprot_xml
r = et.Element(UP+'uniprot', nsmap=NAMESPACE_MAP)
db = et.ElementTree(r)
xml_iterator = et.iterparse(open(uniprot_xml, 'r'), remove_blank_text=True)
for event, el in xml_iterator:
    if el.tag == UP+"entry": 
        r.append(condense_xml_entry(el))

dbRef_interpro = [e for e in r.findall('.//'+UP+'dbReference') if e.get('type') == 'InterPro']
xml_seqs = [e for e in r.findall('.//'+UP+'sequence') if e.text is not None]
print "Added " + str(len(xml_seqs)) + " sequences from the UniProt-XML."

### ANNOTATION
# iterate through all lines of the expression data and annotate with the information in the references
line_ct = sum(1 for line in open(expression_data))
transcript_id_attrib = 'transcript_id' if options.gene_model_gtf else 'transcript'
header = ''
for i, line in enumerate(open(expression_data)):
    if i % 1000 == 0: print "Processed " + str(i) + " lines out of " + str(line_ct) + "."
    line = line.strip()
    
    if i == 0: # handle the header of the expression data
        header = line
        out.write('transfrag_id\t' + header + '\t' + 'transfrag_specifications\ttransfrag_attributes\tref_transcript_specs\tref_transcript_attributes\tUniprot_accessions\tUniProt_names\tUniprot_relatedInterProAccessions\tGO-id;GO-term;GO-evidence;GO-project\n')
        continue
    transcript_id = line.split('\t')[0]
    transfrag = cufflinks_assembly.transcripts[transcript_id]
    nearest_ref = transfrag.attributes["nearest_ref"] if "nearest_ref" in transfrag.attributes else None
    if not nearest_ref:
        out.write('\t'.join([line, str(transfrag)]) + '\n')
        continue
    ref_transcript = ref_transcripts[nearest_ref]

    # get the uniprot entries with the same sequence, if any
    uniprot_entries = [seq.getparent() for seq in xml_seqs if seq.text.replace('\n','').replace('\r','') == ref_transcript.protein_seq]
    if len(uniprot_entries) == 0:
        out.write('\t'.join([line, str(transfrag), str(ref_transcript)]) + '\n')
        continue

    uniprot_accessions = ';'.join(entry.find(UP+'accession').text for entry in uniprot_entries)
    uniprot_names = ';'.join(entry.find('.//'+UP+'fullName').text for entry in uniprot_entries)
    related_dbRefs = []
    for entry in uniprot_entries:
        for dbRef in entry.findall('.//'+UP+'dbReference'):
            related_dbRefs.append(dbRef)

    # get the related interpro terms
    related_interpro_dbRefs = [e for e in related_dbRefs if e.get('type') == 'InterPro']
    related_interpro_ids = ';'.join([','.join(e.get('id') for e in related_interpro_dbRefs)])

    # get the related go terms
    go_terms = [e for e in related_dbRefs if e.get('type') == 'GO']
    go_properties = { term.get('id') : term.findall('.//'+UP+'property') for term in go_terms }
    go_str = '\t'.join(term + ';' + ';'.join(p.get('value') for p in go_properties[term]) for term in go_properties)

    out.write('\t'.join([line, str(transfrag), str(ref_transcript), uniprot_accessions, uniprot_names, related_interpro_ids, go_str]) + '\n')
out.close()