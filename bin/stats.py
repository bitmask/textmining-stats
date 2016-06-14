import argparse
from collections import defaultdict
import re


################################################################################
# This code is written by Helen Cook, with contributions from ...
# 
# Distributed under the BSD license. 
################################################################################

def parse_annotations(inputfile, annotations):
    # Input format is one annotation per line
    with open(inputfile, 'r') as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 5:
                (annotid, doctype, docid, annottype, annotdict, normid, blank, user, unused, text, unused2, start, end) = cols
                if annotdict == "NCBITax" or annotdict == "Uniprot":
                    annot = {'text': text, 'norm': normid, 'start': start, 'end': end, 'user': user, 'docid': docid, 'annottype': annottype}
                    annotations[docid][user].append(annot)
    return annotations


def same_boundaries(annot, annot_list):
    # check whether an annotation with the same boundaries as annot exists on the annot_list
    for an in annot_list:
        if an['end'] == annot['end'] and an['start'] == annot['start']:
            return an
    return None


def same_normalization(annot, an, uniprot):
    if an['norm'] == annot['norm']:
        return True
    else:
        if annot['norm'] == '???' or an['norm'] == '???' or annot['norm'] == 'UNKN' or an['norm'] == 'UNKN':
            # if a normalization was not entered by the user, then say they are not the same
            return False
        if an['annottype'] == "e_2":
            if blast(uniprot[an['norm']]['sequence'], uniprot[annot['norm']]['sequence']):
                return True
    return False


def blast(seq1, seq2):
    # run blast on these sequences and if they have above 90% identity on 90% of the sequence, then consider them to be the same
    # TODO
    return False


def inter_annotator(annotations, n_annotators, uniprot):
    # Print interannotator agreement for all documens 
    for document in annotations.keys():
        if len(annotations[document].keys()) == n_annotators:
            print document
            for user1 in annotations[document].keys():
                for user2 in annotations[document].keys():
                    if user1 == user2:
                        continue
                    for annot in annotations[document][user1]:
                        #print annot
                        match = same_boundaries(annot, annotations[document][user2])
                        if match:
                            if same_normalization(annot, match, uniprot):
                                # print for debugging
                                print "agree " + user1 + " " + user2 + " " + annot['text'] + " " + annot['start'] + " " + annot['end']
                            else:
                                print "norm diff " + user1 + " " + user2 + " " + annot['text'] + " " + annot['start'] + " " + annot['end'] + " " + annot['norm'] + " " + match['norm']
                        else:
                            print "no match for " + user1 + " " + annot['text'] + " " + annot['start'] + " " + annot['end'] + " " + annot['norm']

    # in above code need to calculate precision and recall of each user against each other user, and populate stats
    # TODO
    stats = {} # dictionary user -> user -> recall/precision
    return stats;


def print_stats(stats):
    for user1 in stats:
        for user2 in stats[user1]:
            recall = stats[user1][user2]['recall']
            precision = stats[user1][user2]['precision']
            fscore = 2 * recall * precision / (precision + recall)
            print user1 + "\t" + user2 + "\t" + recall + "\t" + precision + "\t" + fscore
    return True


def parse_uniprot_promeome_dir(dirname):
    # relies on having only one isolate for each species
    import glob
    allproteins = {}
    for xml in glob.glob(dirname + "/*.xml"):
        result = parse_uniprot_single_xml(xml)
        for entry_name, this_protein in result.iteritems():
            # merge the dictionaries
            allproteins[entry_name] = this_protein
    return allproteins


def parse_uniprot_single_xml(filename):
    '''
    Read one inidividual xml files (representing one protein) that has been downloaded to a dir, eg as part of a proteome
    '''

    #sys.stderr.write("parsing " + filename + "\n")
    
    import xml.etree.cElementTree as ET
    try:
        tree = ET.parse(filename)
    except:
        raise

    root = tree.getroot()
    allproteins = {}

    for entry in root:
        source = "" # trEMBL or swiss prot
        entry_name = ""
        short_name = ""
        entry_ids = []
        isolate_name = ""
        isolate_name_details = ""
        isolate_abbrev = ""
        isolate_taxid = 0
        protein_name = ""
        aliases = []
        component_aliases = {}
        gene_name = ""
        date = ""
        sequence = ""
        hosts = [] # unused
        hosts_taxids = []
        chains = []
        pfam = {}
        go = []
        complete_prot = False
        reference_prot = False
        polyprotein = False
        species = 0
        string_id = ""
        embl_ids = []
        gene_ids = []


        if "entry" in entry.tag:
            date = entry.attrib['created']
            source = entry.attrib['dataset'] 
        else:
            continue
        for child in entry:
            if "name" in child.tag:
                entry_name = child.text
            if "accession" in child.tag:
                entry_ids.append(child.text)
            if "organism" in child.tag and "Host" not in child.tag:
                for taxentry in child:
                    if "name" in taxentry.tag:
                        if taxentry.attrib['type'] == "scientific":
                            isolate_name = taxentry.text
                        if taxentry.attrib['type'] == "common":
                            isolate_name_details = taxentry.text
                        if taxentry.attrib['type'] == "synonym":
                            isolate_abbrev = taxentry.text
                    if "dbReference" in taxentry.tag:
                        if taxentry.attrib['type'] == "NCBI Taxonomy":
                            isolate_taxid = int(taxentry.attrib['id'])
            if "protein" in child.tag:
                for rname in child:
                    for name in rname:
                        if "fullName" in name.tag:
                            if 'recommendedName' in rname.tag:
                                protein_name = name.text
                                if re.search('olyprotein', protein_name):
                                    polyprotein = True
                            else:
                                aliases.append(name.text)
                                pass
                        if 'shortName' in name.tag:
                            aliases.append(name.text)
                            pass

                    # get the recommened and alternate names for the chains
                    if 'component' in rname.tag:
                        rec = ""
                        for r in rname:
                            for n in r:
                                if 'recommendedName' in r.tag:  # this logic is right, it is complicated
                                    if 'fullName' in n.tag:
                                        rec = n.text
                                    if 'shortName' in n.tag:
                                        if rec in component_aliases:
                                            component_aliases[rec].append(n.text)
                                        else:
                                            component_aliases[rec] = [n.text]
                                if 'alternativeName' in r.tag:
                                    if 'Name' in n.tag:
                                        if rec in component_aliases:
                                            component_aliases[rec].append(n.text)
                                        else:
                                            component_aliases[rec] = [n.text]
            if "gene" in child.tag:
                for name in child:
                    gene_name = name.text
            if "sequence" in child.tag:
                sequence = child.text.replace('\n', '')
            if "organismHost" in child.tag:
                for ref in child:
                    if "dbReference" in ref.tag:
                        hosts_taxids.append(ref.attrib['id'])
            if "feature" in child.tag:
                if child.attrib['type'] == "chain":
                    chain = {}
                    chain['id'] = child.attrib['id']
                    if 'status' not in child.attrib: # some chains have status='potential' and these seem to be crap
                        chain['name'] = child.attrib['description']
                        if chain['name'] in component_aliases:
                            chain['aliases'] = list(set(component_aliases[chain['name']]))
                        else:
                            chain['aliases'] = []
                        complete = True
                        for loc in child:
                            for pos in loc:
                                if "begin" in pos.tag:
                                    if "status" in pos.attrib and pos.attrib['status'] == "unknown":
                                        complete = False
                                    else:
                                        chain['start'] = pos.attrib['position']
                                if "end" in pos.tag:
                                    if 'status' in pos.attrib and pos.attrib['status'] == 'unknown':
                                        complete = False
                                    else:
                                        chain['end'] = pos.attrib['position']
                        if complete:
                            chains.append(chain)
            if "dbReference" in child.tag:
                if child.attrib['type'] == "Pfam":
                    pfam[child.attrib['id']] = {} 
                if child.attrib['type'] == "GO":
                    go.append(child.attrib['id'])
                if child.attrib['type'] == "EMBL":
                    for p in child:
                        if p.attrib['type'] == 'protein sequence ID':
                            embl_ids.append(p.attrib['value'])
                if child.attrib['type'] == "RefSeq":
                    for p in child:
                        if 'type' in p.attrib:
                            if p.attrib['type'] == 'nucleotide sequence ID':
                                embl_ids.append(p.attrib['value'])
                if child.attrib['type'] == "GeneID":
                    gene_ids.append(child.attrib['id'])


            if "keyword" in child.tag:
                if child.attrib['id'] == 'KW-0181':
                    complete_prot = True
                if child.attrib['id'] == 'KW-1185':
                    reference_prot = True

        if isolate_name_details != "":
            isolate_name += " (" + isolate_name_details + ")"
        if isolate_abbrev == "":
            isolate_abbrev = isolate_name
        else:
            isolate_name += " (" + isolate_abbrev + ")"

        this_protein = {'entry_id': entry_ids[0], 'entry_name': entry_name, 'isolate_name': isolate_name, 'isolate_abbrev': isolate_abbrev, 'gene_name': gene_name, 'protein_name': protein_name, 'species': species, 'isolate_taxid': isolate_taxid, 'reference_proteome': reference_prot, 'complete_proteome': complete_prot, 'polyprotein': polyprotein, 'date': date, 'hosts_taxids': list(set(hosts_taxids)), 'sequence': sequence, 'chain': chains, 'pfam': pfam, 'go': go, 'aliases': list(set(aliases)), 'embl_ids': embl_ids, 'gene_ids': gene_ids, 'source': source}

        allproteins[entry_name] = this_protein        
        
    return allproteins


def main():
    parser = argparse.ArgumentParser(description="Calculate interannotator etc stats")
    # TODO arguments should be the same as the tagger
    parser.add_argument('-f', '--files',
                required=True,
                default = [],
                nargs='+',
                dest='inputfile',
                help="list of annotation files")

    parser.add_argument('-u', '--uniprot',
                required=True,
                dest='uniprot',
                help="uniprot xml dump")

    args=parser.parse_args()

    n_annotators = 0
    annotations = defaultdict(lambda: defaultdict(list))
    if args.inputfile:
        for inputfile in args.inputfile:
            # add all annotations into one dictionary
            annotations = parse_annotations(inputfile, annotations)
            n_annotators += 1

    uniprot = parse_uniprot_single_xml(args.uniprot)
    stats = inter_annotator(annotations, n_annotators, uniprot)
    print_stats(stats)


if __name__ == "__main__":
    main()
