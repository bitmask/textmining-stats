import argparse
from collections import defaultdict
import re
import Bio
import shlex
from subprocess import Popen, PIPE
import os.path
import ncbi_taxonomy.ncbi_taxonomy
import pprint
import numpy as np

################################################################################
# This code is written by Helen Cook, with contributions from ...
# 
# Distributed under the BSD license. 
################################################################################

def get_user_from_inputfile(inputfile):
    field = os.path.basename(inputfile).split(".")
    if field[1] == "test2":
        return "Helen"
    user = field[1].split("-")
    if user[0] == "Kerstin":
        return "Cristina"
    return user[0]

def parse_annotations(inputfile, annotations):
    # Input format is one annotation per line
    user = get_user_from_inputfile(inputfile)
    with open(inputfile, 'r') as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 5:
                (annotid, doctype, docid, annottype, annotdict, normid, blank, userX, unused, text, unused2, start, end) = cols
                if annotdict == "NCBITax" or annotdict == "Uniprot":
                    annot = {'text': text, 'norm': normid, 'start': start, 'end': end, 'user': user, 'docid': docid, 'annottype': annottype}
                    annotations[docid][user].append(annot)
    return annotations

def convert_annottype(annottype):
    if int(annottype) == -2 or int(annottype) == -3:
        return "e_1" # species
    else:
        return "e_2" # protein
    
def parse_tagger(tagger_output_file, entities_file, annotations):
    #Parse the tagger output file and also include the x_entities entries

    entities = {}
    with open(entities_file, "r") as f:
        for line in f:
            (serialno, annottype, norm) = line.rstrip("\n").split("\t")
            entities[serialno] = [annottype, norm]
            if serialno == 10531173 or serialno == 10531174:
                print "serialno: " + serialno + " " + annottype + " " + norm
                print stop

    user = 'tagger' #For identification from other users.
    with open(tagger_output_file, "r") as f:
        for line in f:
            (docid, unused, unused2, start, end, annot, taxid, serialno) = line.rstrip("\n").split("\t")
            (annottype, normid) = entities[serialno]
            annottype = convert_annottype(annottype)
            annot = {'text': annot, 'norm': normid, 'start': start, 'end': str(int(end)+1), 'user': user, 'docid': docid, 'annottype': annottype}
            annotations[docid][user].append(annot)
                    
    return annotations

def same_boundaries(annot, annot_list):
    # check whether an annotation with the same boundaries as annot exists on the annot_list
    for an in annot_list:
        if an['end'] == annot['end'] and an['start'] == annot['start']:
            return an
    return None


def same_normalization(annot, an, uniprot, taxtree):
    if an['norm'] == annot['norm']:
        return True
    else:
        if annot['norm'] == '???' or an['norm'] == '???' or annot['norm'] == 'UNKN' or an['norm'] == 'UNKN':
            # if a normalization was not entered by the user, then say they are not the same
            return False
        if an['annottype'] == "e_2":
            # for proteins, check that they are close enough to each other
            if an['norm'] in uniprot and annot['norm'] in uniprot:
                #print an['norm'] + " " + annot['norm']
                if blast(uniprot[an['norm']]['sequence'], uniprot[annot['norm']]['sequence']):
                    return True
            else:
                return False
        elif an['annottype'] == "e_1":
            # for species, check that they are nearby in the tree
            if int(an['norm']) in climb_tax_tree(int(annot['norm']), taxtree) or int(annot['norm']) in climb_tax_tree(int(an['norm']), taxtree):
                return True
            else:
                return False
    return False


def blast(seq1, seq2):
    # run blast on these sequences and if they have above 90% identity on 90% of the sequence, then consider them to be the same
    tmpdir = "/tmp"

    # this is not threadsafe, but we're only using one thread anyway
    tmp1 = tmpdir + "/tmp1.fna"
    tmp2 = tmpdir + "/tmp2.fna"
    with open(tmp1, 'w') as f1:
        f1.write(seq1 + "\n")
    with open(tmp2, 'w') as f2:
        f2.write(seq2 + "\n")

    command = "blastp -query " + tmp1 + " -subject " + tmp2 + " -outfmt '6 pident'"  # XXX check the match length
    process = Popen(shlex.split(command), stdout=PIPE)
    out, err = process.communicate()
    exit_code = process.wait()

    if out:
        o = out.split("\n")
        out = o[0]
        #print "pident " + str(out)
        if float(out) > 90:
            return True
        return False
    return False

def write_results(phase, user1, user2, document, annottype, annot):
    outfile = "results." + phase + "." + annottype
    with open(outfile, 'a') as f:
        f.write(document + "\t" + user1 + "\t" + user2 + "\t" + annot['text'] + "\t" + annot['norm'] + "\t" + annot['start'] + "\t" + annot['end'] + "\n")

def inter_annotator(annotations, uniprot, taxtree):
    # Print interannotator agreement for all documens 
    stats = {}
    for document in annotations.keys():
        for user1 in annotations[document].keys():
            if not user1 in stats:
                stats[user1] = {}
            for user2 in annotations[document].keys():
                if user1 == user2:
                    continue
                if not user2 in stats[user1]:
                    stats[user1][user2] = {}
                if not document in stats[user1][user2]:
                    stats[user1][user2][document] = {}
                    stats[user1][user2][document]['proteins'] = {}
                    stats[user1][user2][document]['species'] = {}

                    stats[user1][user2][document]['proteins']['tp'] = 0
                    stats[user1][user2][document]['proteins']['fp'] = 0
                    stats[user1][user2][document]['proteins']['fn'] = 0

                    stats[user1][user2][document]['proteins']['n_tp'] = 0
                    stats[user1][user2][document]['proteins']['n_fp'] = 0
                    stats[user1][user2][document]['proteins']['n_fn'] = 0

                    stats[user1][user2][document]['species']['tp'] = 0
                    stats[user1][user2][document]['species']['fp'] = 0
                    stats[user1][user2][document]['species']['fn'] = 0

                    stats[user1][user2][document]['species']['n_tp'] = 0
                    stats[user1][user2][document]['species']['n_fp'] = 0
                    stats[user1][user2][document]['species']['n_fn'] = 0

                for annot in annotations[document][user1]:
                    if annot['annottype'] == "e_2":
                        annottype = 'proteins'
                    else:
                        annottype = 'species'

                    match = same_boundaries(annot, annotations[document][user2])
                    if match:
                        stats[user1][user2][document][annottype]['n_tp'] += 1

                        if same_normalization(annot, match, uniprot, taxtree):
                            #print "agree " + user1 + " " + user2 + " " + annot['text'] + " " + annot['start'] + " " + annot['end']
                            stats[user1][user2][document][annottype]['tp'] += 1

                            write_results("tp", user1, user2, document, annottype, annot)

                        else:
                            #print "norm diff " + user1 + " " + user2 + " " + annot['text'] + " " + annot['start'] + " " + annot['end'] + " " + annot['norm'] + " " + match['norm']
                            stats[user1][user2][document][annottype]['fp'] += 1
                            stats[user1][user2][document][annottype]['fn'] += 1
                            write_results("wrongnorm", user1, user2, document, annottype, annot)

                    else:
                        #print "no match for " + user1 + " " + annot['text'] + " " + annot['start'] + " " + annot['end'] + " " + annot['norm']
                        stats[user1][user2][document][annottype]['fp'] += 1
                        stats[user1][user2][document][annottype]['n_fp'] += 1
                        write_results("wrongboundaries", user1, user2, document, annottype, annot)


                for annottype in ['proteins', 'species']:
                    stats[user1][user2][document][annottype]['fn'] += get_count(annotations[document][user2], annottype) - stats[user1][user2][document][annottype]['tp'] - stats[user1][user2][document][annottype]['fn']

                    stats[user1][user2][document][annottype]['n_fn'] += get_count(annotations[document][user2], annottype) - stats[user1][user2][document][annottype]['n_tp'] - stats[user1][user2][document][annottype]['n_fn']

                    stats[user1][user2][document][annottype]['precision'] = calc_precision(stats[user1][user2][document][annottype]['tp'], stats[user1][user2][document][annottype]['fp'])
                    stats[user1][user2][document][annottype]['recall'] = calc_recall(stats[user1][user2][document][annottype]['tp'], stats[user1][user2][document][annottype]['fn'] )

                    stats[user1][user2][document][annottype]['n_precision'] = calc_precision(stats[user1][user2][document][annottype]['n_tp'], stats[user1][user2][document][annottype]['n_fp'])
                    stats[user1][user2][document][annottype]['n_recall'] = calc_recall(stats[user1][user2][document][annottype]['n_tp'], stats[user1][user2][document][annottype]['n_fn'])

    #pprint.pprint(stats)
    return stats;

def get_count(annot, annottype):
    if annottype == "proteins":
        annottype = "e_2"
    if annottype == "species":
        annottype = "e_1"
    count = 0
    for a in annot:
        if a['annottype'] == annottype:
            count += 1
    return count

def calc_precision(tp, fp):
    if tp + fp == 0:
        return 0;
    else:
        return float(tp) / (tp + fp)

def calc_recall(tp, fn):
    if tp + fn == 0:
        return 0
    else:
        return float(tp) / (tp + fn)

def calc_fscore(r, p):
    if r+p == 0:
        return 0
    else:
        return 2*r*p / (r+p)

def print_stats(stats):
    #print "annottype\tuser1\tuser2\tdocument\tprecision\trecall\tn_precision\tn_recall"
    print "annottype\tuser1\tuser2\tdocument\tf score\t unnorm f score"
    for user1 in stats:
        for user2 in stats[user1]:
            
            for annottype in ['proteins', 'species']:
                ac_tp = 0
                ac_fp = 0
                ac_fn = 0

                ac_n_tp = 0
                ac_n_fp = 0
                ac_n_fn = 0

                for document in stats[user1][user2]:
                    recall = stats[user1][user2][document][annottype]['recall']
                    precision = stats[user1][user2][document][annottype]['precision']
                    n_recall = stats[user1][user2][document][annottype]['n_recall']
                    n_precision = stats[user1][user2][document][annottype]['n_precision']
                    tp = stats[user1][user2][document][annottype]['tp']
                    fp = stats[user1][user2][document][annottype]['fp']
                    fn = stats[user1][user2][document][annottype]['fn']
                    ac_tp += tp
                    ac_fp += fp
                    ac_fn += fn

                    ac_n_tp += stats[user1][user2][document][annottype]['n_tp']
                    ac_n_fp += stats[user1][user2][document][annottype]['n_fp']
                    ac_n_fn += stats[user1][user2][document][annottype]['n_fn']

                    fscore = 0
                    n_fscore = 0
                    if precision + recall > 0:
                        fscore = 2 * recall * precision / (precision + recall)
                        n_fscore = 2 * n_recall * n_precision / (n_precision + n_recall)
                    #print user1 + "\t" + user2 + "\t" + document + "\t" + str(tp) + "\t" + str(fp) + "\t" + str(fn) 
                    #print annottype + "\t" + user1 + "\t" + user2 + "\t" + document + "\t" + str(fscore) + "\t" + str(n_fscore)

                ac_precision = calc_precision(ac_tp, ac_fp)
                ac_recall = calc_recall(ac_tp, ac_fn)
                ac_fscore = calc_fscore(ac_recall, ac_precision)

                ac_n_precision = calc_precision(ac_n_tp, ac_n_fp)
                ac_n_recall = calc_recall(ac_n_tp, ac_n_fn)
                ac_n_fscore = calc_fscore(ac_n_recall, ac_n_precision)

                #print annottype + " overall " + user1 + "\t" + user2 + "\t" + str(ac_precision) + "\t" + str(ac_recall) + "\t" + str(ac_n_precision) + "\t" + str(ac_n_recall)
                print annottype + "\toverall\t" + user1 + "\t" + user2 + "\t" + str(ac_fscore) + "\t" + str(ac_n_fscore)
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

def read_project_list(filename):
    projects = []
    with open(filename, 'r') as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            print cols[0]
            projects.append(cols[0])
    return projects

# todo get these from ncbi_taxonomy module.... ???
def parse_taxtree(filename, returnlevel=False):
    with open(filename) as f:
        # read in the ncbi tax nodes.dmp file
        taxtree = {}
        taxlevel = {}
        for line in f:
            array = line.rstrip("\n").split("|")
            taxid = int(array[0].strip("\t"))
            parent = int(array[1].strip("\t"))
            level = array[2].strip("\t")
            taxtree[taxid] = parent
            taxlevel[taxid] = level
        if returnlevel:
            return (taxtree, taxlevel)
        else:
            return taxtree

def climb_tax_tree(taxid, taxtree):
    # climb the tree to the root and return the lineage
    taxid = int(taxid)

    if taxid in taxtree:
        parent = taxtree[taxid]
    else:
        # TODO this is a problem
        #print "id " + str(taxid) + " not found in taxtree.  this is bad"
        return []

    lineage = [taxid]
    while parent != 1:
        lineage.insert(0, parent)
        if parent in taxtree:
            parent = taxtree[parent]
        else:
            print "parent " + str(parent) + " not found in taxtree.  this is very bad"
            return -1
    return lineage

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

    parser.add_argument('-t', '--taxonomy',
                required=True,
                dest='taxonomy',
                help="ncbi taxonomy nodes.dmp")

    parser.add_argument('-o', '--output',
                required=False,
                dest='output',
                help="tagger output")
    
    parser.add_argument('-e', '--entities',
                required=False,
                dest='entities',
                help="entities file for textmining")

    args=parser.parse_args()

    annotations = defaultdict(lambda: defaultdict(list))

    if args.inputfile:
        for inputfile in args.inputfile:
            # add all annotations into one dictionary
            annotations = parse_annotations(inputfile, annotations)

    if args.output and args.entities:
        annotations = parse_tagger(args.output, args.entities, annotations)

    uniprot = parse_uniprot_single_xml(args.uniprot)
    taxtree = parse_taxtree(args.taxonomy)
    stats = inter_annotator(annotations, uniprot, taxtree)
    print_stats(stats)


if __name__ == "__main__":
    main()
