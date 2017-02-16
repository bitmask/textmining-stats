import argparse
from collections import defaultdict
import re
import Bio
import shlex
from subprocess import Popen, PIPE
import ncbi_taxonomy.ncbi_taxonomy as tax
import pprint
import numpy as np
import glob
import json
import os


# RUN WITH Python 2.0 or fix '' == u'' comparisons of strings to unicode strings

################################################################################
# This code is written by Helen Cook, with contributions from ...
# 
# Distributed under the BSD license. 
################################################################################

def parse_corpus(inputfile):
    corpus = {}
    with open(inputfile, "r") as f:
        for line in f:
            pmid, author, journal, year, text = line.rstrip("\n").split("\t")
            pmid = pmid.lstrip("PMID:")
            corpus[pmid] = text
    return corpus

def get_user_from_inputfile(inputfile):
    field = os.path.basename(inputfile).split(".")
    if field[1] == "test2":
        return "Helen"
    user = field[1].split("-")
    if user[0] == "Kerstin":
        return "Cristina"
    return user[0]

def parse_tagtog_document(pmid, path):
    # get the number of characters in the first paragraph -- this is the offset to add to the coords in the second paragraph to make global coords
    document = ""
    for docfile in glob.glob(path + "/doc.*-" + pmid + '.txt'):
        # should only be one file
        with open(docfile, "r") as f:
            for line in f:
                document += line
                # is already missing the newline
    # split on abstract separator
    p = document.split("Abstract ")
    if len(p) == 1:
        return (p[0], "") # document may be title only
   
    # omg
    paragraphs = []
    paragraphs.append(p[0] + "Abstract ") # add the title
    remaining_text = p[1]

    if 'BACKGROUND' in remaining_text:
        q = remaining_text.split("BACKGROUND ")
        remaining_text = q[1]
        paragraphs[-1] = paragraphs[-1] + "BACKGROUND " # update the title block to include this header
    if 'OBJECTIVES' in remaining_text:
        q = remaining_text.split("OBJECTIVES ")
        paragraphs.append(q[0] + "OBJECTIVES ")
        remaining_text = q[1]
    if 'STUDY DESIGN' in remaining_text:
        q = remaining_text.split("STUDY DESIGN ")
        paragraphs.append(q[0] + "STUDY DESIGN ")
        remaining_text = q[1]
    if 'RESULTS' in remaining_text:
        q = remaining_text.split("RESULTS ")
        paragraphs.append(q[0] + "RESULTS ")
        remaining_text = q[1]
    if 'CONCLUSION' in remaining_text:
        if 'CONCLUSIONS' in remaining_text:
            q = remaining_text.split("CONCLUSIONS ")
            paragraphs.append(q[0] + "CONCLUSIONS ")
            remaining_text = q[1]
        else:
            q = remaining_text.split("CONCLUSION ")
            paragraphs.append(q[0] + "CONCLUSION ")
            remaining_text = q[1]
    if 'UNLABELLED' in remaining_text:
        q = remaining_text.split("UNLABELLED ")
        remaining_text = q[1]
        paragraphs[-1] = paragraphs[-1] + "UNLABELLED "
    if 'IMPORTANCE' in remaining_text:
        q = remaining_text.split("IMPORTANCE ")
        paragraphs.append(q[0] + "IMPORTANCE ")
        remaining_text = q[1]

    paragraphs.append(remaining_text)
    return paragraphs

def get_position_offset(pmid, para_name, path):
    if para_name == "s1h1":
        return 0

    paragraphs = parse_tagtog_document(pmid, path)

    if len(paragraphs) == 2:
        return len(paragraphs[0])

    else:
        r = re.compile('s2s(.)p1')
        which = r.findall(para_name)[0] # this is the index of the paragraph the match is in

        idx = 0
        offset = 0
        while(idx < int(which)):
            offset += len(paragraphs[idx]) # add up the lengths of all preceeding paragraphs
            idx += 1
        return offset

def update_deprecated_taxids(taxid):
    # polyoma viruses were reclassified after tagging started, and old entries aren't in nodes.dmp
    if taxid == '10633':
        return '1891767'
    if taxid == '10624':
        return '36362'
    if taxid == '10634':
        return '1891730'
    if taxid == '10629':
        return '1891762'
    return taxid

def parse_annotations_json(inputfile, annotations):
    user=get_user_from_inputfile(inputfile)
    path = os.path.dirname(inputfile)
    with open(inputfile) as f:
        data = json.load(f)

        pmid = data['sources'][0]['id']
        
        for an in data['entities']:
            annottype = an['classId']
            text = an['offsets'][0]['text'] # TODO index?

            normid = "???"
            for k,n in an['normalizations'].iteritems():
                normid = n['source']['id']

            norms = re.compile(", ?").split(normid)
            for normid in norms:

                normid = update_deprecated_taxids(normid)

                para_name = an['part']
                offset = get_position_offset(pmid, para_name, path)
                start = an['offsets'][0]['start']
                end = int(start) + len(text) # TODO off by n?

                
                annot = {'text': text, 'norm': normid, 'start': str(int(start) + offset), 'end': str(int(end) + offset), 'user': user, 'pmid': pmid, 'annottype': annottype}
                annotations[pmid][user].append(annot)
    return annotations

        


def parse_annotations(inputfile, annotations):
    # Input format is one annotation per line
    user = get_user_from_inputfile(inputfile)
    with open(inputfile, 'r') as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 5:
                (annotid, doctype, pmid, annottype, annotdict, normid, blank, userX, unused, text, para_name, start, end) = cols
                if annotdict == "NCBITax" or annotdict == "Uniprot":
                    path = os.path.dirname(inputfile)
                    offset = get_position_offset(pmid, para_name, path)
                    #print "offset " + str(offset) + " for " + para_name + " " + start + " " + end + " " + text
                    
                    # fix deprecated ncbi tax entries
                    normid = update_deprecated_taxids(normid)

                    annot = {'text': text, 'norm': normid, 'start': str(int(start) + offset), 'end': str(int(end) + offset), 'user': user, 'pmid': pmid, 'annottype': annottype}
                    annotations[pmid][user].append(annot)
    return annotations

def convert_annottype(annottype):
    if int(annottype) == -2 or int(annottype) == -3:
        return "e_1" # species
    else:
        return "e_2" # protein
    if annottype == 0:
        return False
    
def parse_tagger(tagger_output_file, entities_file, annotations):
    #Parse the tagger output file and also include the x_entities entries
    entities = {}
    with open(entities_file, "r") as f:
        for line in f:
            (serialno, annottype, norm) = line.rstrip("\n").split("\t")
            entities[serialno] = [annottype, norm]

    user = 'tagger' #For identification from other users.
    with open(tagger_output_file, "r") as f:
        for line in f:
            (pmid, unused, unused2, start, end, annot, taxid, serialno) = line.rstrip("\n").split("\t")
            if serialno in entities:
                (annottype, normid) = entities[serialno]
                annottype = convert_annottype(annottype)
                if annottype:
                    # fix deprecated ncbi tax entries
                    normid = update_deprecated_taxids(normid)

                    annot = {'text': annot, 'norm': normid, 'start': str(int(start) + 1), 'end': str(int(end)+0), 'user': user, 'pmid': pmid, 'annottype': annottype}
                    annotations[pmid][user].append(annot)
            else:
                # it's a group
                pass
                    
    return annotations

def dedup_tagger(annotations, uniprot, taxtree, taxlevel, unreviewed):
    # take the tagger output and for all annotations that have the same boundaries, dedup normalizations
    for pmid in annotations.keys():
        for user in annotations[pmid].keys():
            if user == 'tagger':
                # get all annotations at each position
                annot_at_pos = {}
                for new in annotations[pmid][user]:
                    same_bdy = False
                    for key,existing in annot_at_pos.iteritems():
                        if same_boundaries(new, existing):
                            annot_at_pos[key].append(new)
                            same_bdy = True
                    if not same_bdy:
                        key = new['start'] + "-" + new['end']
                        annot_at_pos[key] = [new]
                new_annot = []
                for pos,thing in annot_at_pos.iteritems():
                    new_annot += collapse_same_norm(thing, uniprot, taxtree, taxlevel, unreviewed)

                annotations[pmid][user] = new_annot

    return annotations

def convert_tagger_bytes_to_char(annotations, corpus):
    for document in annotations:
        for user in annotations[document]:
            if user == "tagger":
                mapping = get_b2c_mapping(corpus[document])
                converted_annotations = []
                for annot in annotations[document][user]:
                    annot['start'] = str(mapping[int(annot['start'])] + 1)
                    annot['end'] = str(mapping[int(annot['end'])] + 1)
                    converted_annotations.append(annot)
                annotations[document][user] = converted_annotations
    return annotations

def get_b2c_mapping(document):
    # get byte to character mapping for document
    u_results = []
    # for each byte in the document, determine if it belongs to a multi byte character
    mapping = {} # byte to char
    byte = 0
    character = 0
    u_document = document.decode("utf-8") # turn bytes into characters
    for b in u_document:
        u = b.encode('utf-8') # back to bytes
        char_bytes = len(u) # how many bytes does this character consist of
        for i in range(0, char_bytes):
            mapping[byte+i] = character
        byte += char_bytes
        character += 1  
    return mapping

def same_boundaries(annot, annot_list):
    # check whether any annotations with the same boundaries as annot exists on the annot_list
    returnlist = []
    for an in annot_list:
        if an['end'] == annot['end'] and an['start'] == annot['start']:
            returnlist.append(an)
        # allow fuzzy boundaries, but insist that one matches
        elif an['end'] == annot['end']:
            returnlist.append(an)
        elif an['start'] == annot['start']:
            returnlist.append(an)
        # allow any overlap
        elif int(an['start']) < int(annot['start']) < int(an['end']) < int(annot['end']):
            returnlist.append(an)
        elif int(annot['start']) < int(an['start']) < int(annot['end']) < int(an['end']):
            returnlist.append(an)
        # subset
        elif int(annot['start']) < int(an['start']) < int(an['end']) < int(annot['end']):
            returnlist.append(an)
        # superset
        elif int(an['start']) < int(annot['start']) < int(annot['end']) < int(an['end']):
            returnlist.append(an)
    return returnlist

def get_unreviewed_sequence(norm, unreviewed):
    if norm in unreviewed:
        return unreviewed[norm]['sequence']
    else:
        return ''

def where_trees_meet(tree1, tree2):
    lowest_same = 0
    for t,u in zip(tree1, tree2):
        if t == u:
            lowest_same = t
    return int(lowest_same)

def same_normalization(annot, an, uniprot, taxtree, taxlevel, unreviewed):
    if an['norm'] == annot['norm']:
        return True
    else:
        if annot['norm'] == '???' or an['norm'] == '???' or annot['norm'] == 'UNKN' or an['norm'] == 'UNKN' or annot['norm'] == "INVALID_ID" or an['norm'] == "INVALID_ID":
            # if a normalization was not entered by the user, then say they are not the same
            return False
        if an['annottype'] == "e_2" and annot['annottype'] == "e_2":
            # for proteins, check that they are close enough to each other
            sequence1 = ''
            sequence2 = ''

            if an['norm'] in uniprot:
                sequence1 = uniprot[an['norm']]['sequence']
            else:
                sequence1 = get_unreviewed_sequence(an['norm'], unreviewed)

            if annot['norm'] in uniprot:
                sequence2 = uniprot[annot['norm']]['sequence']
            else:
                sequence2 = get_unreviewed_sequence(annot['norm'], unreviewed)


            if sequence1 and sequence2:
                #print an['norm'] + " " + annot['norm']
                if blast(sequence1, sequence2, 90):
                    return True
            else:
                #print "not in uniprot " + an['norm'] + " or " + annot['norm'] + " " + an['pmid']
                return False
        elif an['annottype'] == "e_1" and annot['annottype'] == "e_1":
            # for species, check that they are nearby in the tree
            try:
                tree1 = tax.climb_tax_tree(int(annot['norm']), taxtree)
                tree2 = tax.climb_tax_tree(int(an['norm']), taxtree)

                if 10509 in tree1 and 10509 in tree2:
                    # exception for Adenovirus
                    return True

                level = where_trees_meet(tree1, tree2)
                below_species = is_below_species(level, taxtree, taxlevel)
                if taxlevel[level].lower() == 'species' or below_species == True:
                    return True

                return False
            except:
                return False # if someone has put a protein id where a taxid should go ??
    return False

def is_below_species(taxid, taxtree, taxlevel):
    below_species = False
    if taxlevel[taxid].lower() == 'no rank':
        tree = tax.climb_tax_tree(taxid, taxtree)
        tree.reverse() # go up the tree
        for t in tree:
            if taxlevel[t] == 'species':
                below_species = True
    return below_species



def blast(seq1, seq2, pident):
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
        if float(out) > pident:
            return True
        return False
    return False

def clean_up_old_results():
    for phase in ["tp", "wrongnorm", "wrongboundaries"]:
        for who in ["humans", "tagger", "consensus"]:
            for annottype in ["proteins", "species"]:
                out_file = "results." + who + "." + phase + "." + annottype
                if os.path.exists(out_file):
                    os.remove(out_file)

def write_results(phase, user1, user2, document, annottype, annot):

    if user1 == "tagger":
        who = "tagger"
    elif user2 == "tagger":
        who = "consensus"
    else:
        who = "humans"

    outfile = "results." + who + "." + phase + "." + annottype

    with open(outfile, 'a') as f:
        try:
            f.write(document.encode('utf-8') + "\t")
            f.write(user1 + "\t")
            f.write(user2 + "\t")
            f.write(annot['text'].encode('utf-8') + "\t")
            f.write(annot['norm'].encode('utf-8') + "\t")
            f.write(annot['start'] + "\t")
            f.write(annot['end'] + "\n")
        except:
            print "could not write"
            pprint.pprint(annot)

def report_to_user(a):
    report_file = "results.report"
    with open(report_file, "a") as f:
        f.write(a['user'] + "\t")
        f.write(str(a['pmid']) + "\t")
        f.write(a['text'].encode('utf-8') + "\t")
        f.write(a['start'] + "\t")
        f.write(a['end'] + "\t")
        f.write(a['norm'] + "\t")
        f.write(a['annottype'] + "\n")

def filter_annotations(annotations, taxtree, taxlevel):
    # remove annotations to things at a higher tax level than species
    # and check that annotations are well formed
    for document in annotations.keys():
        for user1 in annotations[document].keys():
            clean_annotations = []
            for a in annotations[document][user1]:

                if a['annottype'] == "e_1":
                    if a['norm'] == '???' or a['norm'] == 'UNKN' or a['norm'] == "INVALID_ID" or a['norm'] == '0':
                        # add without checking
                        clean_annotations.append(a)
                        report_to_user(a)
                    else:

                        # there should not actually be any multinorms at this point
                        norms = re.compile(", ?").split(a['norm'])
                        add = []
                        for n in norms:
                            try:
                                int(n)
                            except:
                                n = 0
                                report_to_user(a)
                            if n:
                                t = taxlevel[int(n)].lower()
                                below_species = is_below_species(int(n), taxtree, taxlevel)
                                if t == "species" or below_species:
                                    add.append(str(n))
                                else:
                                    #t == "genus" or t == "subfamily" or t == "family" or t == "order":
                                    pass

                        if add:
                            add.sort()
                            a['norm'] = ",".join(add)
                            clean_annotations.append(a)

                elif a['annottype'] == "e_2": 
                    # add all the proteins
                    clean_annotations.append(a)

                    # but report problems for some of them
                    is_int = True
                    try:
                        int(a['norm'])
                    except:
                        is_int = False
                    if is_int:
                        report_to_user(a)

                    norms = re.compile(", ?").split(a['norm'])
                    a['norm'] = ",".join(norms)

                    if a['norm'] == '???' or a['norm'] == 'UNKN' or a['norm'] == "INVALID_ID" or a['norm'] == '0':
                        report_to_user(a)

            annotations[document][user1] = clean_annotations
    return annotations

def get_positions(annot):
    positions = defaultdict(lambda: [])
    for user in annot.keys():
        for an in annot[user]:
            start = an['start']
            if not start in positions:
                positions[start] = []
            positions[start].append(an)
    return positions

def get_norms(annot):
    # get all the normalizations at this position
    norms = defaultdict(lambda: [])
    for an in annot:
        norms[ an['norm'] ].append(an)
    return norms

def different_annotators(an1, an2):
    if an1[0]['user'] == an2[0]['user']: #TODO
        return False
    return True

def get_most_support(new):
    norms = get_norms(new)
    sort = sorted(norms.values(), key=lambda x : len(x))
    if len(sort) == 1:
        return sort[-1][0]['norm']
    else:
        if len(sort[-2]) == len(sort[-1]):
            # if there is a tie for first place
            if different_annotators(sort[-2], sort[-1]):
                # and they don't come from the same person
                return 'tie'
            else:
                return sort[-1][0]['norm']
        else:
            # else take the norm from the first norm (they are all the same anyway) in the biggest list
            return sort[-1][0]['norm']

def get_max_end(annot):
    end = 0
    for an in annot:
        if an['end'] > end:
            end = an['end']
    return end

def get_annotators(array):
    annot = []
    for a in array:
        annot.append(a['user'])
    return list(set(annot))

def get_disagreeing_norms(array):
    annot = []
    for a in array:
        annot.append(a['norm'])
    return list(set(annot))

def collapse_same_norm(items, uniprot, taxtree, taxlevel, unreviewed):
    # assumes no user has normalized an entity to two normalizations that are equivalent (why would you do this anyway?)
    deduped = []
    for new in items:
        on_list = False
        for existing in deduped:
            if same_normalization(new, existing, uniprot, taxtree, taxlevel, unreviewed):
                on_list = True
        if not on_list:
            deduped.append(new)
    return deduped

def human_consensus(annotations, uniprot, taxtree, taxlevel, unreviewed):
    con = {}
    for document in annotations.keys():
        #print "document: " + str(document)
        consensus = []

        if "tagger" in annotations[document]:
            annotations[document].pop("tagger")

        # spread each annotation over all the positions in the document that it covers
        docarray = defaultdict(lambda: [])
        for user in annotations[document]:
            for annot in annotations[document][user]:
                for i in range(int(annot['start']), int(annot['end'])):
                    docarray[i].append(annot)

        normlist = []
        laststart = ''

        for idx in range(0, 3000): # no document is longer than 3000 bytes
            if idx in docarray:
                # in a block of annotations
                for li in docarray[idx]:
                    normlist.append(li)
                if not laststart:
                    laststart = docarray[idx][0]['start']
            else:
                if laststart:
                    # this is the first position after the end of a contiguous block
                    prev = idx - 1
                    representative = docarray[prev][0]
                    annot = {}
                    annot['annottype'] = representative['annottype']
                    annot['end'] = str(idx)
                    annot['start'] = laststart
                    annot['pmid'] = representative['pmid']
                    annot['text'] = representative['text'] # TODO fix end position, but this isn't used for anything
                    annot['user'] = "consensus"

                    # get unique hashes, because hashes are not hashable in python
                    tuples = tuple(frozenset(d.iteritems()) for d in normlist)
                    unique = set(tuples)
                    new = [dict(pairs) for pairs in unique]

                    # norms that are fuzzily the same are collapsed here
                    deduped = collapse_same_norm(new, uniprot, taxtree, taxlevel, unreviewed)

                    annot['norm'] = get_most_support(deduped)

                    annot['consensus'] = []
                    for n in new:
                        annot['consensus'].append({'user': n['user'], 'norm': n['norm']})


                    consensus.append(annot)
                    laststart = ''
                    normlist = []
   

        con[document] = consensus
    return con

def inter_annotator(annotations, uniprot, taxtree, taxlevel, unreviewed):
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

                    stats[user1][user2][document]['proteins']['u_n_tp'] = 0 # normalization only of proteins in uniprot
                    stats[user1][user2][document]['proteins']['u_n_fp'] = 0
                    stats[user1][user2][document]['proteins']['u_n_fn'] = 0

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

                    matches = same_boundaries(annot, annotations[document][user2])
                    if matches:
                        # check all matches to see if one has the same annotation
                        # report tp if one does
                        # report fp if no match has the same annotation
                        same_annot = False
                        stats[user1][user2][document][annottype]['n_tp'] += 1

                        for match in matches:
                            if same_normalization(annot, match, uniprot, taxtree, taxlevel, unreviewed):
                                same_annot = True

                        if same_annot:
                            #print "agree " + user1 + " " + user2 + " " + annot['text'] + " " + annot['start'] + " " + annot['end']
                            stats[user1][user2][document][annottype]['tp'] += 1
                            if annottype == 'proteins':
                                stats[user1][user2][document][annottype]['u_n_tp'] += 1
                            write_results("tp", user1, user2, document, annottype, annot)

                        else:
                            #print "norm diff " + user1 + " " + user2 + " " + annot['text'] + " " + annot['start'] + " " + annot['end'] + " " + annot['norm'] + " " + match['norm']
                            stats[user1][user2][document][annottype]['fp'] += 1
                            #stats[user1][user2][document][annottype]['fn'] += 1
                            write_results("wrongnorm", user1, user2, document, annottype, annot)

                            if annot['annottype'] == 'proteins' and annot['norm'] in uniprot:
                                stats[user1][user2][document]['proteins']['u_n_fp'] += 1


                    else:
                        stats[user1][user2][document][annottype]['fp'] += 1
                        stats[user1][user2][document][annottype]['n_fp'] += 1
                        if annottype == 'proteins':
                            stats[user1][user2][document][annottype]['u_n_fp'] += 1
                        write_results("wrongboundaries", user1, user2, document, annottype, annot)



                for annottype in ['proteins', 'species']:
                    stats[user1][user2][document][annottype]['fn'] += get_false_negs(get_count(annotations[document][user2], annottype), stats[user1][user2][document][annottype]['tp'])

                    if annottype == 'proteins':
                        stats[user1][user2][document][annottype]['u_n_fn'] += get_false_negs(get_count(annotations[document][user2], annottype), stats[user1][user2][document][annottype]['u_n_tp'])

                    stats[user1][user2][document][annottype]['n_fn'] += get_false_negs(get_count(annotations[document][user2], annottype), stats[user1][user2][document][annottype]['n_tp'])

                    stats[user1][user2][document][annottype]['precision'] = calc_precision(stats[user1][user2][document][annottype]['tp'], stats[user1][user2][document][annottype]['fp'])
                    stats[user1][user2][document][annottype]['recall'] = calc_recall(stats[user1][user2][document][annottype]['tp'], stats[user1][user2][document][annottype]['fn'] )

                    stats[user1][user2][document][annottype]['n_precision'] = calc_precision(stats[user1][user2][document][annottype]['n_tp'], stats[user1][user2][document][annottype]['n_fp'])
                    stats[user1][user2][document][annottype]['n_recall'] = calc_recall(stats[user1][user2][document][annottype]['n_tp'], stats[user1][user2][document][annottype]['n_fn'])

                stats[user1][user2][document]['proteins']['u_n_precision'] = calc_precision(stats[user1][user2][document]['proteins']['u_n_tp'], stats[user1][user2][document]['proteins']['u_n_fp'])
                stats[user1][user2][document]['proteins']['u_n_recall'] = calc_recall(stats[user1][user2][document]['proteins']['u_n_tp'], stats[user1][user2][document]['proteins']['u_n_fn'] )

    return stats;

def get_false_negs(count, tp):
    # get the number of false negatives from the number of annotations the other user made, and the number of true positives
    if count > tp:
        return count - tp
    else:
        return 0

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

def print_consensus(con):
    outputfile = "results.consensus"
    with open(outputfile, "w") as f:
        for pmid in con:
            for an in con[pmid]:
                if an['norm'] == 'tie':
                    t = an['text'].encode('utf-8')
                    f.write(pmid + "\t" + str(t) + "\t" + an['start'] + "\t" + an['end'] + "\t" + str(an['consensus']) + "\n")


def print_stats(stats, report_p_r):
    #print "annottype\tuser1\tuser2\tdocument\tprecision\trecall\tn_precision\tn_recall"
    print "annottype\tuser1\tuser2\tdocument\tf score\t unnorm f score"
    store = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))
    for user1 in stats:
        for user2 in stats[user1]:
            
            for annottype in ['proteins', 'species']:
                ac_tp = 0
                ac_fp = 0
                ac_fn = 0

                ac_n_tp = 0
                ac_n_fp = 0
                ac_n_fn = 0

                ac_u_n_tp = 0
                ac_u_n_fp = 0
                ac_u_n_fn = 0

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

                    if annottype == 'proteins':
                        ac_u_n_tp += stats[user1][user2][document][annottype]['u_n_tp']
                        ac_u_n_fp += stats[user1][user2][document][annottype]['u_n_fp']
                        ac_u_n_fn += stats[user1][user2][document][annottype]['u_n_fn']


                # across all documents
                ac_precision = calc_precision(ac_tp, ac_fp)
                ac_recall = calc_recall(ac_tp, ac_fn)
                ac_fscore = calc_fscore(ac_recall, ac_precision)

                ac_n_precision = calc_precision(ac_n_tp, ac_n_fp)
                ac_n_recall = calc_recall(ac_n_tp, ac_n_fn)
                ac_n_fscore = calc_fscore(ac_n_recall, ac_n_precision)

                store[user1][user2][annottype]['tp'] = ac_tp
                store[user1][user2][annottype]['n_tp'] = ac_n_tp
                store[user1][user2][annottype]['fp'] = ac_fp
                store[user1][user2][annottype]['n_fp'] = ac_n_fp
                store[user1][user2][annottype]['fn'] = ac_fn
                store[user1][user2][annottype]['n_fn'] = ac_n_fn

                if annottype == 'proteins':
                    store[user1][user2][annottype]['u_n_tp'] = ac_u_n_tp
                    store[user1][user2][annottype]['u_n_fp'] = ac_u_n_fp
                    store[user1][user2][annottype]['u_n_fn'] = ac_u_n_fn


    for user1 in store:
        for user2 in store[user1]:
            for annottype in ['proteins', 'species']:
                recall = 0
                prec = 0
                n_recall = 0
                n_prec = 0
                u_n_recall = 0
                u_n_prec = 0
                if annottype not in store[user1][user2]:
                    if annottype not in store[user2][user1]:
                        ac_fscore = 1
                        ac_n_fscore = 1
                    else:
                        ac_fscore = "NA"
                        ac_n_fscore = "NA"
                else:

                    # TODO check that there is no document that user1 has tagged that the tagger has not tagged not manually
                    if user1 == 'tagger' and user2 == 'consensus':
                        store['tagger']['consensus']['species']['fn'] += 9

                    if annottype not in store[user2][user1]:
                        ac_fscore = "NA"
                        ac_n_fscore = "NA"
                    else:
                        if report_p_r:
                            prec = calc_precision( store[user1][user2][annottype]['tp'], store[user1][user2][annottype]['fp']) 
                            n_prec = calc_precision( store[user1][user2][annottype]['n_tp'], store[user1][user2][annottype]['n_fp'])

                            recall = calc_recall( store[user1][user2][annottype]['tp'], store[user1][user2][annottype]['fn']) 
                            n_recall = calc_recall( store[user1][user2][annottype]['n_tp'], store[user1][user2][annottype]['n_fn'])
                            tp = store[user1][user2][annottype]['n_tp']
                            fn = store[user1][user2][annottype]['n_fn']
                            print annottype + " " + user1 + " " + user2 + " n_tp " + str(tp) + " n_fn " + str(fn)
                            print "n_recall " + str(n_recall)

                            ac_fscore = calc_fscore(prec, recall)
                            ac_n_fscore = calc_fscore(n_prec, n_recall)

                        else:
                            prec = calc_precision( store[user1][user2][annottype]['tp'] + store[user2][user1][annottype]['tp'] , store[user1][user2][annottype]['fp'] + store[user2][user1][annottype]['fp'])
                            n_prec = calc_precision( store[user1][user2][annottype]['n_tp'] + store[user2][user1][annottype]['n_tp'] , store[user1][user2][annottype]['n_fp'] + store[user2][user1][annottype]['n_fp'])
                            ac_fscore = calc_fscore(prec, prec)
                            ac_n_fscore = calc_fscore(n_prec, n_prec)
                    #print annottype + " overall " + user1 + "\t" + user2 + "\t" + str(ac_precision) + "\t" + str(ac_recall) + "\t" + str(ac_n_precision) + "\t" + str(ac_n_recall)
                if annottype == "proteins":
                    if report_p_r:
                        u_n_prec = calc_precision(store[user1][user2][annottype]['u_n_tp'], store[user1][user2][annottype]['u_n_fp'])
                        u_n_recall = calc_precision(store[user1][user2][annottype]['u_n_tp'], store[user1][user2][annottype]['u_n_fn'])
                        ac_u_n_fscore = calc_fscore(u_n_prec, u_n_recall)

                    else:
                        u_n_prec = calc_precision( store[user1][user2][annottype]['u_n_tp'] + store[user2][user1][annottype]['u_n_tp'], store[user1][user2][annottype]['u_n_fp'] + store[user2][user1][annottype]['u_n_fp'])
                        ac_u_n_fscore = calc_fscore(u_n_prec, u_n_prec)
                if annottype == "species":
                    ac_u_n_fscore = 0
                    u_n_prec = 0
                    u_n_recall = 0

                if report_p_r:
                    if user1 == "tagger":
                        print annottype + "\t" + user1 + "\t" + user2 + "\twith_normalization\t" + str(prec) + "\t" + str(recall) 
                        print annottype + "\t" + user1 + "\t" + user2 + "\tboundaries_only\t" + str(n_prec) + "\t" + str(n_recall) 
                        print annottype + "\t" + user1 + "\t" + user2 + "\tnorm_in_uniprot\t" + str(u_n_prec) + "\t" + str(u_n_recall)
                    #print annottype + "\toverall\t" + user1 + "\t" + user2 + "\t" + str(ac_fscore) + "\t" + str(ac_n_fscore) + "\t" + str(ac_u_n_fscore) + "\t" + str(prec) + "\t" + str(recall) + "\t" + str(n_prec) + "\t" + str(n_recall) + "\t" + str(u_n_prec) + "\t" + str(u_n_recall)
                else:
                    print annottype + "\toverall\t" + user1 + "\t" + user2 + "\t" + str(ac_fscore) + "\t" + str(ac_n_fscore) + "\t" + str(ac_u_n_fscore)
    return True


def parse_uniprot_proteome_dir(dirname):
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

def combine_annot(consensus, tagger):
    for doc, item in consensus.iteritems():
        tagger[doc]['consensus'] = item
    return tagger

def print_unreviewed_proteins(annotations, uniprot):
    with open("results.unreviewed_proteins", "w") as f:
        for document in annotations.keys():
            for user in annotations[document].keys():
                for annot in annotations[document][user]:
                    if annot['annottype'] == 'e_2':
                        if not annot['norm'] in uniprot:
                            f.write(annot['norm'].encode('utf-8') + "\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate interannotator etc stats")
    # TODO arguments should be the same as the tagger

    parser.add_argument('-f', '--files',
                required=True,
                default = [],
                nargs='+',
                dest='inputfile',
                help="list of annotation files")

    parser.add_argument('-c', '--corpus',
                required=False,
                default = "",
                dest='corpus',
                help="corpus file that tagger reads")

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
    tagger_annot = defaultdict(lambda: defaultdict(list))

    if args.inputfile:
        for inputfile in args.inputfile:
            # add all annotations into one dictionary
            annotations = parse_annotations_json(inputfile, annotations)

    clean_up_old_results()

    uniprot = parse_uniprot_single_xml(args.uniprot)

    unrev_path = "../../../uniprot_to_payload/data_in/proteomes_unreviewed/"
    unreviewed = parse_uniprot_proteome_dir(unrev_path)

    taxtree, taxlevel = tax.parse_taxtree(args.taxonomy, True)

    if args.output and args.entities and args.corpus:
        tagger_annot = parse_tagger(args.output, args.entities, tagger_annot)
        tagger_annot = dedup_tagger(tagger_annot, uniprot, taxtree, taxlevel, unreviewed)
        corpus = parse_corpus(args.corpus)
        tagger_annot = convert_tagger_bytes_to_char(tagger_annot, corpus)
        tagger_annot = filter_annotations(tagger_annot, taxtree, taxlevel)

    annotations = filter_annotations(annotations, taxtree, taxlevel) # remove annotations for species that are above species level

    # get the list of unreviewed proteins to download and parse above
    #print_unreviewed_proteins(annotations, uniprot)

    consensus = human_consensus(annotations, uniprot, taxtree, taxlevel, unreviewed)
    #pprint.pprint(annotations['3040055']['Helen'])
    #pprint.pprint(consensus['3040055'])
    #pprint.pprint(annotations['9191870']['Rudolfs'])
    #pprint.pprint(consensus['9191870'])
    #print " Rudolfs ***************************************"
    #pprint.pprint(annotations['16227217']['Rudolfs'])
    #print " Juanmi ***************************************"
    #pprint.pprint(annotations['16227217']['Juanmi'])
    #print " Consensus ***************************************"
    #pprint.pprint(consensus['16227217'])
    #print " cristina ***************************************"
    #pprint.pprint(annotations['15680420']['Cristina'])
    #print " helen ***************************************"
    #pprint.pprint(annotations['15680420']['Helen'])
    #print " Consensus ***************************************"
    #pprint.pprint(consensus['15680420'])
    iaa = inter_annotator(annotations, uniprot, taxtree, taxlevel, unreviewed)
    print_stats(iaa, False)
    print_consensus(consensus)
    print "combined results"
    combined = combine_annot(consensus, tagger_annot)
    #pprint.pprint(combined)
    stats = inter_annotator(combined, uniprot, taxtree, taxlevel, unreviewed)
    print_stats(stats, True)


if __name__ == "__main__":
    main()
