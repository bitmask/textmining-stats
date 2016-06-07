import argparse
from collections import defaultdict


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
                    annot = {'text': text, 'norm': normid, 'start': start, 'end': end, 'user': user, 'docid': docid}
                    annotations[docid][user].append(annot)
    return annotations

def inter_annotator(annotations):
    for document in annotations.keys():
        if len(annotations[document].keys()) > 1:
            print document
            print annotations[document].keys()
    return True;

def print_stats(stats):
    return True

def main():
    parser = argparse.ArgumentParser(description="Calculate interannotator etc stats")
    # TODO arguments should be the same as the tagger
    parser.add_argument('-f', '--files',
                required=True,
                default = [],
                nargs='+',
                dest='inputfile',
                help="list of annotation files")

    args=parser.parse_args()

    annotations = defaultdict(lambda: defaultdict(list))
    if args.inputfile:
        for inputfile in args.inputfile:
            # add all annotations into one dictionary
            annotations = parse_annotations(inputfile, annotations)

    stats = inter_annotator(annotations)
    print_stats(stats)


if __name__ == "__main__":
    main()
