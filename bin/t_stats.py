import unittest
import ncbi_taxonomy.ncbi_taxonomy as tax

from collections import defaultdict
import pprint

import os, sys
lib_path = os.path.abspath('.')
sys.path.append(lib_path)

from stats import *

consensus_testcases = {
    'agree': { 'humans': {'user1': [{'start': '10', 'end': '20', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': 'humans'}],
                          'user2': [{'start': '10', 'end': '20', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': 'humans'}],
                         },
               'expected': [{'start': '10', 'end': '20', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': 'consensus', 'pmid': 'humans'}],
             },

    'di_n2': { 'humans': {'user1': [{'start': '10', 'end': '20', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': 'humans'}],
                          'user2': [{'start': '10', 'end': '20', 'norm': '11234', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': 'humans'}],
                         },
               'expected': [{'start': '10', 'end': '20', 'norm': '0', 'annottype': 'e_1', 'text': 'foo', 'user': 'consensus', 'pmid': 'humans'}],
             },
    'di_n3': { 'humans': {'user1': [{'start': '10', 'end': '20', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': 'humans'}],
                          'user2': [{'start': '10', 'end': '20', 'norm': '11234', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': 'humans'}],
                          'user3': [{'start': '10', 'end': '20', 'norm': '11234', 'annottype': 'e_1', 'text': 'foo', 'user': '3', 'pmid': 'humans'}],
                         },
               'expected': [{'start': '10', 'end': '20', 'norm': '11234', 'annottype': 'e_1', 'text': 'foo', 'user': 'consensus', 'pmid': 'humans'}],
             },

    'di_p1': { 'humans': {'user1': [{'start': '10', 'end': '20', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': 'humans'}],
                          'user2': [{'start': '10', 'end': '30', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': 'humans'}],
                         },
               'expected': [{'start': '10', 'end': '30', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': 'consensus', 'pmid': 'humans'}],
             },
    'di_p2': { 'humans': {'user1': [{'start': '10', 'end': '20', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': 'humans'}],
                          'user2': [{'start': '15', 'end': '30', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': 'humans'}],
                         },
               'expected': [{'start': '10', 'end': '30', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': 'consensus', 'pmid': 'humans'}],
             },

    'd_p_n': { 'humans': {'user1': [{'start': '10', 'end': '20', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': 'humans'}],
                          'user2': [{'start': '15', 'end': '30', 'norm': '11234', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': 'humans'}],
                          'user3': [{'start': '15', 'end': '30', 'norm': '11234', 'annottype': 'e_1', 'text': 'foo', 'user': '3', 'pmid': 'humans'}],
                         },
               'expected': [{'start': '10', 'end': '30', 'norm': '11234', 'annottype': 'e_1', 'text': 'foo', 'user': 'consensus', 'pmid': 'humans'}],
             },
}

fscore_testcases = [
    # test name, document, recall, precision, n_recall, n_precision
    ["proteins", "same_bound_same_norm_100", "100", 1, 1, 1, 1],
    ["proteins", "diff_bound_same_norm_101", "101", 1, 1, 1, 1],
    ["proteins", "same_bound_diff_norm_102", "102", 0, 0, 1, 1],
    #["proteins", "only_in_one", "103_103", 0, 0, 0, 0],
    ["proteins", "diff_bound_diff_norm_104", "104", 0, 0, 1, 1],
    ["proteins", "same_bound_similar_norm_105", "105", 1, 1, 1, 1],
    ["proteins", "same_bound_similar_norm2_106", "106", 1, 1, 1, 1],


    ["species", "same_bound_same_norm_3388767", "3388767", 1, 1, 1, 1],
    ["species", "diff_bound_same_norm_3388768", "3388768", 1, 1, 1, 1],
    ["species", "same_bound_diff_norm_3388769", "3388769", 0, 0, 1, 1],
    #["species", "only_in_one_3388770", "3388770", 0, 0, 0, 0],
    ["species", "diff_bound_diff_norm_3388771", "3388771", 0, 0, 1, 1],
    ["species", "same_bound_similar_norm_3388772", "3388772", 1, 1, 1, 1],
    ["species", "same_bound_similar_norm_adeno_3001", "3001", 1, 1, 1, 1],

]

testcases_boundaries = [
    # same
    [{'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': '1'}, 
     {'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '1'},
     {'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '1'},
     True
     ],

    # end same
    [{'start': '25', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': '2'}, 
     {'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '2'},
     {'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '2'},
     True
     ],

    # start same
    [{'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': '3'}, 
     {'start': '20', 'end': '25', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '3'},
     {'start': '20', 'end': '25', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '3'},
     True
     ],

    # ooodddlll
    [{'start': '20', 'end': '25', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': '4'}, 
     {'start': '23', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '4'},
     {'start': '23', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '4'},
     True
     ],

    # llldddooo
    [{'start': '22', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': '5'}, 
     {'start': '20', 'end': '27', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '5'},
     {'start': '20', 'end': '27', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '5'},
     True
     ],

    # subset
    [{'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': '6'}, 
     {'start': '23', 'end': '25', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '6'},
     {'start': '23', 'end': '25', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '6'},
     True
     ],

    # superset
    [{'start': '22', 'end': '27', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': '7'}, 
     {'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '7'},
     {'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '7'},
     True
     ],

    # no match
    [{'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': '8'}, 
     {'start': '10', 'end': '19', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '8'},
     None,
     False
     ],

    # no match
    [{'start': '20', 'end': '29', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '1', 'pmid': '9'}, 
     {'start': '30', 'end': '39', 'norm': '11676', 'annottype': 'e_1', 'text': 'foo', 'user': '2', 'pmid': '9'},
     None,
     False
     ],

    
]

parse_tagtog = {
        "24933691": [ "The role of matrix in HIV-1 envelope glycoprotein incorporation. Abstract ",
                "Incorporation of the viral envelope (Env) glycoprotein is a critical requirement for the production of infectious HIV-1 particles. It has long been appreciated that the matrix (MA) domain of the Gag polyprotein and the cytoplasmic tail of Env are central players in the process of Env incorporation, but the precise mechanisms have been elusive. Several recent developments have thrown light on the contributions of both proteins, prompting a re-evaluation of the role of MA during Env incorporation. The two domains appear to play distinct but complementary roles, with the cytoplasmic tail of Env responsible for directing Env to the site of assembly and the matrix domain accommodating the cytoplasmic tail of Env in the Gag lattice."],

        "25352622": [ "Human cytomegalovirus-encoded pUL7 is a novel CEACAM1-like molecule responsible for promotion of angiogenesis. Abstract UNLABELLED ",
                "Persistent human cytomegalovirus (HCMV) infection has been linked to several diseases, including atherosclerosis, transplant vascular sclerosis (TVS), restenosis, and glioblastoma. We have previously shown that factors secreted from HCMV-infected cells induce angiogenesis and that this process is due, at least in part, to increased secretion of interleukin-6 (IL-6). In order to identify the HCMV gene(s) responsible for angiogenesis promotion, we constructed a large panel of replication-competent HCMV recombinants. One HCMV recombinant deleted for UL1 to UL10 was unable to induce secretion of factors necessary for angiogenesis. Fine mapping using additional HCMV recombinants identified UL7 as a viral gene required for production of angiogenic factors from HCMV-infected cells. Transient expression of pUL7 induced phosphorylation of STAT3 and ERK1/2 MAP kinases and production of proangiogenic factors, including IL-6. Addition of recombinant pUL7 to cells was sufficient for angiogenesis and was again associated with increased IL-6 expression. Analysis of the UL7 structure revealed a conserved domain similar to the immunoglobulin superfamily domain and related to the N-terminal V-like domain of carcinoembryonic antigen-related cell adhesion molecule 1 (CEACAM1). Our report therefore identifies UL7 as a novel HCMV-encoded molecule that is both structurally and functionally related to cellular CEACAM1, a proangiogenic factor highly expressed during vasculogenesis. IMPORTANCE ",
                "A hallmark of cytomegalovirus (CMV) infection is its ability to modulate the host cellular machinery, resulting in the secretion of factors associated with long-term diseases such as vascular disorders and cancer. We previously demonstrated that HCMV infection alters the types and quantities of bioactive proteins released from cells (designated the HCMV secretome) that are involved in the promotion of angiogenesis and wound healing. A key proangiogenic and antiapoptotic factor identified from a proteomic-based approach was IL-6. In the present report, we show for the first time that HCMV UL7 encodes a soluble molecule that is a structural and functional homologue of the CEACAM1 proangiogenic cellular factor. This report thereby identifies a critical component of the HCMV secretome that may be responsible, at least in part, for the vascular dysregulation associated with persistent HCMV infection."],

        "25821226": ["Ebola Virus Glycoprotein Promotes Enhanced Viral Egress by Preventing Ebola VP40 From Associating With the Host Restriction Factor BST2/Tetherin. Abstract BACKGROUND ",
                "BST2/tetherin is an innate immune molecule with the unique ability to restrict the egress of human immunodeficiency virus (HIV) and other enveloped viruses, including Ebola virus (EBOV). Coincident with this discovery was the finding that the HIV Vpu protein down-regulates BST2 from the cell surface, thereby promoting viral release. Evidence suggests that the EBOV envelope glycoprotein (GP) also counteracts BST2, although the mechanism is unclear. RESULTS ",
                "We find that total levels of BST2 remain unchanged in the presence of GP, whereas surface BST2 is significantly reduced. GP is known to sterically mask surface receptors via its mucin domain. Our evaluation of mutant GP molecules indicate that masking of BST2 by GP is probably responsible for the apparent surface BST2 down-regulation; however, this masking does not explain the observed virus-like particle egress enhancement. We discovered that VP40 coimmunoprecipitates and colocalizes with BST2 in the absence but not in the presence of GP. CONCLUSIONS ",
                "These results suggest that GP may overcome the BST2 restriction by blocking an interaction between VP40 and BST2. Furthermore, we have observed that GP may enhance BST2 incorporation into virus-like particles. Understanding this novel EBOV immune evasion strategy will provide valuable insights into the pathogenicity of this deadly pathogen."],
        }

offset_testcases = {
        "24933691": [ ["s1h1", 0],
                      ["s2p1", 74],
                    ],
        "25352622": [ ["s1h1", 0],
                      ["s2s1p1", len(parse_tagtog["25352622"][0])],
                      ["s2s2p1", len(parse_tagtog["25352622"][1]) + len(parse_tagtog["25352622"][0])],
                    ],
        "25821226": [ ["s1h1", 0],
                      ["s2s1p1", len(parse_tagtog["25821226"][0])],
                      ["s2s2p1", len(parse_tagtog["25821226"][1]) + len(parse_tagtog["25821226"][0])],
                      ["s2s3p1", len(parse_tagtog["25821226"][2]) + len(parse_tagtog["25821226"][1]) + len(parse_tagtog["25821226"][0])],
                    ],
    }

def test_generator_fscore(annottype, document, recall, precision, n_recall, n_precision, uniprot, taxtree):
    def test(self):
        files = ["data.1." + document, "data.2." + document]

        annotations = defaultdict(lambda: defaultdict(list))
        for f in files:
            annotations = parse_annotations(self.data_dir + f, annotations)
        stats = inter_annotator(annotations, uniprot, taxtree)

        self.assertEqual( stats["1"]["2"][document][annottype]["recall"], recall, 'recall')
        self.assertEqual( stats["1"]["2"][document][annottype]["precision"], precision, 'precision')
        self.assertEqual( stats["1"]["2"][document][annottype]["n_recall"], n_recall, 'n_recall' )
        self.assertEqual( stats["1"]["2"][document][annottype]["n_precision"], n_precision, 'n_precision' )
    return test

def test_generator_offset(pmid, para_name, expected):
    def test(self):
        path = "../test_data/tagtog/"
        offset = get_position_offset(pmid, para_name, path)
        self.assertEqual(offset, expected)
    return test

def test_generator_consensus(case, uniprot, taxtree):
    def test(self):
        actual = human_consensus({'humans': case['humans']}, uniprot, taxtree)
        self.assertEqual(actual['humans'], case['expected'])
    return test

class TestIAA(unittest.TestCase):
    data_dir = "../test_data/"

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # tests are generated by test_generator_*() and added here at runtime

    # extra tests that aren't generated by the generators
    def test_parse_tagtog(self):
        path = "../test_data/tagtog"

        for pmid, expected in parse_tagtog.iteritems():
            documents = parse_tagtog_document(pmid, path)
            for i in range(0, len(documents)):
                self.assertEqual(documents[i], expected[i])

    def test_tagger_parsing(self):
        annotations = defaultdict(lambda: defaultdict(list))

        inputfile = "../test_data/tagtog/data.test2.azQqteWeLBbFAwJOP4HRFDEwv1l8-24933691.tsv"
        annotations = parse_annotations(inputfile, annotations)

        tagger_output_file = "../../tagtog/tagger_output/virus.matches"
        entities_file = "../../tagtog/tagger_output/virus_entities.tsv"
        annotations = parse_tagger(tagger_output_file, entities_file, annotations)

        tagger_1 = annotations["24933691"]['tagger'][0]
        tagger_2 = annotations["24933691"]['tagger'][2]

        helen_1 = annotations["24933691"]['Helen'][1]
        helen_2 = annotations["24933691"]['Helen'][6]

        self.assertEqual(tagger_1['start'], helen_1['start'])
        self.assertEqual(tagger_1['end'], helen_1['end'])
        self.assertEqual(tagger_2['start'], helen_2['start'])
        self.assertEqual(tagger_2['end'], helen_2['end'])
        
    def test_boundaries(self):
        for t in testcases_boundaries:
            an = t[0]
            alist = [t[1]]
            expected = [t[2]]
            if t[2] is None:
                expected = []
            actual = same_boundaries(an, alist)
            self.assertEqual(actual, expected)
            


if __name__ == '__main__':
    filename_uniprot = "../../../uniprot_to_payload/data_in/uniprot_viruses.xml"
    filename_taxonomy =  "../../../taxonomy/nodes.dmp"
    uniprot = parse_uniprot_single_xml(filename_uniprot)
    taxtree = tax.parse_taxtree(filename_taxonomy)

    for t in fscore_testcases:
        annottype, testname, document, recall, precision, n_recall, n_precision = t
        test_name = 'test_fscore_' + annottype + '_%s' % "_".join(testname.split(" "))
        test = test_generator_fscore(annottype, document, recall, precision, n_recall, n_precision, uniprot, taxtree)
        setattr(TestIAA, test_name, test) # assign the correct name to the test

    for t in offset_testcases.iteritems():
        pmid, cases = t
        for c in cases:
            para_name, expected = c
            test_name = 'test_offset_' + pmid + '_' + para_name
            test = test_generator_offset(pmid, para_name, expected)
            setattr(TestIAA, test_name, test)

    for t in consensus_testcases.iteritems():
        name, case = t
        test_name = 'test_consensus_' + name
        test = test_generator_consensus(case, uniprot, taxtree)
        setattr(TestIAA, test_name, test)

    unittest.main()

suite = unittest.TestLoader().loadTestsFromTestCase(TestIAA)
unittest.TextTestRunner(verbosity=2).run(suite)

