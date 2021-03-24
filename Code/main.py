#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# Algorithm for combined retrieve of proteomes and
# prediction of their binding epitopes.
# 24.03.2021: Product.
#################################################


import argparse
import time
import uniProt
import NetMHCpan as Nmp
import Syfpeithi as Syfp
import utilities as ut


start_o = time.time()
parser = argparse.ArgumentParser(description="Collection of proteome data and associated predicted epitopes.")
parser.add_argument('-f', '--full', help="Executes the algorithm for individual data. You will be asked to input the"
                                         " parameters for organism, location(s) (split by comma), status of review and"
                                         " output format of the UniProt data. Organism parameter must exactly occure in"
                                         " the last entry of taxonomy in UniProt. E.x. choose '(strain USA300)' instead"
                                         " of 'strain USA300' if the strain USA300 / TH1516is not wanted. Locations are"
                                         " optional to set. Only reviewed entries are certainly reliable. Tab format"
                                         " will collect proteom data and epitope predictions whereas fasta format only"
                                         " produces the .fsa files of the proteomes. At the end of protein retrieval"
                                         " deleting redundant data for the same protein but different strains is asked."
                                         " The allele parameter must exactly match the international nomenclature until"
                                         " allele position like HLA-A*02:01. Epitope length is restrict by the"
                                         " prediction algorithms, look it up beforehand SYFPEITHI's score ranges from"
                                         " 0 up to a maximum value that is dependent on the allel. For HLA-A*02:01 e.x."
                                         " it is 36. NetMHCpan's score is the percentile rank of an epitope compared to"
                                         " a set of random natural peptides and thus is affected by bias and ranges"
                                         " from 0 to 100. At the end of each step a spreadsheet file named by you or"
                                         " standard will be produced.",
                    action="store_true",
                    )
parser.add_argument('-e', '--example', help="Executes the algorithm for Organism: Staphylococcus Aureus; Locations:"
                                            " cell membrane, cell wall, cytoplasm and secreted; Reviewed: yes; Format:"
                                            " tab; Delete redundant: yes; Allele: HLA-A*02:01; Length: 9; SYFPEITHI"
                                            " theshold: 20; NetMHCpan thershold: 3.0; Output names: example_sa_uniprot"
                                            ".xlsx, example_sa_syfpeithi.xlsx, example_sa_netmhcpan.xlsx.",
                    action="store_true")
args = parser.parse_args()

# Test if argument was specified.
if args.example or args.full:
    if args.example:
        example = True
    else:
        example = False

    # Execute retrieval of proteome data from uniprot
    protein_data = uniProt.exec_uniprot(example)

    # Set prediction parameters
    params = ut.pred_params(example)

    # Execute epitope prediction by SYFPEITHI
    start = time.time()
    syfp_tmp = Syfp.exec_syfp(protein_data, example, params)
    end = time.time()
    print("SYFPEITHI done in " + ut.calc_time(end-start))

    # Execute epitope prediction by NetMHCpan
    start = time.time()
    nmp_tmp = Nmp.exec_netmhcpan(protein_data, example, params)
    end = time.time()
    print("NetMHCpan done  in " + ut.calc_time(end-start))

    # Print due to redundance deleted sequences
    ut.print_popped(protein_data[2])

    end_o = time.time()
    print("Overall done in: " + ut.calc_time(end_o-start_o))
else:
    exit("Please see --help [-h] and specify the argument!")
