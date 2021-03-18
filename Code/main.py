#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# Algorithm for combined retrieve of proteins and
# prediction of their binding epitopes.
# 11.02.2021: Development.
#################################################


import argparse
import time
import uniProt
import NetMHCpan as Nmp
import Syfpeithi as Syfp
import utilities as ut


start_o = time.time()
parser = argparse.ArgumentParser(description="Collect peptide sequences at locations for one organism.")
parser.add_argument('-f', '--full', help="Executes the full algorithm.\nYou will be asked for the organism, "
                                         "location(s), status of review and output format for UniProt data. 'tab'"
                                         "format will produce an overview over all proteins and also execute the "
                                         "predictions. 'fasta' only produces the .fsa files.\n"
                                         "Subsequently, unless format is 'fasta', for every peptide the epitopes with "
                                         "top 10 binding affinity will be predicted by NetMHCPan_4.0 and SYFPEITHI.\n"
                                         "Therefore you will be asked about the allel and the epitope length.\n"
                                         "Results will be written in a spreadsheet for each algorithm.",
                    action="store_true")
parser.add_argument('-e', '--example', help="Executes the algorithm for Staphylococcus Aureus at Cell Wall, "
                                            "Cytoplasm and Secreted, only reviewed in tab format.\n"
                                            "The epitopes are predicted for HLA-A*02:01 and a length of 9 amino acids."
                                            "\n The outputs are example_SA_uniprot.xlsx, example_SA_netMHCPan.xlsx"
                                            " and example_SA_SYFPEITHI.xlsx.",
                    action="store_true")
args = parser.parse_args()

if args.example or args.full:
    if args.example:
        example = True
    else:
        example = False

    protein_data = uniProt.exec_uniprot(example)
    if isinstance(protein_data, tuple):
        params = ut.pred_params(example)
        start = time.time()
        syfp_tmp = Syfp.exec_syfp(protein_data, example, params)
        end = time.time()
        print("SYFPEITHI done in " + ut.calc_time(end-start))
        start = time.time()
        nmp_tmp = Nmp.exec_netmhcpan(protein_data, example, params)
        end = time.time()
        print("NetMHCpan done  in " + ut.calc_time(end-start))

        ut.print_popped(protein_data[2])

    else:
        print("Error in UniProt.")

    end_o = time.time()
    print("Overall done in: " + ut.calc_time(end_o-start_o))
else:
    exit("Please see --help [-h] and specify the argument!")
