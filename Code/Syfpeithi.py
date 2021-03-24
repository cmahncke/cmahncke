#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# Algorithm for combined retrieve of proteomes and
# prediction of their binding epitopes.
# 24.03.2021: Product.
#################################################


import time
import operator
import requests as req
import utilities as ut
from collections import Counter


# Main function to execute prediction by SYFPEITHI.
def exec_syfp(uniprot, example, params):
    print("SYFPEITHI:")
    # Dict containing sequences.
    seq_dict = uniprot[0]
    # Dict containing signal peptide ranges.
    sig_dict = uniprot[1]
    # Dict containing keys that have the same sequence.
    same_seq = uniprot[3]
    # Boolean about deletion of redundant sequences.
    redundant_del = uniprot[4]
    # Dict with indices of epitopes located in the signal peptide.
    in_signal = {}
    # List with epitopes for counting the frequency.
    epi_list = []
    # Storage for prediction data by location.
    epi_data = {}
    # Storage for whole prediction data.
    data = {}
    # Set parameters.
    allele, length, dumb, threshold = params
    del dumb
    params = {"Allel": allele, "Length": length, "Threshold": str(threshold)}

    # Calculate amount of predictions for providing progress time and percentage.
    # Prediction for redundant sequences does not affect progress time since their data is just copied.
    count_seqs = sum(len(seqs) for seqs in seq_dict.values())
    count_same_seqs = sum(len(seqs) for seqs in same_seq.values())
    if not redundant_del:
        count_all_predictions = count_seqs - count_same_seqs
    else:
        count_all_predictions = count_seqs
    count = 0
    start = time.time()

    # Iterate over locations.
    for loc_key in seq_dict:
        epi_data[loc_key] = {}
        in_signal[loc_key] = {}

        # Iterate over IDs of the proteins.
        for gene_id in seq_dict[loc_key]:
            epi_data[loc_key][gene_id] = []
            in_signal[loc_key][gene_id] = []
            protein_name = seq_dict[loc_key][gene_id][1]

            # Check if an ID leading to an equal sequence exists.
            # If so, paste its prediction data to this entry.
            if gene_id in same_seq[loc_key]:
                equal_seq_id = same_seq[loc_key][gene_id]
                if data[equal_seq_id]:
                    epi_list.extend(data[equal_seq_id][0])
                    epi_data[loc_key][gene_id] = data[equal_seq_id][1]
                    in_signal[loc_key][gene_id].extend(data[equal_seq_id][2])
                    count += 1
                    ut.progress(count, count_seqs, count_all_predictions, start)
                    continue

            # Cut the sequences at AA 2000 due to SYFPEITHI webserver input limit.
            seq = seq_dict[loc_key][gene_id][0]
            too_long = False
            if len(seq) > 2000:
                seq = seq[:1999]
                too_long = True

            # Fill the URL with parameters POST it for prediction.
            url = 'http://www.syfpeithi.de/bin/MHCServer.dll/EpitopePrediction?Motif=' + allele + '&amers=' \
                  + length + '&SEQU=' + seq + '&DoIT=++Run++'
            try:
                response = req.post(url)
            except ConnectionError:
                print("Connection Error at: " + gene_id)
                break

            if response.ok:
                # Filter important data, calculate densities and add everything to storages.
                ut.request_handle(method="SYFPEITHI", request=response, indices=[loc_key, gene_id],
                                  pred_vars=[len(seq), length, threshold, too_long],
                                  epi_storages=[epi_list, epi_data, data], signal_storages=[sig_dict, in_signal],
                                  protein_name=protein_name)

            else:
                # Print note at console and output file if the request is not okay.
                print(loc_key, gene_id, "\tBad Answer.")
                print(response.text)
                epi_data[loc_key][gene_id] = [protein_name, "N/A", "N/A", "BAD ANSWER."]
                data[gene_id] = None

            count += 1
            ut.progress(count, count_seqs, count_all_predictions, start)

    # Count appearences of each epitope.
    frequencies = dict(sorted(Counter(epi_list).items(), key=operator.itemgetter(1), reverse=True))

    # Filename input
    if example:
        filename = "example_sa_syfpeithi.xlsx"
    else:
        name = input("\nLeave clear for standard 'SYFPEITHI_data.xlsx' or enter custom filename: ")
        if name:
            filename = name
        else:
            filename = "SYFPEITHI_data.xlsx"

    header = 'id,protein,density,sig density,position; epitope; score'
    ut.print_sheet(header, params, epi_data, frequencies, filename, "SYFPEITHI", in_signal)
    return epi_data
