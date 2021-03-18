#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# Algorithm for combined retrieve of proteins and
# prediction of their binding epitopes.
# 17.03.2021: Development.
#################################################


import time
import operator
import requests as req
import utilities as ut
from collections import Counter


# To fill data dict with given strings
def fill_data(s, a, length):
    data = {
        'method': 'netmhcpan',
        'sequence_text': s,
        'allele': a,
        'length': length
    }
    return data


def exec_netmhcpan(uniprot, example, params):
    # Dict with sequences
    seq_dict = uniprot[0]
    # Dict with signal peptides to look if epitope is in signal range.
    sig_dict = uniprot[1]
    # Dict with keys that have the same sequence
    same_seq = uniprot[3]
    #
    redundant_del = uniprot[4]
    # Dict with indices of epitopes located in the signal peptide
    in_signal = {}
    # List with epitopes for counting the frequency.
    epi_list = []
    epi_data = {}
    # Storage for whole prediction data
    data = {}
    print("NetMHCpan:")
    # Set parameters.
    allele, length, threshold, dumb = params
    del dumb
    params = {"Allel": allele, "Length": length, "Threshold": str(threshold)}

    header = 'id,protein,density,sig density,position; epitope; score'
    count_seqs = sum(len(seqs) for seqs in seq_dict.values())
    count_same_seqs = sum(len(seqs) for seqs in same_seq.values())
    if not redundant_del:
        count_all_predictions = count_seqs - count_same_seqs
    else:
        count_all_predictions = count_seqs
    count = 0
    start = time.time()

    for loc_key in seq_dict:
        epi_data[loc_key] = {}
        in_signal[loc_key] = {}

        for gene_id in seq_dict[loc_key]:
            epi_data[loc_key][gene_id] = []
            temp_epi_list = []
            temp_epis = []
            in_signal[loc_key][gene_id] = []

            if gene_id in same_seq[loc_key]:
                equal_seq_id = same_seq[loc_key][gene_id]
                if data[equal_seq_id]:
                    epi_list.extend(data[equal_seq_id][0])
                    epi_data[loc_key][gene_id] = data[equal_seq_id][1]
                    in_signal[loc_key][gene_id].extend(data[equal_seq_id][2])
                    count += 1
                    ut.progress(count, count_seqs, count_all_predictions, start)
                    continue

            # Request.
            seq = seq_dict[loc_key][gene_id][0]
            pred_params = fill_data(seq, allele, length)
            try:
                response = req.post('http://tools-cluster-interface.iedb.org/tools_api/mhci/', data=pred_params)
            except ConnectionError:
                print("Connection Error at: " + gene_id)
                break

            # Information.
            protein_name = seq_dict[loc_key][gene_id][1]
            seq_len = len(seq)

            if response.ok:
                # Remove header
                lines = response.text.splitlines()[1:]

                # Filter important information.
                for line in lines:
                    pred_data = line.split("\t")
                    if len(pred_data) >= 8:
                        position, epi, score = pred_data[2], pred_data[5], pred_data[9]
                    else:
                        break

                    # Only take top scoring epi_datas. Threshold is 3.0.
                    if float(score) <= threshold:
                        # To list index of epitopes which are located in the signal peptide of the protein.
                        # Set signal range for density calc.
                        if sig_dict[loc_key][gene_id].isnumeric() and int(position) <= int(sig_dict[loc_key][gene_id]):
                            in_signal[loc_key][gene_id].append(lines.index(line)+3)

                        # Add filtered information to list and dictionary about epitopes.
                        temp_epi_list.append(epi)
                        temp_epis.append(position + "; " + epi + "; " + score)

                # Append protein name for information.
                epi_data[loc_key][gene_id].append(protein_name)

                # Calculate epitope density.
                # Overall.
                epitopes_bound = len(temp_epis)
                epitopes_total = seq_len - int(length) + 1
                if epitopes_total > 0:
                    density = epitopes_bound / epitopes_total
                else:
                    density = 0
                # Store density with epitope data.
                epi_data[loc_key][gene_id].append(density)

                # In signal.
                if sig_dict[loc_key][gene_id].isnumeric():
                    # If signal is given as number, not 'N/A'.
                    signal_length = sig_dict[loc_key][gene_id]
                    signal_epitopes_bound = len(in_signal[loc_key][gene_id])
                    signal_epitopes_total = int(signal_length) - int(length)
                    if signal_epitopes_total > 0:
                        signal_density = signal_epitopes_bound / signal_epitopes_total
                    else: signal_density = 0
                    # Store signal density with epitope data.
                    epi_data[loc_key][gene_id].append(signal_density)
                else:
                    epi_data[loc_key][gene_id].append("N/A")

                epi_data[loc_key][gene_id].extend(temp_epis)
                epi_list.extend(temp_epi_list)
                data[gene_id] = [temp_epi_list, epi_data[loc_key][gene_id], in_signal[loc_key][gene_id]]

            else:
                print(loc_key, gene_id, "\tBad Answer!")
                print(response.text)
                epi_data[loc_key][gene_id] = [protein_name, "N/A", "N/A", "BAD ANSWER."]
                data[gene_id] = None

            count += 1
            ut.progress(count, count_seqs, count_all_predictions, start)

    # Counted appearance of every epitope.
    frequencies = dict(sorted(Counter(epi_list).items(), key=operator.itemgetter(1), reverse=True))

    # Filename input
    if example:
        filename = "example_sa_netmhcpan.xlsx"
    else:
        name = input("\nLeave clear for standard 'netMHCPan_data.xlsx' or enter custom filename: ")
        if name:
            filename = name
        else:
            filename = "netMHCpan_data.xlsx"

    ut.print_sheet(header, params, epi_data, frequencies, filename, "NetMHCPan-4.0", in_signal)
    return epi_data
