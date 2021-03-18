#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# Algorithm for combined retrieve of proteins and
# prediction of their binding epitopes.
# 11.02.2021: Development.
#################################################


import re
import time
import operator
import requests as req
import utilities as ut
from collections import Counter


# To ask for the values to setup SYFPEITHI.
def input_data():
    allel = input("Allele: ")
    length = input("Length: ")
    return allel, length


# Response from SYFPEITHI is just HTML file. This function filters the relevant data from HTML code and splits it into
# groups for subsequent use.
def request_handle(request, threshold, coords, protein_name, len_data, sig_dict):
    # clean the request test from html fragments.
    stripped = re.sub('<[^<]+?>|&nbsp;|\r|[ ]', '', request.text)
    # only use the top 10.
    lines = stripped.splitlines()[2:-1]

    # Rename variables for homogeneity
    loc_key, gene_id = coords
    length, seq_len = len_data

    epi_list = []
    epi_info = []
    temp_epis = []
    in_signal_list = []

    for line in lines:
        # Split the data in position, sequence and score.
        groups = re.match(r'([0-9]+)([A-Z]+)(-?[0-9]+)', line, re.I).groups()

        # Only take top scoring epi_datas. Threshold is 20.
        if int(groups[2]) > threshold:
            # To ist index of epitopes which are located in the signal peptide of the protein.
            # Set signal range for density calc.
            position = groups[0]
            if sig_dict[loc_key][gene_id].isnumeric() and int(position) <= int(sig_dict[loc_key][gene_id]):
                in_signal_list.append(lines.index(line)+3)

            # Add filtered information to list and dictionary about epitopes.
            epi_list.append(groups[1])
            temp_epis.append("; ".join(groups))

    # Append protein name for information.
    epi_info.append(protein_name)

    # Calculate epitope density.
    # Overall.
    epitopes_bound = len(epi_list)
    epitopes_total = seq_len - int(length) +1
    if epitopes_total > 0:
        density = epitopes_bound / epitopes_total
    else:
        density = 0
    # Store density with epitope data.
    epi_info.append(density)

    # In signal.
    if sig_dict[loc_key][gene_id].isnumeric():
        # If signal is given as number, not 'N/A'.
        signal_length = sig_dict[loc_key][gene_id]
        signal_epitopes_bound = len(in_signal_list)
        signal_epitopes_total = int(signal_length) - int(length)
        if signal_epitopes_total > 0:
            signal_density = signal_epitopes_bound / signal_epitopes_total
        else:
            signal_density = 0
        # Store signal density with epitope data.
        epi_info.append(signal_density)
    else:
        epi_info.append("N/A")

    epi_info.extend(temp_epis)

    return epi_list, epi_info, in_signal_list


# Main function for executing SYFPEITHI on the data from uniprot.
def exec_syfp(uniprot, example, params):
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
    print("SYFPEITHI:")
    # Set parameters.
    allele, length, dumb, threshold = params
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
            in_signal[loc_key][gene_id] = []

            if gene_id in same_seq[loc_key]:
                equal_seq_id = same_seq[loc_key][gene_id]
                epi_list.extend(data[equal_seq_id][0])
                epi_data[loc_key][gene_id] = data[equal_seq_id][1]
                in_signal[loc_key][gene_id].extend(data[equal_seq_id][2])
                count += 1
                ut.progress(count, count_seqs, count_all_predictions, start)
                continue

            # Request.
            sequence_text = seq_dict[loc_key][gene_id][0]
            too_long = False
            if len(sequence_text) > 2000:
                sequence_text = sequence_text[:1999]
                too_long = True

            url = 'http://www.syfpeithi.de/bin/MHCServer.dll/EpitopePrediction?Motif=' + allele + '&amers=' \
                  + length + '&SEQU=' + sequence_text + '&DoIT=++Run++'
            try:
                response = req.post(url)
            except ConnectionError:
                print("Connection Error at: " + gene_id)
                break

            # Information.
            protein_name = seq_dict[loc_key][gene_id][1]

            if response.ok:
                # Only take top 10 epi_datas and filter important information.
                clean = request_handle(response, threshold, [loc_key, gene_id], protein_name,
                                       [length, len(sequence_text)], sig_dict)
                if too_long:
                    clean[1].append("Sequence was cut at pos 2000 due to SYFPEITHI Webserver limit.")
                # Store filtered information.
                data[gene_id] = clean
                epi_list.extend(clean[0])
                epi_data[loc_key][gene_id] = clean[1]
                in_signal[loc_key][gene_id].extend(clean[2])

            else:
                print(loc_key, gene_id, "\tBad Answer.")
                epi_data[loc_key][gene_id] = [protein_name, "N/A", "N/A", "BAD ANSWER."]

            count += 1
            ut.progress(count, count_seqs, count_all_predictions, start)

    # Counted appearences of every epitope.
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

    ut.print_sheet(header, params, epi_data, frequencies, filename, "SYFPEITHI", in_signal)
    return epi_data
