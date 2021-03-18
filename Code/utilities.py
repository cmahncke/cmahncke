#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# Algorithm for combined retrieve of proteins and
# prediction of their binding epitopes.
# 11.02.2021: Development.
#################################################

import xlsxwriter
import numpy as np
import re
import time
import sys


def print_file(text, file):
    with open(file, 'w') as fi:
        fi.write(text)
    fi.close()


def extend_file(text, file):
    with open(file, 'a') as f:
        f.write(text)
    f.close()


# To write requested data into xlsx files.
def print_sheet(header, params, str_dict, freq_dict, filename, method, *signal_mark):
    if signal_mark:
        in_signal = signal_mark[0]
    else:
        in_signal = None

    # Create workbook and add first sheet for overview.
    wb = xlsxwriter.Workbook(filename)
    sheet = wb.add_worksheet("method")
    form = wb.add_format()

    # Overview.
    form.set_bold()
    sheet.write(0, 0, method, form)
    # Parameters that has been set.
    sheet.write(2, 0, "Parameter:", form)
    for param in params:
        col_index = list(params.keys()).index(param)
        sheet.write(1, col_index + 1, param, form)
        if params[param].isnumeric():
            sheet.write(2, col_index + 1, int(params[param]), form.set_num_format("0"))
        else:
            sheet.write(2, col_index + 1, params[param])

    sheet.write(4, 0, "Frequencies over all loci:", form)
    for gene_id in freq_dict:
        row_index = list(freq_dict.keys()).index(gene_id)
        sheet.write(row_index + 5, 0, gene_id)
        if isinstance(freq_dict[gene_id], int):
            sheet.write(row_index + 5, 1, int(freq_dict[gene_id]))
        else:
            sheet.write(row_index + 5, 1, freq_dict[gene_id])

    # General data.
    base = "http://www.uniprot.org/uniprot/"
    header = header.split(",")
    header.insert(-1, "link")

    # For every location there is a sheet added.
    for loc_key in str_dict:
        sheet = wb.add_worksheet(loc_key)
        form = wb.add_format()

        # The header is written on each sheets top.
        for col_index in range(len(header)):
            # Filter the argument from columns like feature(SIGNAL) or linage(PHYLUM).
            if re.findall(r'\(', header[col_index]):
                if re.findall(r'\((.*)\)', header[col_index]):
                    string = re.findall(r'\((.*)\)', header[col_index])[0]
                    if string.lower() != "signal":
                        string = "lineage"
            else:
                string = header[col_index]
            form.set_bold()
            sheet.write(0, col_index, string.upper(), form)

        # Calculate 3rd quantile of epitope density.
        if method.lower() != "uniprot":
            try:
                dens_string = list(list(zip(*list(str_dict[loc_key].values())))[1])
                densities = remove_nas(dens_string)
                q3_density = np.quantile(densities, .75)
            except IndexError:
                print("No Densities")
        else:
            q3_density = None

        # For every gene in the current location the ID marks the first column.
        for gene_id in str_dict[loc_key]:
            form = wb.add_format()
            row_index = list(str_dict[loc_key].keys()).index(gene_id)

            # Color ID if gene density is above Q3 of epitope density.
            if q3_density:
                try:
                    if float(str_dict[loc_key][gene_id][1]) > q3_density:
                        form.set_bg_color('#00FF00')
                except (TypeError, ValueError) as err:
                    print(gene_id + " No Density!")
                    print(err)

            sheet.write(row_index + 1, 0, gene_id, form)

            # Color sequenze note in red
            if str_dict[loc_key][gene_id][-1] == "Sequence has been cut at pos 2000 due to SYFPEITHI online limit.":
                form.set_bg_color('#FF4500')

            # Subsequently the columns are filled with data stored at gene ID.
            result_len = len(str_dict[loc_key][gene_id])
            link = 0
            for col_index in range(result_len):
                form = wb.add_format()
                if method.lower() == "uniprot":

                    # Protein names in italic.
                    if header[col_index + 1] == 'protein_names':
                        form.set_italic()
                    # Numbers in number format.
                    if header[col_index + 1] == ('length' or 'mass' or 'feature(SIGNAL)'):
                        form.set_num_format('0')
                    # Convert mass to number.
                    if header[col_index + 1] == 'mass':
                        str_dict[loc_key][gene_id][col_index] = re.sub(",", "", str_dict[loc_key][gene_id][col_index])
                        form.set_num_format('0')
                    # Print only last lineage entry, not whole lineage
                    if header[col_index + 1] == 'lineage(ALL)':
                        str_dict[loc_key][gene_id][col_index] = str_dict[loc_key][gene_id][col_index].split(",")[-1]
                    # Write Link at the end of each row in link format
                    if header[col_index + 1] == 'link':
                        sheet.write_url(row_index + 1, col_index + 1, base + gene_id, string='UniProt')
                        link += 1
                else:
                    # Protein names in italic.
                    if str_dict[loc_key][gene_id][-1] == "Sequence was cut at pos 2000 due to SYFPEITHI Webserver limit.":
                        if col_index in {0, 1, 2}:
                            form.set_bg_color('#FF4500')
                        if col_index == len(str_dict[loc_key][gene_id]) - 1:
                            form.set_font_color('red')
                    if col_index == 0:
                        form.set_italic()
                    # Numbers in number format.
                    if col_index in {1, 2}:
                        form.set_num_format('0.##0')
                    # Add border to epitopes
                    # Write link in link format and add border to epitope data
                    if col_index == 3:
                        form.set_left()
                        sheet.write_url(row_index + 1, col_index + 1, base + gene_id, string='UniProt')
                        link += 1
                    if col_index == 2 and result_len == 3:
                        sheet.write_url(row_index + 1, col_index + 2, base + gene_id, string='UniProt')
                    # Mark the epitopes inside signal.
                    if in_signal and col_index in in_signal[loc_key][gene_id]:
                        form.set_font_color('green')

                if isinstance(str_dict[loc_key][gene_id][col_index], int):
                    sheet.write(row_index + 1, col_index + 1 + link, int(str_dict[loc_key][gene_id][col_index]), form)
                elif isinstance(str_dict[loc_key][gene_id][col_index], float):
                    sheet.write(row_index + 1, col_index + 1 + link, float(str_dict[loc_key][gene_id][col_index]), form)
                elif str_dict[loc_key][gene_id][col_index].isnumeric():
                    sheet.write(row_index + 1, col_index + 1 + link, int(str_dict[loc_key][gene_id][col_index]), form)
                else:
                    sheet.write(row_index + 1, col_index + 1 + link, str_dict[loc_key][gene_id][col_index], form)

            # At the end of each row there is an url guiding to the UniProt page of the gene
    wb.close()


# Deletes redundant sequences from two dictionaries. In particular, the one printed to spreadsheets and the one given
# to the predition algorithms.
def del_redundant(data_dict, redundant_ids):
    popped = {}

    for key in list(data_dict):
        if key in redundant_ids:
            sequence = data_dict[key][-1]
            popped[key] = sequence
            data_dict.pop(key)
    return popped


def pred_params(bool):
    if bool:
        allel = "HLA-A*02:01"
        length = "9"
        nmp = 3
        spt = 20
    else:
        allel = input("Allel: ")
        length = input("Length: ")
        print("Scores: highest -> lowest")
        spt = int(input("SYFPEITHI-Threshold [30 -> 0]: "))
        nmp = float(input("NetMHCpan-Threshold [0.0 -> 100.0]: "))
    return allel, length, nmp, spt


def calc_time(passed_seconds):
    seconds = passed_seconds % 60
    minutes = passed_seconds / 60 % 60
    hours = passed_seconds / 3600
    time = "%02d:%02d:%02d" % (hours, minutes, seconds)
    return time


def print_popped(seq_dict):
    answer = input("Print popped sequences? [yes/no] ")
    if answer.lower() == "yes":
        for loc_key in seq_dict:
            print("\n" + loc_key + "\n")
            for gene_id in seq_dict[loc_key]:
                print(gene_id)
                print("\n" + seq_dict[loc_key][gene_id] + "\n")
        return 0
    else:
        return 0


def remove_nas(str_list):
    return [element for element in str_list if element != 'N/A']


def if_remove_redundant(example):
    if example:
        answer = 'yes'
    else:
        answer = input("Delete redundant sequences? [yes/no]: ")
    if answer.lower() == 'yes':
        return True
    else:
        return False


# Provide estimated time and percentage of completed tasks
def progress(count, count_all, count_all_preds, start_time):
    percentage_old = (count-1) / count_all * 100
    if count == 3:
        end = time.time()
        estimated_seconds = (end - start_time) * count_all_preds / 3
        print("\nEstimated time:", calc_time(estimated_seconds))

    percentage_new = count / count_all * 100
    if int(percentage_new) > int(percentage_old):
        sys.stdout.write("\r{0}%".format(int(percentage_new)))
        sys.stdout.flush()
