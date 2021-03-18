#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# Algorithm for combined retrieve of proteins and
# prediction of their binding epitopes.
# 11.02.2021: Development.
#################################################


import requests as req
import utilities as ut
from collections import Counter
import operator


# To ask for the arguments via the console.
def make_query_input():
    try:
        org = input("\nOrganism: ")
        locs = input("Locations: ").split()
        rev = input("Reviewed [yes/no]: ")
        queries = {}

        # Dictionary of queries is keyed by the locations.
        if locs:
            locs_str = ",".join(locs)
            for loc in locs:
                queries[loc] = 'organism:"' + org + '" locations:(location:"' + loc + '") AND reviewed:' + rev
        else:
            queries["All"] = 'organism:"' + org + '" AND reviewed:' + rev
            locs_str = "All"

        params = {'Organism': org, 'Locations': locs_str, 'Reviewed': rev}

        return queries, params
    except IOError:
        print("IOError occurred in query input! ...try again.")


# To make the query for the request by given arguments.
def make_query_args(org, locs, rev):
    try:
        # Dictionary of queries is keyed by the locations.
        queries = {}
        if isinstance(locs, list):
            locs_str = ",".join(locs)
            for loc_key in locs:
                queries[loc_key] = 'organism:"' + org + '" locations:(location:"' + loc_key + '") AND reviewed:' + rev

        # If there is only one location and it is not given as list.
        else:
            locs_str = locs
            queries = 'organism:' + org + ' AND locations:(location:' + locs + ') AND reviewed:' + rev

        params = {'Organism': org, 'Locations': locs_str, 'Reviewed': rev}
        return queries, params
    except ValueError:
        print("Value error occurred when making query! ...try again.")


# To ask for the output format.
def form_input():
    try:
        #print("'tab' format will create an overview and subsequently execute the prediction\n"
        #      "'fasta' produces .fsa FASTA-files without following prediction.")
        form = input("Output format: ")
        return form
    except IOError:
        print("IOError occurred in form input! ...try again.")


# To ask for the output columns.
def columns_input():
    try:
        standard = "genes,length,mass,lineage(ALL)"
        #print("Standard columns: 'gene name'\t'length'\t'mass'\t'lineage'\n"
        #      "ID and protein names are attached at the beginning, "
        #      "signal pepides and sequences at the end automatically.\n"
        #      "For more and how to enter see https://www.uniprot.org/help/uniprotkb_column_names\n"
        #      "Leave clear for standard or type like 'id,sequence,...' without whitespace.")
        cols = None
        #cols = input("Columns [Enter for standard]: ")

        if cols:
            columns = cols
        else:
            columns = standard

        return columns
    except IOError:
        print("IOError occurred in columns input! ...try again.")


# Based on the entries made, the function does the request from uniprot.org.
def search_uniprot(queries, columns, form, params):
    try:
        base = 'http://www.uniprot.org'
        resource = '/uniprot/'
        result = {}
        same_seq = {}

        # For every location there is a request loaded into a dictionary keyed by locations.
        for loc_key in queries:

            # For storing keys to redundant sequences
            sequences = []
            seq_key = {}

            payload = {'query': queries[loc_key],
                       'format': form,
                       'columns': columns}
            request = req.get(base + resource, params=payload)
            result[loc_key] = {}
            same_seq[loc_key] = {}

            if request.ok:
                if form == "tab":
                    # In tab format the data is stored in a nested dictionary with the first column, id in
                    # example, as secondary key.
                    # Delete header and tail of the request.
                    entries = request.text.splitlines()[1:]
                    for entry in entries:
                        split = entry.split("\t")

                        # Filter not wanted strains
                        #strain = split[-3].split(",")[-1]
                        #if params['Organism'] not in strain:
                        #    continue

                        # Store keys with equal sequences
                        seq = split[-1]
                        gene_id = split[0]
                        if seq in sequences:
                            same_seq[loc_key][gene_id] = seq_key[seq]
                        else:
                            sequences.append(seq)
                            seq_key[seq] = gene_id

                        # Filter signal range
                        signal_evidence = split[-2].split(";")
                        if signal_evidence != ['']:
                            signal = signal_evidence[0].split()[1]
                            signal_evidence[0] = signal
                            split[-2] = signal.split(".")[-1]
                        else: split[-2] = "NA"

                        # Nested dictionary with locations as primary keys and the first column, e.g. id, as secondary.
                        result[loc_key][gene_id] = split[1:]

                else:
                    result[loc_key] = request.text

            else:
                print("Error: Request not okay!")
                return -1

        return result, same_seq
    except ConnectionError:
        print("Connection error occurred in uniprot request function! ...try again.")


# To execute the input, search and handle the request for the output
def exec_uniprot(example):
    if example:

        # Example for Staphylococcus aureus at cell_wall, cytoplasm and secreted. Only reviewed peptides.
        # Output format is 'tab' and the columns in output are entry id, gene name, length, mass, signal
        # peptide and amino acid sequence.
        q, params = make_query_args('Staphylococcus aureus', ['cell_membrane', 'cell_wall', 'cytoplasm','secreted'], 'yes')
        f = "tab"
        c = "genes,length,mass,lineage(ALL)"

    else:
        # Universal algorithm for asking the customer about query, form and column entries.
        q, params = make_query_input()
        f = form_input()
        # Columns are only used in tab format, not in fasta.
        if f == "tab":
            c = columns_input()

    # Data will be presented in easy tab format.
    if f == "tab":
        columns = "id,protein_names," + c + ",feature(SIGNAL),sequence"
        data, same_seq = search_uniprot(q, columns, f, params)

        # Container for passing sequences and signal peptides to NetMHCpan and SYFPEITHI (*_dict) and for
        # counting the sequences to give an overview in the output.
        seq_dict = {}
        seq_list = []
        sig_dict = {}
        popped_seqs = {}
        remove_redundant_bool = ut.if_remove_redundant(example)

        # Nested dictionaries are always ordered by the location given in the query and subsequently by the
        # gene_id of the gene given by uniProt.
        for loc_key in data:

            if remove_redundant_bool:
                popped_seqs[loc_key] = ut.del_redundant(data[loc_key], same_seq[loc_key])

            seq_dict[loc_key] = {}
            sig_dict[loc_key] = {}

            for gene_id in data[loc_key]:
                # Fill sequence dict
                seq_dict[loc_key][gene_id] = [data[loc_key][gene_id][-1], data[loc_key][gene_id][0]]

                # Fill signal dict.
                sig_dict[loc_key][gene_id] = data[loc_key][gene_id][-2]

                # Fill sequence list
                seq_list.append(data[loc_key][gene_id][-1])

        # With the sequences as keys the dict contains the number of appearances sorted decreasingly.
        info = {"Warning:": "If redundant have been deleted only the first appearence of each sequence is stored."}
        frequencies = {**info, **dict(sorted(Counter(seq_list).items(), key=operator.itemgetter(1), reverse=True))}

        # Filename input
        if example:
            filename = "example_sa_uniprot.xlsx"
        else:
            name = input("Leave clear for standard 'uniProt_data.xlsx' or enter custom filename: ")
            if name:
                filename = name
            else:
                filename = "uniProt_data.xlsx"

        # Function to print data in a spreadsheet is defined in 'utilities.py'.
        ut.print_sheet(columns, params, data, frequencies, filename, "UniProt")
        return seq_dict, sig_dict, popped_seqs, same_seq, remove_redundant_bool

    # If you want the data in FASTA format
    elif f == "fasta":
        # In this case c is not used.
        fasta_data = search_uniprot(q, columns="", form="fasta")

        # Again every request is stored in a single dict entry by location and each written in a single .fsa
        # FASTA file for later use.
        for loc_key in fasta_data:
            name = input("Leave clear for standard 'uniProt_" + loc_key + ".fsa' or enter custom filename without ending: ")
            if name:
                filename = name + ".fsa"
            else:
                filename = "uniProt_" + loc_key + ".fsa"
            ut.print_file(fasta_data[loc_key], filename)

        return 0

    else:
        print("Error: Unknown file format!\nOnly 'tab' and 'fasta' until now.")
        return -1
