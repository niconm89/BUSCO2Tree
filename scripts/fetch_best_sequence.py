#!/usr/bin/env python3
# coding: utf-8
"""
Extract the best scoring BUSCO for the provided busco id
"""
import os

def fetch(busco_id, mode, folder, run_name, gene_set=None):
    """
    :param busco_id: the busco id
    :param mode: the mode, tran or prot
    :param folder: the run folder
    :param run_name: the run name
    :param gene_set: the path to the gene set fasta file
    :return: the sequence
    """

    if folder[-1] != '/':
        folder += '/'

    # load the full table
    full_table = {}
    full_table_file = open('%sfull_table_%s.tsv' % (folder, run_name), 'r')
    for line in full_table_file:
        if not line.startswith('#'):
            if line.strip().split()[0] not in full_table:
                full_table.update({line.strip().split()[0]: [line.strip().split()[1:]]})
            else:
                full_table[line.strip().split()[0]].append(line.strip().split()[1:])

    if mode == 'prot':
        return _fetch_prot(busco_id, folder, full_table, gene_set)
    elif mode == 'tran':
        return _fetch_tran(busco_id, folder, full_table)
    else:
        print('Wrong mode specified')


def _fetch_best(busco_id, folder, full_table):
    """
    :return: the sequence id to retrieve, with the best score
    """
    result_scores = []
    good_results = []
    for entry in full_table[busco_id]:
        result_scores.append(float(entry[2]))
    # open all hmm result file for this busco id
    for file in os.listdir('%shmmer_output/' % folder):
        if file.startswith(busco_id):
            for line in open('%shmmer_output/%s' % (folder, file), 'r'):
                if not line.startswith('#'):
                    if float(line.split()[7]) in result_scores:
                        good_results.append(line)

    # Recheck that the best result has the longest protein, should be the best score as well. 
    # Edit: in fact no, for transcriptomes it is not since the 6 frames translation are evaluated and only one is correct, so the warning below are not really useful, don't focus too much on it.
    id_to_return = None
    best_length = 0
    best_score = 0
    result_with_best_length = None
    result_with_best_score = None

    for result in good_results:
        if float(result.split()[7]) > best_score:
            id_to_return = result.split()[0]
            best_score = float(result.split()[7])
            result_with_best_score = result.split()[0]
        if int(result.split()[16]) - int(result.split()[15]) > best_length:
            best_length = int(result.split()[16]) - int(result.split()[15])
            result_with_best_length = result.split()[0]

    #if result_with_best_length != result_with_best_score:
    #    _logger.warning('best score and best length are not the same for %s in %s' % (busco_id, folder))
    #else:
    #    _logger.debug('best score and best length are the same for %s in %s' % (busco_id, folder))

    return id_to_return


def _fetch_tran(busco_id, folder, full_table):
    """
    :param busco_id:
    :param folder:
    :param full_table:
    :return: the best sequence for a BUSCO id, extracted from a transcriptome BUSCO run
    """

    id_to_return = _fetch_best(busco_id, folder, full_table)

    # Now return the correct sequence
    sequences = open('%stranslated_proteins/%s.faa' % (folder, '_'.join(id_to_return.split('_')[:-1])), 'r')
    sequence = ''
    found = False
    for line in sequences:
        if found and line.startswith('>'):
            break
        elif line.strip() == '>%s' % id_to_return:
            found = True
            sequence += '%s\n' % line.strip()
        elif found:
            sequence += '%s' % line.strip()
    return sequence


def _fetch_prot(busco_id, folder, full_table, gene_set):
    """
    :param busco_id:
    :param folder:
    :param full_table:
    :gene_set:
    :return: the best sequence for a BUSCO id, extracted from a protein BUSCO run     """
    id_to_return = _fetch_best(busco_id, folder, full_table)
    # Now return the correct sequence
    sequences = open('%s' % gene_set, 'r')
    sequence = ''
    found = False
    for line in sequences:
        if found and line.startswith('>'):
            break
        elif line.strip().split(' ')[0] == '>%s' % id_to_return:
            found = True
            sequence += '%s\n' % line.strip()
        elif found:
            sequence += '%s' % line.strip()
    return sequence

