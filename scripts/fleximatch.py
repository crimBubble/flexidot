import numpy as np
import logging
import regex
from Bio.Seq import Seq

from scripts import flexiutil


def find_match_pos(seq1, seq2, wordsize, type_nuc, logger,
                   rc_option=True, report_lcs=False, convert_wobbles=False,
                   max_N_percentage=49, substitution_count=0):
    """
    Finds all matching positions between two sequences where matches are >= word size, using a regular expression search.
    This function supports fuzzy matching by allowing up to a specified number of substitutions.
    It can also convert matching points into lines representing the length of the match, and handle ambiguities if enabled.

    Parameters:
    - seq1 (str): The first sequence.
    - seq2 (str): The second sequence.
    - wordsize (int): The minimum word size (kmer length) for matching positions.
    - type_nuc (bool): Indicates if sequences are nucleotide (True) or protein (False).
    - logger (logging.Logger): Logger instance for logging messages.
    - rc_option (bool, optional): If True, considers reverse complement matching for nucleotide sequences.
      Defaults to True.
    - report_lcs (bool, optional): If True, calculates and reports the length of the longest common substring (LCS) as
      part of the return. Defaults to False.
    - convert_wobbles (bool, optional): If True, handles wobbles/ambiguities in matching positions. Defaults to False.
    - max_N_percentage (int, optional): The maximum allowed percentage of ambiguous bases (Ns for nucleotides, Xs for
      proteins) in a word. Defaults to 49.
    - substitution_count (int, optional): The number of substitutions (mismatches) allowed per word for fuzzy matching.
      Defaults to 0.

    Returns:
    - A tuple containing numpy arrays for x and y coordinates of matching positions in the forward and
    reverse direction, optionally followed by LCS lengths for forward and reverse matches if report_lcs is True.
    Specifically, returns (x1, y1, x2, y2) or (x1, y1, x2, y2, lcs_for, lcs_rev) based on report_lcs.

    Note:
    - This function processes diagonal lines considering only matches within a certain distance to the middle
    diagonal if 'narrow_diagonal_interval' is specified (currently commented out and reserved for future optimization).
    """

    # read sequences
    seq_one, seq_two = seq1.upper(), (seq2.upper())
    # len_one,
    len_two = len(seq_two)

    # look for Ns in DNA or Xs in proteins (minimum word size)
    any_residue = "N" if type_nuc == 'nucl' else "X"

    # set ambiguity code for wobble replacement
    general_ambiguity_code = flexiutil.alphabets(type_nuc)[2]
    ambiguity_match_dict = flexiutil.alphabets(type_nuc)[3]
    ambig_residues = f"[{''.join(sorted(general_ambiguity_code.keys()))}]"

    # Check for wobble presence in either sequence
    wobble_found = any(residue in seq_one or residue in seq_two for residue in ambig_residues) \
        if convert_wobbles else False

    # See where the wobbles are in a sequence (only debugging)
    if logger.getEffectiveLevel() == logging.DEBUG:

        def find_and_print_wobble_indexes(seq, residues, seq_name):
            for index, char in enumerate(seq):
                if char in residues:
                    logger.debug(f"Found '{char}' in {seq_name} at index {index}")

        # Check both sequences
        find_and_print_wobble_indexes(seq_one, ambig_residues, "seq_one")
        find_and_print_wobble_indexes(seq_two, ambig_residues, "seq_two")

    # dictionary for matches
    diag_dict_for = {}
    diag_dict_rc = {}
    counter = [0, 0]

    # Start with the default "forward" matching case
    data_list = [(str(seq_one), str(seq_two), diag_dict_for, 0, "for")]
    # Append the "reverse complement" matching case if rc_option is enabled
    if rc_option:
        data_list.append((str(seq_one), str(seq_two.reverse_complement()), diag_dict_rc, 1, "rc"))

    # Loop over forward and reverse data_lists
    for (seq_query, seq_target, diag_dict, counter_pos, diag_ori) in data_list:

        # Exit loop if counter position 1 is reached but reverse option is turned off
        if not rc_option and counter_pos == 1:
            break

        # Split query sequence into kmers
        for idx in range(len(str(seq_query)) - wordsize + 1):  # loop over whole query sequence

            kmer = str(seq_query)[idx:idx + wordsize]  # current kmer

            # Skip kmer if to many Ns or gaps
            if (kmer.count(any_residue) * 100. / wordsize >= max_N_percentage or
                    kmer.count("-") * 100. / wordsize >= max_N_percentage):
                continue

            # Update kmer if wobbles allowed and found
            if convert_wobbles and wobble_found:
                kmer_string = ""
                # Replace each residue with matching residues or wobbles
                for jdx in range(len(kmer)):
                    kmer_string += ambiguity_match_dict[kmer[jdx]]
            else:
                kmer_string = kmer

            # Convert to regular expression tolerating substitution errors
            if type(substitution_count) is int and substitution_count != 0:
                kmer_string = f"({kmer_string}){{s<={substitution_count}}}"

            # Search for regular expression in target sequence
            kdx = 0
            # Initial search
            match = regex.search(kmer_string, seq_target[kdx:])

            # Search for regex until no more matches are found or end of target seq is reached
            while match is not None:

                counter[counter_pos] += 1

                # Process the match, positions are relative to kdx
                match_start = match.start()
                match_end = match.end()

                kmer2 = match.group()

                # skip excessive N/X stretches (big black areas)
                if kmer2.count(any_residue) * 100. / wordsize <= max_N_percentage:
                    diag = idx - (kdx + match_start)
                    points = set(range(idx + 1, idx + wordsize + 1))

                    # Add diagonal or elongate existing
                    if diag not in diag_dict:
                        diag_dict[diag] = points
                    else:
                        diag_dict[diag].update(points)

                kdx += match_start + 1

                match = regex.search(kmer_string, seq_target[kdx:])

    if rc_option:
        text = f"[matches: {counter[0]} for; {counter[1]} rc]"
    else:
        text = f"[matches: {counter[0]} for; no rc]"
    logger.debug(text)

    logger.debug(f'Diagonal dict forward ({len(diag_dict_for)} matches)')
    logger.debug(f'Diagonal dict reverse ({len(diag_dict_rc)} matches)')

    # Convert diagonals from dict to separate numpy arrays for x and y
    x1, y1 = flexiutil.process_diagonals(diag_dict_for)
    x2, y2 = flexiutil.process_diagonals(diag_dict_rc, len_two + 1) if rc_option else (
        np.array([], dtype=object), np.array([], dtype=object))

    # Handling LCS part outside and simplifying the return statement
    # TODO: LCS calculation needs to be moved for multithreading and tiling
    lcs_for = flexiutil.lcs_from_x_values(x1) if report_lcs else None
    lcs_rev = flexiutil.lcs_from_x_values(x2) if report_lcs else None

    return_value = (x1, y1, x2, y2, lcs_for, lcs_rev) if report_lcs else (x1, y1, x2, y2)

    print(f"X1: {x1}")
    print(f"Y1: {y1}")
    return x1, y1, x2, y2


# test
test = False

if test:

    seq1 = Seq("CCAAAGTGACAGGATATCGAGTCCGG")
    seq2 = Seq("AATCAGTGACAGGATATCCTGTCACTGATTAGTGA")
    wordsize = 5
    type_nuc = "nucl"
    rc_option = True
    max_N_percentage = 49
    report_lcs = False
    convert_wobbles = False
    substitution_count = 0
    logger = helpers.get_dummy_logger(log_to_console=True)

    find_match_pos(seq1, seq2, wordsize, type_nuc, logger)
