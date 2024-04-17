# helper functions
# added in version 2
import glob
import os
import random
import sys
import time

import numpy as np
import pandas as pd
import logging
import argparse

import unicodedata
from Bio import SeqIO
from matplotlib import colors as mcolors

# variables
fasta_extensions = [".fasta", ".fa", ".fna", ".faa", ".fas"]


def setup_logger(args):
    logger = logging.getLogger(__name__)
    logger.setLevel(args.log_level)

    # Create a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(args.log_level)
    logger.addHandler(console_handler)

    # Create a file handler if log file is provided
    if args.log_file:
        file_handler = logging.FileHandler(args.log_file)
        file_handler.setLevel(args.log_level)
        logger.addHandler(file_handler)

    # Set log format
    log_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(log_format)
    if args.log_file:
        file_handler.setFormatter(log_format)

    return logger


def get_dummy_logger(log_to_console=False):
    # Logger for testing
    logger = logging.getLogger("dummy_logger")
    if log_to_console:
        # Configure logger to output to console
        logger.setLevel(logging.DEBUG)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        logger.addHandler(console_handler)
        # Set log format
        log_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(log_format)
        logger.setLevel(logging.DEBUG)
    else:
        # Set logger level to higher than CRITICAL so nothing gets logged
        logger.setLevel(logging.CRITICAL + 1)
    return logger


def time_track(starting_time, logger, step=str):
    """
    calculate time passed since last time measurement
    """
    now = time.time()
    delta = now - starting_time
    logger.info(f'{step} took {delta} seconds.')
    return now


def interrupt_script(logger):
    user_input = input("Do you want to continue? (y/n): ")

    if user_input.lower() == 'y':
        logger.info("Continuing ...")
        # Continue with the rest of your script here
    elif user_input.lower() == 'n':
        logger.info("Exiting FlexiDot.")
        sys.exit(1)  # Exit the script
    else:
        logger.info("Invalid input. Please enter 'y' to continue or 'n' to exit.")
        interrupt_script(logger)


def check_positive_int(value):
    ivalue = int(value)
    if ivalue < 1:
        raise argparse.ArgumentTypeError(f"{value} is not a positive integer")
    return ivalue


def check_input_files(input_fasta, logger):
    """
    check if files exist, check file extensions, get files as list

    @param logger: logger
    @param inputs: user input string (--input_fasta), might be file(s), dir or auto
    """

    extensions = fasta_extensions

    def file_is_accessible(file_to_check):
        try:
            simple_file_check = open(file_to_check, "r")
            simple_file_check.close()
        except FileNotFoundError as e:
            logger.error(f"{e}")
            sys.exit(1)
        except PermissionError as e:
            logger.error(f"{e}")
            sys.exit(1)

        return True

    input_fasta_checked = []

    logger.debug(f"Checking file inputs: {input_fasta}")

    # check if input exists and is a valid option
    for file in input_fasta:

        # Option 1: file is a file
        if os.path.exists(file) and os.path.isfile(file):
            logger.debug(f'Is file: {file}')

            file_extension = os.path.splitext(file)[1]

            if file_extension in extensions:

                logger.debug(f'From {file}: Is {file} accessible?')
                if file_is_accessible(file):
                    input_fasta_checked.append(file)

            else:
                e = f'{file}: incorrect format. (Missing supported file extension)\n' \
                    f'Supported file extensions: {extensions}'
                logger.error(e)
                raise ValueError(e)

        # Option 2: file is a directory
        elif os.path.isdir(file):
            logger.debug(f'Is directory: {file}')

            cwd = file  # input is not a file but a directory

            logger.info(f"Searching {cwd} for fasta files.")

            matching_files = [get_file for ext in extensions for get_file in glob.glob(os.path.join(cwd, '*' + ext))]

            logger.info(f"Found: {', '.join(matching_files)}.")

            for found_file in matching_files:
                logger.debug(f'From {file}: Is {found_file} accessible?')
                if file_is_accessible(found_file):
                    input_fasta_checked.append(found_file)

        # Option 3: file is auto
        elif file.lower() == 'auto':
            logger.debug(f'Input is auto.')

            cwd = os.getcwd()

            logger.info(f"Searching {cwd} for fasta files.")

            matching_files = [get_file for ext in extensions for get_file in glob.glob(os.path.join(cwd, '*' + ext))]

            logger.info(f"Found: {', '.join(matching_files)}.")

            for found_file in matching_files:
                logger.debug(f'From {file}: Is {found_file} accessible?')
                if file_is_accessible(found_file):
                    input_fasta_checked.append(found_file)

        else:
            logger.error(f'Invalid input {file}.')
            raise FileNotFoundError(f'Invalid input {file}.')

    logger.info(f'Files to be used for plotting:\n'
                    f'{", ".join(input_fasta_checked)}')

    return input_fasta_checked


def check_wobbles(args, defaults, logger):
    # Check input variables: wobbles
    logger.info('Checking wobble conversion.')

    if args.wobble_conversion and defaults['max_N_percentage'] > 49:
        defaults['max_N_percentage'] = 49
        if args.type_nuc:
            defaults['ambiq_res'] = "N"
        else:
            defaults['ambiq_res'] = "X"

    return args, defaults


def csv_to_dict(csv_data, mode):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_data, header=1, sep='\t')

    # Convert the DataFrame to a dictionary

    # default values
    if mode == 'values':
        data_dict = df.set_index("parameter")["value"].to_dict()

    if mode == 'colors':
        data_dict = df.set_index("feature")[["color", "transparency", "zoom"]].apply(tuple, axis=1).to_dict()

    return data_dict


def gff_to_dict(gff_data):
    # Read the GFF3 file into a pandas DataFrame
    df = pd.read_csv(gff_data, sep='\t', comment='#', header=None,
                     names=["seqid", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"])

    # Process the DataFrame to create the desired dictionary
    data_dict = {}

    for index, row in df.iterrows():
        seqid = row["seqid"]
        feature = row["feature"]

        if row["strand"] == '-':
            feature += "_rev"

        start = row["start"]
        end = row["end"]

        if seqid not in data_dict:
            data_dict[seqid] = []

        data_dict[seqid].append((feature, start, end))

    return data_dict


def _degap_fasta(input_fasta):
    """
    remove gaps from fasta - new degapped sequence file created
    """

    # degap all sequence files
    output_fastas = []
    if type(input_fasta) != list:
        input_fasta = list(input_fasta)

    for input_fas in input_fasta:
        base_name, ext = os.path.splitext(input_fas)
        output_fas = base_name + "_degapped" + ext

        with open(input_fas, 'r') as in_file:

            with open(output_fas, 'w') as out_file:

                for line in in_file:
                    if line.startswith(">"):
                        out_file.write(line.strip() + os.linesep)
                    else:
                        out_file.write(line.strip().replace("-", "") + os.linesep)

        output_fastas.append(output_fas)

    return output_fastas


def read_seq(input_fasta, logger):
    """
    read fasta sequences from (all) file(s)

    @param input_fasta: list of fasta file(s), from check_input_files
    """

    # Name for the SQLite index database file
    logger.debug(f'Indexing fasta files... {input_fasta}')
    db_name = ":memory:"  # create the SQLite database in memory

    try:
        # Indexes the sequences from all input FASTA files without concatenating them
        # This creates or updates the .idx file as necessary
        seq_dict = SeqIO.index_db(db_name, input_fasta, "fasta")

        logger.debug(f'Done indexing. Indexed {len(seq_dict)} sequences.')

    except ValueError as e:
        logger.error(f"Error reading fasta sequences - please check input files, e.g. for duplicate names! {e}")
        return {}, []
    except Exception as e:  # It's a good practice to capture specific exceptions
        logger.error(f"Error reading fasta sequences - please check input files! {e}")
        return {}, []

    # NOTE: gaps are now handled similar to Ns and produce gaps in the lines
    # for seq in seq_dict:
    #    if "-" in seq_dict[seq].seq:
    #        # ungapped = seq_dict[seq].seq.ungap("-") # cannot be assigned back to sequence record
    #        logger.info(f'Sequences are degapped prior analysis!!!')
    #        return read_seq(_degap_fasta(input_fasta), logger)

    # get ordered sequence names
    sequences = []
    for seq in seq_dict:
        sequences.append(seq_dict[seq].id)

    return seq_dict, sequences


def read_annotation_colors(in_gff, in_gff_config, logger):
    if in_gff_config is not None:
        logger.info(f'Reading GFF color configuration file: {in_gff_config.name}')
        gff_feat_colors = csv_to_dict(in_gff_config, mode='colors')

    else:
        logger.info(f'Loading default GFF color configuration.')
        current_script_dir = os.path.dirname(os.path.abspath(__file__))
        csv_file_path = os.path.join(current_script_dir, "..", "configs", "config_gff_colors.csv")
        gff_feat_colors = csv_to_dict(csv_file_path, mode='colors')

    logger.info(f'Validating gff_color_config values.')
    for key, values in gff_feat_colors.items():
        color_value = values[0]

        # Validate color
        if not mcolors.is_color_like(color_value):
            # Replace invalid color with 'grey'
            gff_feat_colors[key] = ('grey',) + values[1:]
            logger.info(f'Invalid color. Defaulted {key} to color grey.')

        # Validate alpha value
        try:
            alpha_value = float(values[1])
        except ValueError:
            gff_feat_colors[key] = (values[0], 0.75, values[2])
            logger.info(f'Invalid alpha. Defaulted {key} to transparency of 75%.')

        # Validate transparency value
        try:
            zoom_value = float(values[2])
        except ValueError:
            gff_feat_colors[key] = values[:2] + (0,)
            logger.info(f'Invalid zoom. Defaulted {key} to no zoom.')

    # default coloring for unknown annotations
    if "others" not in gff_feat_colors.keys():
        gff_feat_colors["others"] = ("grey", 0.5, 0)

    logger.debug(f'{gff_feat_colors}')

    return gff_feat_colors


def read_gffs(input_gff_files, color_dict, logger):
    """"
    create feature dictionary from input_gff
    sequence name as key and (feature type, start, stop) as value
    """
    logger.debug(f'Read GFFs input: \n{input_gff_files}\n'
                 f'Read GFFs colors: \n{color_dict}')

    gff_filenames = []

    # GFF files are TextIOWrapper objects from reading with argparse
    for item in input_gff_files:
        print(item.name)
        gff_filenames.append(item.name)

    unknown_feats = set([])
    used_feats = set([])
    feat_dict = {}
    for input_gff in gff_filenames:
        logger.info(f'Reading GFF file: {input_gff}')

        new_entries = gff_to_dict(input_gff)

        feat_dict.update(new_entries)

    # Check if for all features in feat_dict colors are defined in color_dict
    for seqid, features in feat_dict.items():  # Iterate through the items in feat_dict

        logger.debug(f'Checking colors for {seqid} and {features}.')

        for feat_type, _, _ in features:  # get feature type, ignore position info

            if feat_type not in color_dict.keys():  # compare

                logger.info(f'Missing annotation specification for {feat_type}: '
                            f'Defaulting {feat_type} to pink, 0.5, 0')
                # TODO: Add default color, transparency and zoom to config_values.csv
                color_dict[feat_type] = ('purple', 0.5, 0)  # if missing, add default values
                unknown_feats.add(feat_type)

            used_feats.add(feat_type)

    if len(unknown_feats) != 0:
        logger.info(f'Missing shading information for {len(unknown_feats)} feature types. '
                    f'Using default values.')

    logger.debug(f'Annotation colors: {color_dict}\n'
                 f'Used features: {used_feats}')

    return feat_dict, color_dict, used_feats


def alphabets(type_nuc):
    """
    provide ambiguity code for sequences
    """

    nucleotide_alphabet = ["A", "C", "G", "T"]

    nucleotide_alphabet_full = ["A", "C", "G", "T", "N", "B", "D", "H",
                                "V", "Y", "R", "W", "S", "K", "M"]

    nucleotide_ambiguity_code = {"N": ["A", "C", "G", "T"],  # any
                                 "B": ["C", "G", "T"],  # not A
                                 "D": ["A", "G", "T"],  # not C
                                 "H": ["A", "C", "T"],  # not G
                                 "V": ["A", "C", "G"],  # not T
                                 "Y": ["C", "T"],  # pyrimidine
                                 "R": ["A", "G"],  # purine
                                 "W": ["A", "T"],  # weak
                                 "S": ["C", "G"],  # strong
                                 "K": ["G", "T"],  # keto
                                 "M": ["A", "C"]}  # amino

    nucleotide_match_dict = {"N": "[ACGTNBDHVYRWSKM]",  # any
                             "B": "[CGTNBDHVYRWSKM]",  # not A
                             "D": "[AGTNBDHVYRWSKM]",  # not C
                             "H": "[ACTNBDHVYRWSKM]",  # not G
                             "V": "[ACGNBDHVYRWSKM]",  # not T
                             "K": "[GTNBDHVYRWSK]",  # keto - not A,C,M
                             "M": "[ACNBDHVYRWSM]",  # amino - not G,T,K
                             "W": "[ATNBDHVYRWKM]",  # weak - not C,G,S
                             "S": "[CGNBDHVYRSKM]",  # strong - not A,G,W
                             "Y": "[CTNBDHVYWSKM]",  # pyrimidine - not A,G,R
                             "R": "[AGNBDHVRWSKM]",  # purine - not C,T,Y
                             "A": "[ANDHVRWM]",
                             "C": "[CNBHVYSM]",
                             "G": "[GNBDVRSK]",
                             "T": "[TNBDHYWK]"}

    aminoacid_alphabet = ["A", "R", "N", "D", "C", "E", "Q", "G",
                          "H", "I", "L", "K", "M", "F", "P", "S",
                          "T", "W", "Y", "V", "U", "O", "*"]

    aminoacid_alphabet_full = ["A", "R", "N", "D", "C", "E", "Q", "G",
                               "H", "I", "L", "K", "M", "F", "P", "S",
                               "T", "W", "Y", "V", "U", "O", "*", "J",
                               "Z", "B", "X"]

    aminoacid_ambiguity_code = {"J": ["I", "L"],
                                "Z": ["Q", "E"],
                                "B": ["N", "D"],
                                "X": ["A", "R", "N", "D", "C", "E", "Q", "G",
                                      "H", "I", "L", "K", "M", "F", "P", "S",
                                      "T", "W", "Y", "V", "U", "O", "*"]}  # any

    aminoacid_match_dict = {"J": "[ILJ]",
                            "Z": "[QEZ]",
                            "B": "[NDB]",
                            # "X": ".",
                            "X": "[ARNDCEQGHILKMFPSTWYVUO*XBZJ]",
                            "A": "[AX]",
                            "R": "[RX]",
                            "N": "[NXB]",
                            "D": "[DXB]",
                            "C": "[CX]",
                            "E": "[EXZ]",
                            "Q": "[QXZ]",
                            "G": "[GX]",
                            "H": "[HX]",
                            "I": "[IXJ]",
                            "L": "[LXJ]",
                            "K": "[KX]",
                            "M": "[MX]",
                            "F": "[FX]",
                            "P": "[PX]",
                            "S": "[SX]",
                            "T": "[TX]",
                            "W": "[WX]",
                            "Y": "[YX]",
                            "V": "[VX]",
                            "U": "[UX]",
                            "O": "[OX]",
                            "*": "[*X]"}

    aa_only = {'E', 'F', 'I', 'J', 'L', 'O', 'Q', 'P', 'U', 'X', 'Z', '*'}

    if type_nuc == 'nucl':
        return nucleotide_alphabet, nucleotide_alphabet_full, nucleotide_ambiguity_code, nucleotide_match_dict
    else:
        return aminoacid_alphabet, aminoacid_alphabet_full, aminoacid_ambiguity_code, aminoacid_match_dict


def split_diagonals(data, stepsize=1):
    """
    split array if point difference exceeds stepsize
    data = sorted list of numbers
    """
    return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)


def process_diagonals(diag_dict, factor=0):
    # Convert diagonals from dict to separate numpy arrays for x and y
    x, y = [], []
    for diag, points in diag_dict.items():
        x_values = np.sort(list(points))  # Convert to list if not already, then sort
        x.extend(split_diagonals(x_values))
        if factor == 0:
            y.extend(split_diagonals(x_values - diag))  # regular y values
        else:
            y.extend(split_diagonals((factor + diag) - x_values, -1))  # reverse y values
    return np.array(x, dtype=object), np.array(y, dtype=object)


def lcs_from_x_values(x_values):
    """
    calculate length of longest common substring based on nested list of numbers
    """
    if len(x_values) == 0:
        return 0
    # get lengths of each subarray data
    lengths = np.array([len(i) for i in x_values])
    return max(lengths)


def shorten_name(seq_name, max_len=20, title_clip_pos="B"):  # , delim="_"):
    """
    shorten sequence names (for diagram titles)
    """

    # check title length
    if len(seq_name) <= max_len:
        return seq_name

    # take last characters
    if title_clip_pos == "E":
        name = seq_name[len(seq_name) - max_len:]

    # take first characters
    else:
        name = seq_name[:max_len]

    return name


def unicode_name(name):
    """
    replace non-ascii characters in string (e.g. for use in matplotlib)
    """
    unicode_string = eval(f'u"{name}"')
    return unicodedata.normalize('NFKD', unicode_string)  # .encode('ascii','ignore')


def wobble_replacement(sequence, general_ambiguity_code, logger):
    """
    get all degenerated sequences for sequence with ambiguous residues
    (only residues considered that are keys in wobble_dictionary)
    """

    # get positions of ambiguous residues
    wobble_pos = []
    for idx in range(len(sequence)):
        letter = sequence[idx]
        if letter in general_ambiguity_code.keys():
            wobble_pos.append(idx)

    logger.debug(f'Wobble replacement ({sequence}): {len(wobble_pos)} wobbles.')

    # replace one wobble through each iteration by all possible residues
    # repeat if still wobbles in new kmers
    kmer_variants = [sequence]
    while True:
        logger.debug(f'Wobble replacement ({sequence}): {len(kmer_variants)} kmer variants.')
        temp_kmers = set([])
        for kmer in kmer_variants:
            for idx in wobble_pos:
                letter = kmer[idx]
                if letter in general_ambiguity_code.keys():
                    for base in general_ambiguity_code[kmer[idx]]:
                        new_kmer = kmer[:idx] + base + kmer[idx + 1:]
                        temp_kmers.add(new_kmer)
        wobble = False
        for kmer in temp_kmers:
            for idx in range(len(kmer)):
                letter = kmer[idx]
                if letter in general_ambiguity_code.keys():
                    wobble = True
                    break
            if wobble:
                break
        kmer_variants = set(list(temp_kmers)[:])
        if not wobble:
            break

    return kmer_variants


def longest_common_substring(s1, s2):
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return longest


def start_message(count):
    # Message pools
    messages_for_few_args = [
        "Winging it today!",
        "Keeping it simple, huh?",
        "Quick and easy does it!",
    ]

    messages_for_many_args = [
        "Being extra specific today!",
        "All set with detailed instructions!",
        "Ready for a deep dive with all those parameters!",
    ]

    message = f"Starting flexidot..\n (っ◕‿◕)っ <("

    if count <= 5:
        thingy = random.choice(messages_for_few_args)
    else:
        thingy = random.choice(messages_for_many_args)

    return message + thingy + ")"


def determine_layout(args, sequences, logger):
    """
    Determine the layout for the collage based on input arguments and the number of sequences.

    Parameters:
    - args: Namespace containing command-line arguments, including 'collage_output' and 'plot_size'.
    - sequences: List of sequences to be plotted.
    - logger: Logger instance for logging messages.

    Returns:
    - Tuple of (mcols, nrows, multi) where:
      - mcols: Number of columns in the collage.
      - nrows: Number of rows in the collage.
      - multi: Boolean indicating whether a multi-collage layout is used.
    """
    mcols, nrows = args.collage_output
    multi = mcols > 1 or nrows > 1

    if len(sequences) == 0:
        logger.error('No sequences provided for selfdotplot! Terminating!')
        sys.exit(1)
    elif len(sequences) == 1 and multi:
        logger.info(f'Ignoring --collage_output {mcols} {nrows}. Change to default.')
        mcols, nrows, multi = 1, 1, False
    if multi and mcols > len(sequences):
        mcols = len(sequences)
        nrows = 1
        logger.info(f'Selfdotplot Collage: Few sequences - correcting number of rows and columns: cols={mcols}, rows={nrows}')
    elif multi and mcols * (nrows - 1) > len(sequences):
        nrows = ((len(sequences) - 1) // mcols) + 1
        logger.info(f'Selfdotplot Collage: Few sequences - correcting number of rows and columns: cols={mcols}, rows={nrows}')
    if multi and not (nrows == 1 and mcols == 1) and args.plot_size <= args.label_size // 2:
        args.label_size = args.plot_size * 3 // 2
        logger.info(f'Reducing label size for better visualization to {args.label_size}')

    return mcols, nrows, multi, args.label_size


def set_prefix(prefix):
    # prepare prefix and prefix for legend files, if required

    if prefix is not None:
        if prefix[-1] not in ["-", "_"]:
            prefix = prefix + "_"
    else:
        prefix = ""

    return prefix
