#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
FlexiDot Version 2.02

FlexiDot: Highly customizable ambiguity-aware dotplots for visual sequence investigation

Kathrin M. Seibt, Thomas Schmidt and Tony Heitkam
Institute of Botany, TU Dresden, Dresden, 01277, Germany

Bioinformatics (2018) Vol. 34 (20), 3575–3577, doi 10.1093/bioinformatics/bty395
"""
import os
import sys
import time
import argparse
import yaml
import atexit
from importlib import resources
from scripts import flexiutil, flexiclass, flexidraw, fleximatch


# Load config file
def load_config():
    # Access the package resource using files() and navigate to the config file
    config_path = resources.files("configs") / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def main(args, parameters):
    # ____________Start & log____________

    # Initialize logging
    logger = flexiutil.setup_logger(args=args)

    # Log the exact command line used to start the script
    command_line = ' '.join(sys.argv)
    start_count = sum(1 for arg in sys.argv if arg.startswith('-'))
    logger.info(f"Command invoked: {command_line}")

    logger.info(flexiutil.start_message(start_count))

    # Log the script start and the command line arguments
    args_str = ' '.join(f'{k}={v}\n' for k, v in vars(args).items())
    logger.debug(f"Started with arguments: \n {args_str}")

    # Parse default values from config.yaml file (not user input arguments)
    def_calc_params = parameters['calculation_parameters']

    # Log the parameters
    params_str = ' '.join(f'{k}={v}\n' for k, v in def_calc_params.items())
    logger.debug(f"Started with parameters: \n {params_str}")

    # ____________Check inputs____________

    # Check inputs by user and defaults, update them
    args, def_calc_params = flexiutil.check_wobbles(args=args, defaults=def_calc_params, logger=logger)
    input_fasta_files = flexiutil.check_input_files(input_fasta=args.input_fasta, logger=logger)

    if args.input_gff_files is not None:
        gff_colors = flexiutil.read_annotation_colors(in_gff=args.input_gff_files,
                                                      in_gff_config=args.gff_color_config_file,
                                                      logger=logger)
    else:
        if args.gff_color_config_file is not None:
            logger.info(f'No GFF file is provided but you provided a configuration. '
                        f'Did you forget to provide a GFF file?')
            flexiutil.interrupt_script(logger)
        gff_colors = {}

    # if color is set to white, reverse complementary matches are skipped
    if not args.reverse_complement_search:
        args.line_col_rev = "white"
    elif args.type == "prot":
        logger.info(f'Reverse complement deactivated for proteins!')
        args.line_col_rev = "white"

    # ____________Check and convert mode____________

    # Set modes to match flexidot v1 (old numbers, new strings)
    if isinstance(args.plotting_mode, str):
        plotting_mode = [args.plotting_mode]
    else:
        plotting_mode = args.plotting_mode

    mode_text = []
    modes = []

    # Interchangeable mode identifier
    modes_mapping = {"0": ("self", 0),
                     "self": ("self", 0),
                     "1": ("paired", 1),
                     "paired": ("paired", 1),
                     "2": ("poly", 2),
                     "poly": ("poly", 2)}

    for item in plotting_mode:
        if item in modes_mapping:
            mode_name, mode_value = modes_mapping[item]
            modes.append(mode_value)
            mode_text.append(mode_name)

    modes = list(set(modes))
    mode_text = list(set(mode_text))
    logger.info(f'\n{60 * "="}\n'
                f'{8 * " "}Running plotting modes: {", ".join(mode_text)}'
                f'\n{60 * "="}\n')

    ###################
    # create dotplots #
    ###################

    # ____________create legend____________

    # read gff annotation data if provided for shading
    if args.input_gff_files is not None and len(args.input_gff_files) != 0:
        logger.info(f'Creating figure legend')
        legend_prefix = flexiutil.set_prefix(args.output_file_prefix)
        logger.info(f'Reading {len(args.input_gff_files)} GFF annotation files.')
        feat_dict, color_dict, used_feats = flexiutil.read_gffs(args.input_gff_files,
                                                                color_dict=gff_colors,
                                                                logger=logger)

        # create color legend
        colors, alphas = [], []
        for item in sorted(used_feats):
            colors.append(color_dict[item][0])
            alphas.append(color_dict[item][1])

        flexidraw.legend_figure(colors=colors, lcs_shading_num=len(used_feats),
                                # TODO: fix lcs shading, move legend creation to end
                                bins=sorted(used_feats),
                                alphas=alphas, gff_legend=True,
                                prefix=legend_prefix, filetype=args.filetype,
                                logger=logger)

        logger.info(f'GFF feature types: {", ".join(sorted(used_feats))}\n'
                    f'Colors: {", ".join(sorted(colors))}\n')

    # ____________read sequences____________
    # create a sequence dict in memory using SeqIO.index
    # get sequence ids as list
    seq_dict, sequences = flexiutil.read_seq(input_fasta_files, logger)

    if seq_dict == {}:
        logger.error("Failed to load sequences")
        return list()

    if args.alphabetic_sorting:
        sequences = sorted(sequences)

    # ____________self dotplots____________

    t1 = time.time()

    if 0 in modes:

        # TODO: match sequences and return seq.id + lines

        list_of_fig_names = flexidraw.selfdotplot(fasta_list=input_fasta_files,
                                                  args=args,
                                                  defaults=def_calc_params,
                                                  gff_colors=gff_colors,
                                                  logger=logger)

        if list_of_fig_names is not None and len(list_of_fig_names) != 0:
            text = "-> Image file(s): %s\n" % (", ".join(list_of_fig_names))
        else:
            text = "No image files were created!\n"
        flexiutil.time_track(t1, logger, step='Creating all plots')
        logger.info(text)
    '''
    # paired dotplots
    if 1 in modes:

        # check if multiple dotplots in collage output
        mcols, nrows = args.collage_output

        if mcols > 1 or nrows > 1:
            multi = True
        else:
            multi = False

        if multi:
            list_of_fig_names = pairdotplot(seq_list, wordsize, prefix=prefix, label_size=label_size, title_length=title_length, title_clip_pos=title_clip_pos, plot_size=plot_size, filetype=filetype, type_nuc=type_nuc, convert_wobbles=convert_wobbles, rc_option=rc_option, substitution_count=substitution_count, alphabetic_sorting=alphabetic_sorting, only_vs_first_seq=only_vs_first_seq, multi=multi, ncols=ncols, nrows=nrows, length_scaling=length_scaling, mirror_y_axis=mirror_y_axis, max_N_percentage=max_N_percentage, verbose=verbose)
            t1 = time_track(t1)
        else:
            if not length_scaling:
                    text = "\nPairwise dotplot with individual output files scaled by sequence length automatically!"
                    logprint(text, start=False, printing=True)
            list_of_fig_names = pairdotplot(seq_list, wordsize, prefix=prefix, label_size=label_size, title_length=title_length, title_clip_pos=title_clip_pos, plot_size=plot_size, filetype=filetype, type_nuc=type_nuc, convert_wobbles=convert_wobbles, rc_option=rc_option, substitution_count=substitution_count, alphabetic_sorting=alphabetic_sorting, only_vs_first_seq=only_vs_first_seq, multi=multi, ncols=ncols, nrows=nrows, length_scaling=True, mirror_y_axis=mirror_y_axis, max_N_percentage=max_N_percentage, verbose=verbose)
            t1 = time_track(t1)
        if list_of_fig_names != None and len(list_of_fig_names) != 0:
            text = "-> Image file(s): %s\n" % (", ".join(list_of_fig_names))
        else:
            text = "No image files were created!\n"
        logprint(text, start=False, printing=True)
        logprint(60*"=")
    '''
    '''

    # all-against-all dotplot
    if 2 in modes:
        list_of_fig_names = polydotplot(seq_list, wordsize, prefix=prefix, label_size=label_size, title_length=title_length, title_clip_pos=title_clip_pos, plot_size=plot_size, filetype=filetype, type_nuc=type_nuc, convert_wobbles=convert_wobbles, rc_option=rc_option, substitution_count=substitution_count, alphabetic_sorting=alphabetic_sorting, lcs_shading=lcs_shading, lcs_shading_num=lcs_shading_num, lcs_shading_ref=lcs_shading_ref, lcs_shading_interval_len=lcs_shading_interval_len, lcs_shading_ori=lcs_shading_ori, input_user_matrix_file=input_user_matrix_file, user_matrix_print=user_matrix_print, spacing=spacing, gff_files=gff, gff_color_dict=gff_feat_colors, mirror_y_axis=mirror_y_axis, max_N_percentage=max_N_percentage, verbose=verbose)
        t1 = time_track(t1)
        if list_of_fig_names != None and len(list_of_fig_names) != 0:
            text = "-> Image file(s): %s\n" % (", ".join(list_of_fig_names))
        else:
            text = "No image files were created!\n"
        logprint(text, start=False, printing=True)
        logprint(60*"=")

    '''
    logger.debug(f'\nDefaults:\n{def_calc_params}\nParameters:\n{args}')

    logger.info(f'\n{60 * "="}\n'
                f'{8 * " "}Thank you for using FlexiDot!'
                f'\n{60 * "="}\n')

    ############
    # clean up #
    ############

    #os.remove(helpers.fasta_db_idx)

    sys.exit(0)


if __name__ == "__main__":
    # Load config file
    config = load_config()
    # split out values from config
    default_args = config['cmd_input_arguments']
    default_params = config['misc_parameters']

    # Get user input arguments
    parser = argparse.ArgumentParser(
        description=
        '''
        FlexiDot: Highly customizable ambiguity-aware dotplots for visual sequence investigation
        DEPENDENCIES:
            - python 3.8 or higher with packages
            ''',
        epilog="""""",
        formatter_class=flexiclass.CustomFormatter)

    # Parser groups for user-friendly help, sort arguments into custom groups
    requiredNamed = parser.add_argument_group('Required arguments:')
    outParam = parser.add_argument_group('Output options:')
    calcParam = parser.add_argument_group('Calculation parameters:')
    graphFormat = parser.add_argument_group('Graphic formatting:')
    gffShading = parser.add_argument_group('GFF shading (only plotting modes self, poly):')
    lcsShading = parser.add_argument_group('LCS shading (only plotting mode poly):')
    usrMatrixShading = parser.add_argument_group('User matrix shading (only plotting mode poly):')
    miscParameter = parser.add_argument_group('Misc parameter:')
    logLogging = parser.add_argument_group('Logging:')

    # Input arguments (required fasta file(s))
    requiredNamed.add_argument(
        "-i", "--input_fasta",
        default=None,
        required=True,
        type=str,
        nargs='+',
        help='''Input fasta file(s) (Use multiple files separated by spaces).
        Use "-i/--input auto" to use all files in current directory.
        Use "-i/--input <path_to_directory>" to use all files in specified directory.
        Recognized file-extensions are: ".fasta", ".fa", ".fna", ".faa", ".fas"
        Required!''')

    # Output options
    outParam.add_argument(
        "-o", "--output_file_prefix",
        type=str,
        default=default_args['output_options']['output_file_prefix'],
        help="File prefix to be added to the generated filenames.")
    outParam.add_argument(
        "-f", "--filetype",
        type=str,
        default=default_args['output_options']['filetype'],
        choices=['PNG', 'PDF', "SVG"],
        help="Output file format.")
    outParam.add_argument(
        "-c", "--collage_output",
        nargs=2,
        default=default_args['output_options']['collage_output'],
        type=flexiutil.check_positive_int,
        help='''Multiple dotplots are combined in a collage. 
        Provide two positive integers for number of columns and rows separated by a space.
        Old flexidot default can be recreated using -c 4 5.''')
    outParam.add_argument(
        "-s", "--alphabetic_sorting",
        action='store_true',
        default=default_args['output_options']['alphabetic_sorting'],
        help="Sort sequences alphabetically according to titles.")

    # Calculation parameters
    calcParam.add_argument(
        "-k", "--wordsize",
        type=int,
        default=default_args['calculation_parameters']['wordsize'],
        help="Word size (kmer length) for dotplot comparison.")
    calcParam.add_argument(
        "-p", "--plotting_mode",
        type=str,
        default=default_args['calculation_parameters']['plotting_mode'],
        choices=['0', 'self', '1', 'paired', '2', 'poly'],
        nargs="+",
        help='''Mode of FlexiDot dotplotting (Separate by space to run multiple modes).''')
    calcParam.add_argument(
        "-t", "--type",
        type=str,
        default=default_args['calculation_parameters']['type'],
        choices=['nucl', 'prot'],
        help="Type of residue is nucleotide or amino acid.")
    calcParam.add_argument(
        "-w", "--wobble_conversion",
        action='store_true',
        default=default_args['calculation_parameters']['wobble_conversion'],
        help="Turn on ambiguity handling for relaxed matching. Note: Increases computational time significantly.")
    calcParam.add_argument(
        "-S", "--substitution_count",
        type=int,
        default=default_args['calculation_parameters']['substitution_count'],
        help="Number of substitutions (mismatches) allowed per window for relaxed matching. Note: Increases "
             "computational time significantly.")
    calcParam.add_argument(
        "-r", "--reverse_complement_search",
        action=flexiclass.CustomStoreBooleanAction,
        default=default_args['calculation_parameters']['reverse_complement_search'],
        help="Find reverse complementary matches (only for --type nucl).")
    calcParam.add_argument(
        "-O", "--only_vs_first_seq",
        action='store_true',
        default=default_args['calculation_parameters']['only_vs_first_seq'],
        help="Turn on to limit pairwise comparisons to match all sequences to 1st sequence only "
             "(only if --plotting_mode=1).")

    # narrow_diagonal_interval currently not in use, faster plotting enabled by multi-threading
    '''calcParam.add_argument(
        "-N", "--narrow_diagonal_interval",
        action=custom_classes.StoreIfSet,
        default=0,
        type=int,
        help="Interval size for quicker narrow self-dotplot. "
             "Limits calculation to interval of given size along the main diagonal "
             "(only if --plotting_mode=0).")'''

    # Graphics formatting
    graphFormat.add_argument(
        "-A", "--line_width",
        type=float,
        default=default_args['graphic_formatting']['line_width'],
        help="Line width.")
    graphFormat.add_argument(
        "-B", "--line_col_for",
        type=str,
        default=default_args['graphic_formatting']['line_col_for'],
        help="Line color for forward matches.")
    graphFormat.add_argument(
        "-C", "--line_col_rev",
        type=str,
        default=default_args['graphic_formatting']['line_col_rev'],
        help="Line color for reverse matches.")
    graphFormat.add_argument(
        "-D", "--x_label_pos",
        action='store_true',
        default=default_args['graphic_formatting']['x_label_pos'],
        help="Switch x-lable position (top to bottom). (default = top)")
    graphFormat.add_argument(
        "-E", "--label_size",
        type=int,
        default=default_args['graphic_formatting']['label_size'],
        help="Font size.")
    graphFormat.add_argument(
        "-F", "--spacing",
        type=float,
        default=default_args['graphic_formatting']['spacing'],
        help="Spacing between all-against-all dotplots (only with --plotting_mode poly).")
    graphFormat.add_argument(
        "-M", "--mirror_y_axis",
        action='store_true',
        default=default_args['graphic_formatting']['mirror_y_axis'],
        help="Turn on y-axis flipping (bottom to top).")
    graphFormat.add_argument(
        "-P", "--plot_size",
        type=int,
        default=default_args['graphic_formatting']['plot_size'],
        help="Plot size.")
    graphFormat.add_argument(
        "-R", "--representation",
        default=default_args['graphic_formatting']['representation'],
        choices=[0, 1, 2],
        help='''Region of plot to display (only with --plotting_mode poly):
            0... full
            1... upper
            2... lower\n''')
    graphFormat.add_argument(
        "-L", "--length_scaling",
        action='store_true',
        default=default_args['graphic_formatting']['length_scaling'],
        help='''Turn on axis scaling by sequence length (only with --plotting_mode poly). 
        By default plots are squared.''')
    graphFormat.add_argument(
        "-T", "--title_length",
        nargs="+",
        default=default_args['graphic_formatting']['title_length'],
        action=flexiclass.TitleLengthAction,
        help='''Limit title length (sequence name length in plots).
        Specify the position of selection from beginning (B) or end (E).
        Example: Select the last 10 characters to plot by using -T 10 E.\n''')

    # GFF shading parameters
    gffShading.add_argument(
        "-g", "--input_gff_files",
        type=argparse.FileType('r'),
        nargs='+',
        help="GFF3 file used for markup in self-dotplots (file name or space-separated file list).")
    gffShading.add_argument(
        "-G", "--gff_color_config_file",
        type=argparse.FileType('r'),
        help='''Tab-delimited config file for custom gff shading:
        column 1: feature type
        column 2: color
        column 3: alpha
        column 4: zoom factor (for small regions)''')

    # LCS Shading
    lcsShading.add_argument(
        "-x", "--lcs_shading",
        action='store_true',
        default=default_args['lcs_shading']['lcs_shading'],
        help="Turn on shading of subdotplot based on length of longest common substring (LCS).")
    lcsShading.add_argument(
        "-X", "--lcs_shading_num",
        type=int,
        default=default_args['lcs_shading']['lcs_shading_num'],
        help="Number of shading intervals (hues) for LCS and user matrix shading (-u).")
    lcsShading.add_argument(
        "-y", "--lcs_shading_ref",
        default=default_args['lcs_shading']['lcs_shading_ref'],
        choices=[0, 1, 2],
        help='''Reference for LCS shading:
        0... maximal LCS length
        1... maximally possible length (length of shorter seq in pairwise comparison)
        2... given interval sizes by --lcs_shading_interval_len\n''')
    lcsShading.add_argument(
        "-Y", "--lcs_shading_interval_len",
        action=flexiclass.StoreIfSet,
        type=int,
        default=default_args['lcs_shading']['lcs_shading_interval_len'],
        help="Length of interval for LCS shading (only for --lcs_shading_ref 2). "
             "AA (default: 10); DNA")
    lcsShading.add_argument(
        "-z", "--lcs_shading_ori",
        default=default_args['lcs_shading']['lcs_shading_ori'],
        choices=[0, 1, 2],
        help='''Shade subdotplots according to LCS on:
        0... forward strand
        1... reverse strand
        2... both strands (forward above, reverse below diagonal; for use with -u best LCS is used)\n''')

    # User matrix shading
    usrMatrixShading.add_argument(
        "-u", "--input_user_matrix_file",
        type=argparse.FileType('r'),
        help='''CSV file for shading above diagonal according to values in matrix file:
            column 1: sequence name
            column 2-n: values (e.g., identity matrix from ms-alignment, strings are ignored)''')

    # Misc
    miscParameter.add_argument(
        "-cpu", "--threads",
        type=int,
        default=1,
        help="Set number of cpu cores to use for multi-threaded processing steps."
    )

    # Logging parameters
    logLogging.add_argument(
        "--log_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default=default_args['logging']['log_level'],
        help="Set the logging level.")

    logLogging.add_argument(
        "--log_file",
        default=default_args['logging']['log_file'],
        help="Path to the log file. Existing files will be appended.")

    cmd_args = parser.parse_args()

    main(args=cmd_args, parameters=default_params)
