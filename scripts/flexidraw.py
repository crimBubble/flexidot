import sys
import time

# matplotlib
import matplotlib.collections as cllct
import matplotlib.patches as patches
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, Manager
from matplotlib.backends.backend_pdf import PdfPages

from scripts import flexiutil
from scripts.fleximatch import find_match_pos


def legend_figure(logger, colors, lcs_shading_num,
                  unit="%", filetype="png", max_len=None, min_len=0,
                  bins=[], alphas=[], gff_legend=False, prefix=""):
    """
    create figure color legend
    """
    max_legend_length_row = 8
    max_legend_length_col = 4

    # check if length of information fit
    if not gff_legend and ((len(bins) != 0 and len(colors) != lcs_shading_num + 1) or (
            len(bins) != 0 and len(colors) != len(bins) + 1)):
        if len(bins) != 0 and len(colors) != lcs_shading_num + 1:
            text = f"**Attention**\nlcs_shading_num ({lcs_shading_num}) does not match number of colors ({len(bins)})!\n"
        elif len(bins) != 0 and len(colors) != len(bins) + 1:
            text = f"**Attention**\nnumber of LCS length bins ({len(colors)}) does not match number of colors ({len(bins)})!\n"
        logger.info(text)
    elif gff_legend and len(bins) != len(colors):
        text = f"**Attention**\nnumber of GFF Feature Types ({len(colors)}) does not match number of colors ({len(bins)})!\n"
        logger.info(text)

    # set alpha values to opaque if none are provided
    if len(alphas) == 0:
        for item in colors:
            alphas.append(1)

    # legend data points
    data_points = list(range(len(colors)))
    if not gff_legend:

        # specify intervals, if max_len provided
        if max_len is not None:
            multi_factor = 100  # one digit
            if max_len <= 1:
                multi_factor = 1000  # two digits
            # len_interval_size = (max_len-min_len) * multi_factor *1. // lcs_shading_num * (1./ multi_factor)
            len_interval_size = (max_len - min_len) * 1. / lcs_shading_num
            len_pos = [float("%.2f" % min_len)]
            # calculate interval positions
            for idx in range(lcs_shading_num):
                len_pos.append(float("%.2f" % (len_pos[-1] + len_interval_size)))

            if "custom-matrix" in prefix.lower() and (0 <= max_len <= 100 and 0 <= min_len <= 100):
                unit = "%"
            elif "custom-matrix" in prefix.lower():
                unit = ""

            text = f"{lcs_shading_num + 1} Legend intervals from {min_len:.2f} to {max_len:.2f}:\n\t{', '.join(len_pos)} - number: {len(len_pos)}, step: {len_interval_size:.2f}, unit: {unit}\n"
            logger.info(text)
            pos = len_pos
            interval_size = len_interval_size
        # generate legend labels acc. to standard interval notation
        else:
            # use default max_len = 100 and min_len = 0
            len_interval_size = 100. / lcs_shading_num
            pos = [float("%.2f" % 0)]
            # calculate interval positions
            for idx in range(lcs_shading_num):
                pos.append(float("%.2f" % (pos[-1] + len_interval_size)))

            # interval_size = 100 // lcs_shading_num
            # pos = list(range(interval_size, 101+interval_size, interval_size))

        # remove unnecessary zeros in decimal places (i.e. if x.x00 in all entries)
        while True:
            last_digit_all_zero = True
            no_delim = False
            for idx in range(len(pos)):
                # only process if fraction with decimal places
                if not "." in str(pos[idx]):
                    no_delim = True
                    break
                # only process when all entries end in zero
                elif str(pos[idx])[-1] != "0":
                    last_digit_all_zero = False
                    break
            if not last_digit_all_zero or no_delim:
                break
            # remove last decimal place (== 0) from all entries
            else:
                temp_pos = pos[:]
                for idx in range(len(pos)):
                    if not str(pos[idx])[-2] == ".":
                        pos[idx] = float(str(pos[idx])[:-1])
                    else:
                        pos[idx] = int(str(pos[idx])[:-2])
                logger.info(f'Shortening legend entries: {temp_pos} - {pos}')

        # eliminate fractions if unit == bp/aa
        if unit in ["aa", "bp"]:
            for idx in range(len(pos)):
                temp_pos = pos[:]
                rounded_unit = False
                if "." in str(pos[idx]):
                    rounded_unit = True
                    # round values up to next integer (keep integer, if not a fraction)
                    pos[idx] = int(pos[idx] // 1) + int(pos[idx] % 1 > 0)
                    if idx == len(pos) - 1 and pos[idx] == 101:
                        pos[idx] = 100
            if rounded_unit:
                logger.info(f"Fractions not permitted for unit '{unit}': {temp_pos} -> {pos}")

        if len(bins) != 0:  # labels provided
            legend_labels = bins[:]
            legend_labels.append("max")
            legend_labels_lengths = []
            for item in bins:
                legend_labels_lengths.append(f"[{item - min(bins)} {unit}, {item} {unit})")
            if len(bins) == len(colors) - 1:
                legend_labels_lengths.append(f"[{max(bins)} {unit}, \u221E]")  # infinite

        else:
            legend_labels = []
            legend_labels_lengths = []
            for idx in range(len(pos)):
                num = pos[idx]
                try:
                    legend_labels.append("[%d%%, %d%%)" % (num - interval_size, num))
                except:
                    legend_labels.append("[%d%%, %d%%)" % (num - len_interval_size, num))
                if max_len is not None:
                    num = len_pos[idx]
                    # as int or float
                    if num == int(num) and int(len_interval_size) == len_interval_size:
                        legend_labels_lengths.append("[%d %s, %d %s)" % (num, unit, num + len_interval_size, unit))
                    else:
                        legend_labels_lengths.append("[%.2f %s, %.2f %s)" % (num, unit, num + len_interval_size, unit))
            legend_labels[-1] = "100" + unit
            if max_len is not None:
                if num == int(num) and int(len_interval_size) == len_interval_size:
                    legend_labels_lengths[-1] = u"[%d %s, \u221E]" % (max_len, unit)
                else:
                    legend_labels_lengths[-1] = u"[%.2f %s, \u221E]" % (max_len, unit)

    # set labels and choose file name
    if gff_legend:
        label_text = bins[:]
        edge_col = None
        legend_file_name = "Legend_GFF_Shading_n%d.%s" % (lcs_shading_num, filetype)
    elif max_len is not None:
        label_text = legend_labels_lengths[:]
        edge_col = "black"
        legend_file_name = "Legend_LCS_Shading_max%d%s_n%d.%s" % (max_len, unit, lcs_shading_num + 1, filetype)
    elif len(bins) != 0:
        label_text = legend_labels_lengths[:]
        edge_col = "black"
        legend_file_name = "Legend_LCS_Shading_%d%s_n%d.%s" % (bins[0], unit, lcs_shading_num + 1, filetype)
    else:
        label_text = legend_labels[:]
        edge_col = "black"
        legend_file_name = "Legend_LCS_Shading_%%len_n%d.%s" % (lcs_shading_num + 1, filetype)

    if prefix is not None and prefix != "":
        if prefix[-1] not in ["-", "_"]:
            prefix = prefix + "_"
        if "custom-matrix" in prefix.lower():
            prefix = prefix.replace("custom-matrix", "")[1:]
            legend_file_name = prefix + legend_file_name.replace("LCS", "CustomMatrix")
        elif "GFF" in legend_file_name or "LCS" in legend_file_name:
            legend_file_name = prefix + legend_file_name

    # plot legend figure
    fig, ax = plt.subplots(3, 1, figsize=(len(colors) * 2, len(colors) * 2))
    for idx in range(len(colors)):
        ax[0].bar(data_points[idx] + 1, data_points[idx] + 1, color=colors[idx], label=label_text[idx],
                  alpha=alphas[idx], edgecolor=edge_col)
        ax[1].bar(data_points[idx] + 1, 0, color=colors[idx], label=label_text[idx],
                  alpha=alphas[idx], edgecolor=edge_col)
        ax[2].bar(data_points[idx] + 1, 0, color=colors[idx], label=label_text[idx],
                  alpha=alphas[idx], edgecolor=edge_col)
    ax[1].set_ylim(0, 1)
    ax[2].set_ylim(0, 1)
    ax[1].legend(ncol=((len(colors) - 1) // max_legend_length_row) + 1, framealpha=1)  # vertical legend
    col_num = len(colors)
    if len(colors) > max_legend_length_col:
        remainder = 0
        if len(colors) % max_legend_length_col != 0:
            remainder = 1
        row_num = len(colors) // max_legend_length_col + remainder
        remainder = 0
        if len(colors) % row_num != 0:
            remainder = 1
        col_num = len(colors) // row_num + remainder
    ax[2].legend(ncol=col_num, framealpha=1)  # horizontal legend

    plt.savefig(legend_file_name)

    return legend_file_name


def calc_fig_ratio(mcols, nrows, plot_size):
    """
    calculate size ratio for given number of columns (ncols) and rows (nrows)
    with plot_size as maximum width and length
    """
    ratio = mcols * 1. / nrows
    if mcols >= nrows:
        figsize_x = plot_size
        figsize_y = plot_size * 1. / ratio
    else:
        figsize_x = plot_size * ratio
        figsize_y = plot_size
    return figsize_x, figsize_y


def selfdotplot(fasta_list, args, defaults, gff_colors, logger):
    """
    self-against-self dotplot
    partially from biopython cookbook

    Parameters:
    - fasta_list: list of checked fasta files
    - args: user input options
    - defaults: default values
    - gff_colors: Color dict with feature as key, color, transparency and zoom as values
    - logger: logging
    """

    '''# read sequences
    seq_dict, sequences = flexiutil.read_seq(fasta_list, logger)

    if seq_dict == {}:
        logger.error("Failed to load sequences")
        return list()

    if args.alphabetic_sorting:
        sequences = sorted(sequences)'''

    # check if multiple dotplots in collage output
    mcols, nrows, multi, args.label_size = flexiutil.determine_layout(args, sequences, logger)

    # prepare prefix and prefix for legend files, if required
    prefix = flexiutil.set_prefix(args.output_file_prefix)

    if args.type == 'nucl':
        aa_bp_unit = 'bp'
    else:
        aa_bp_unit = 'aa'
        if hasattr(args, 'lcs_shading_interval_len_set') and args.lcs_shading_interval_len:
            logger.debug(f'Setting unit to aa and shading interval length set to {args.lcs_shading_interval_len}.')
        else:
            args.lcs_shading_interval_len = 10

    logger.info(f'\n{60 * "="}\n'
                f'{8 * " "}Creating {len(sequences)} selfdotplot images.'
                f'\n{60 * "="}\n')

    name_graph = "Selfdotplot"

    # Create suffix for output files
    suffix = ''
    if args.substitution_count != 0:
        suffix += f'_S{args.substitution_count}'
    if args.wobble_conversion:
        suffix += "_wobbles"
    if args.reverse_complement_search:
        suffix += "_rc"
    if multi:
        suffix += "_collage"

    # calculate fig ratios
    figsize_x, figsize_y = calc_fig_ratio(mcols, nrows, args.plot_size)

    # Clear any previous plots
    plt.cla()

    if multi:
        page_counter = 1
        fig = plt.figure(figsize=(figsize_x, figsize_y))

    list_of_fig_names = []

    if args.filetype == 'PDF':
        fig_name = f'{prefix}{name_graph}_ws{args.wordsize}{suffix}.pdf'
        list_of_fig_names.append(fig_name)
        pdf = PdfPages(fig_name)

    counter = 0
    # TODO: add multi-threading (start)
    #  thread over seq_names, return dict/pf-df with seq_name, counter, lines & col_lists (save for user?)
    #  pass to matplotlib
    #  Idea: create plot grid layout before with initial diagonal and annotation,
    #  add additional lines and shading afterwards

    for seq_name in sequences:
        t1 = time.time()
        logger.info(f'{seq_name}')

        counter += 1

        if not multi:
            plt.cla()

        # read sequence
        seq_record = seq_dict[seq_name]
        name_seq = seq_record.id
        seq_one = seq_record.seq
        length_seq = len(seq_one)
        logger.debug(seq_one[1:50])

        # get positions of matches
        logger.debug(f'Matching using regex search.')

        x_lists, y_lists, x_lists_rc, y_lists_rc = (
            find_match_pos(seq_one, seq_one, args.wordsize,
                           type_nuc=args.type, logger=logger,
                           rc_option=args.reverse_complement_search,
                           convert_wobbles=args.wobble_conversion,
                           max_N_percentage=defaults["max_N_percentage"],
                           substitution_count=args.substitution_count))

        t1 = flexiutil.time_track(t1, logger, step='Sequence regex matching')
        # logger.debug(f'X list: {x_lists}')

        # shade annotated regions
        gff_shade_list = []
        if args.input_gff_files is not None and len(args.input_gff_files) != 0:
            if seq_name in feat_dict.keys():
                features = feat_dict[seq_name]
                for item in features:
                    feat_type, start, stop = item
                    feat_color, strength, zoom = color_dict[feat_type]
                    start = max(0, start - zoom - 0.5)
                    stop = min(length_seq + 1, stop + zoom + 0.5)
                    width = stop - start
                    gff_shade_list.append(tuple([feat_type, start, width, feat_color, strength, zoom]))

        # collect lines
        lines = []
        color_list = []
        for (x_lines, y_lines, col) in [(x_lists_rc, y_lists_rc, args.line_col_rev),
                                        (x_lists, y_lists, args.line_col_for)]:
            if col != "white" and len(x_lines) != 0:
                for ldx in range(len(x_lines)):
                    try:
                        lines.append([(x_lines[ldx][0], y_lines[ldx][0]), (x_lines[ldx][-1], y_lines[ldx][-1])])
                        color_list.append(col)
                    except Exception as e:
                        logger.info(f"Proceeding after error with line collection at index {ldx} "
                                    f"(data={x_lines[ldx]}) in Selfdotplot"
                                    f"\n{e}")
        color_list = np.array(color_list)
        # logger.debug(f'\nLines: {lines}'
        #             f'\nColors: {color_list}')
        # TODO: multi-threading (end)

        t1 = flexiutil.time_track(t1, logger, step='Sequence matching')

        # plotting with matplotlib
        #################################

        # combined plotting
        if multi:

            # plotting subplot with matplotlib
            ax = plt.subplot(nrows, mcols, counter)  # rows, columns, plot number

            # gff-based shading of annotated regions
            for item in gff_shade_list:
                feat_type, start, width, feat_color, strength, zoom = item
                ax.add_patch(patches.Rectangle((start, start),  # (x,y)
                                               width, width,  # width, height
                                               edgecolor=None, linewidth=args.line_width + zoom,
                                               fill=True, facecolor=feat_color,
                                               alpha=strength))

            # draw lines
            lc = cllct.LineCollection(lines, colors=color_list, linewidths=args.line_width)
            ax.add_collection(lc)

            # format axes
            # print(P.xticks()[0], P.yticks()[0])
            plt.axis('scaled')  # make images quadratic
            plt.xlim(0, length_seq + 1)
            if args.mirror_y_axis:
                plt.ylim(0, length_seq + 1)  # rotate y-axis (point upwards)
            else:
                plt.ylim(length_seq + 1, 0)  # rotate y-axis (point downwards)
            plt.xlabel(f"[{aa_bp_unit}]", fontsize=args.label_size)
            plt.ylabel(f"[{aa_bp_unit}]", fontsize=args.label_size)
            plt.tick_params(axis='both', which='major', labelsize=args.label_size * .9)

            # # use same tick labels for x and y-axis
            # tick_locs, tick_labels = P.yticks()
            # P.xticks(tick_locs)
            # P.xlim(0, length_seq+1)

            plt.title(flexiutil.unicode_name(
                flexiutil.shorten_name(name_seq,
                                     max_len=args.title_length[0],
                                     title_clip_pos=args.title_length[1])),
                fontsize=args.label_size, fontweight='bold')
            # P.title(unicode_name(name_seq), fontsize=label_size*1.3, fontweight='bold')

            # save figure and re-initiate if page is full
            if counter == mcols * nrows:

                # finalize layout - margins & spacing between plots
                try:
                    plt.tight_layout(h_pad=.02, w_pad=.02)
                except ValueError or RuntimeError or TypeError or AttributeError as e:
                    logger.error(f"{e}\nAttention - pylab.tight_layout failed! "
                                 "Please check sequence names and layout settings!")
                plt.subplots_adjust(hspace=0.5, wspace=0.5)  # space between rows - def 0.4

                # save figure
                if args.filetype == "PDF":
                    plt.suptitle(f"Page {page_counter:03d}", size=args.label_size * 1.5, weight="bold", y=1.02)
                    pdf.savefig(bbox_inches='tight')
                else:
                    # name and create output files (names derived from SEQNAME)
                    fig_name = f"{prefix}{name_graph}_ws{args.wordsize}{suffix}-{page_counter:03d}.{args.filetype}"
                    list_of_fig_names.append(fig_name)
                    plt.savefig(fig_name, bbox_inches='tight')

                plt.close()
                plt.cla()

                counter = 0
                page_counter += 1

                fig = plt.figure(figsize=(figsize_x, figsize_y))

        # plotting separate figure files
        else:  # not multi

            fig = plt.figure(figsize=(args.plot_size, args.plot_size))  # figure size needs to be a square
            ax = plt.subplot(1, 1, 1)  # rows, columns, plotnumber

            # gff-based shading of annotated regions
            for item in gff_shade_list:
                feat_type, start, width, feat_color, strength, zoom = item
                ax.add_patch(patches.Rectangle((start, start),  # (x,y)
                                               width, width,  # width, height
                                               edgecolor=None, linewidth=args.line_width + zoom,
                                               fill=True, facecolor=feat_color,
                                               alpha=strength))

            # draw lines
            lc = cllct.LineCollection(lines, colors=color_list, linewidths=args.line_width)
            ax.add_collection(lc)

            # format axes
            plt.axis('scaled')  # make images quadratic
            plt.xlim(0, length_seq + 1)
            if args.mirror_y_axis:
                plt.ylim(0, length_seq + 1)  # rotate y-axis (point upwards)
            else:
                plt.ylim(length_seq + 1, 0)  # rotate y-axis (point downwards)
            plt.xlabel(f"[{aa_bp_unit}]", fontsize=args.label_size)
            plt.ylabel(f"[{aa_bp_unit}]", fontsize=args.label_size)
            plt.tick_params(axis='both', which='major', labelsize=args.label_size * .9)

            # # use same tick labels for x and y axis
            # tick_locs, tick_labels = P.yticks()
            # P.xticks(tick_locs)
            # P.xlim(0, length_seq+1)

            plt.title(flexiutil.unicode_name(
                flexiutil.shorten_name(name_seq,
                                     max_len=args.title_length[0],
                                     title_clip_pos=args.title_length[1])),
                fontsize=args.label_size * 1.3, fontweight='bold')

            # save figure
            if args.filetype == "PDF":
                plt.suptitle(f"Page {counter:03d}", size=args.label_size * 1.5, weight="bold", y=1.02)
                pdf.savefig(bbox_inches='tight')
            else:
                # name and create output files (names derived from SEQNAME)
                fig_name = (f"{prefix}{name_graph}-{counter}_"
                            f"{flexiutil.shorten_name(name_seq, max_len=args.title_length[0], title_clip_pos=args.title_length[1])}"
                            f"_ws{args.wordsize}{suffix}.{args.filetype}")

                list_of_fig_names.append(fig_name)
                plt.savefig(fig_name, bbox_inches='tight')

            plt.close()
            plt.cla()  # clear any prior graph

        t1 = flexiutil.time_track(t1, logger, step='Drawing')

    if multi and counter >= 1:
        # finalize layout - margins & spacing between plots
        try:
            plt.tight_layout(h_pad=.02, w_pad=.02)
        except ValueError or RuntimeError or TypeError or AttributeError as e:
            logger.error(f"{e}\nAttention - pylab.tight_layout failed! "
                         "Please check sequence names and layout settings!")
        plt.subplots_adjust(hspace=0.5, wspace=0.5)  # space between rows - def 0.4

        # save figure
        if args.filetype == "PDF":
            plt.suptitle(f"Page {page_counter:03d}", size=args.label_size * 1.5, weight="bold", y=1.02)
            pdf.savefig(bbox_inches='tight')
        else:
            # name and create output files (names derived from SEQNAME)
            fig_name = (f"{prefix}{name_graph}"
                        f"_ws{args.wordsize}{suffix}-{page_counter:03d}.{args.filetype}")
            list_of_fig_names.append(fig_name)
            plt.savefig(fig_name, bbox_inches='tight')

        plt.close()
        plt.cla()  # clear any prior graph

    if args.filetype == "PDF":
        pdf.close()

    logger.info(f"Drawing selfdotplots done.")

    seq_dict.close()

    return list_of_fig_names
