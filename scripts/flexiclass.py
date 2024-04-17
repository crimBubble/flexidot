import argparse
import os
import glob
from scripts import flexiutil


# Custom formatter for -h/--help
class CustomFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


# Custom FileType to add recognition of "-i auto" and directory instead of file(s)
class CustomFileType(argparse.FileType):
    def __call__(self, string):
        print("reading input")
        if string == "auto":
            # Get the current working directory
            cwd = os.getcwd()

            # List of file extensions
            extensions = helpers.fasta_extensions

            # Use a list comprehension to create the list of matching file paths
            matching_files = [file for ext in extensions for file in glob.glob(os.path.join(cwd, '*' + ext))]

            # Return the list of opened file objects
            return [super(CustomFileType, self).__call__(file) for file in matching_files]

        elif os.path.isdir(string):
            # List of file extensions
            extensions = [".fasta", ".fa", ".fna", ".faa", ".fas"]

            # Use a list comprehension to create the list of matching file paths
            matching_files = [file for ext in extensions for file in
                              glob.glob(os.path.join(string, '*' + ext))]

            print("reading")

            # Return the list of opened file objects
            return [super(CustomFileType, self).__call__(file) for file in matching_files]

        else:
            return super(CustomFileType, self).__call__(string)


# Custom boolean storage to set True and False more versatile and make store_true as default viable
class CustomStoreBooleanAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values.lower() in ("t", "true", "y", "yes", "1", "j", "on"):
            setattr(namespace, self.dest, True)
        elif values.lower() in ("f", "false", "n", "no", "0", "n", "off"):
            setattr(namespace, self.dest, False)
        else:
            parser.error(f"Invalid value for {option_string}: {values}")


class StoreIfSet(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)
        setattr(namespace, f"{self.dest}_set", True)


class TitleLengthAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):

        if len(values) == 1:
            value = values[0]
            if value.upper() in ['B', 'E']:
                length = 20
                symbol = value.upper()
                setattr(namespace, self.dest, [length, symbol])
            else:
                try:
                    length = int(value)
                    if length >= 0:
                        setattr(namespace, self.dest, [length, 'B'])
                    else:
                        parser.error(f"Invalid -T/--title_length input '{value}'. "
                                     f"Title length needs to be a positive integer. Set as 'int B/E'.")
                except ValueError or TypeError as e:
                    parser.error(f"{e}. Invalid -T/--title_length input '{value}'. Set as 'int B/E'.")
        elif len(values) == 2:
            try:
                length = int(values[0])
                if length < 0:
                    parser.error(f"Invalid -T/--title_length input '{length}'. "
                                 f"Title length needs to be a positive integer. Set as 'int B/E'.")
            except ValueError or TypeError as e:
                parser.error(f"{e}. Invalid -T/--title_length input length value '{values[0]}'. Set as 'int B/E'.")

            symbol = values[1].upper()
            if symbol not in ['B', 'E']:
                parser.error(f"Invalid -T/--title_length input '{symbol}'. "
                             f"Title symbol needs to be 'B' or 'E'. Set as 'int B/E'.")
            setattr(namespace, self.dest, [length, symbol])
        else:
            raise argparse.ArgumentTypeError("Invalid -T/--title_length input. Set as 'int B/E'")




