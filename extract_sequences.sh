#!/usr/bin/python3
import sys
import re
import argparse


# overloaded error method so it prints help info when user gives wrong arguments
class NewParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('ERROR: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def create_parser():
    """
    Creates a help parser (using custom argparser class)
    :return:
    """
    p = NewParser()

    p.add_argument('pattern', type=str,
                   help = "Regular expression to match headers with.")

    p.add_argument('fasta', type=str,
                   help = "Fasta file to extract sequences from.")

    p.description = 'This program extracts sequences whose headers match a pattern specified by the pattern argument.' \
                    '\nFor instance, the pattern "ribosomal RNA" will extract all sequences with ribosomal RNA in their header.'


    args = p.parse_args(sys.argv[1:])
    return args



if __name__ == '__main__':
    args = create_parser()
    f = open(args.fasta, "r")
    p = args.pattern.replace("(", "\\(").replace(")", "\\)")
    
    read_flag = False
    for line in f:
        if (line[0] == ">"):
            read_flag = re.search(p, line) is not None
    
        if read_flag:
            print(line.rstrip())

