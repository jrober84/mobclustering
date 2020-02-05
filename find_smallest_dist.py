
import pandas as pd
from argparse import (ArgumentParser, FileType)

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(
        description="Filter fasta file with list of contigs")
    parser.add_argument('-i', '--infile', type=str, required=True, help='Input mash file')

    return parser.parse_args()



def main():
    args = parse_args()
    input_data = pd.read_csv(args.infile,header=None,sep="\t")
    max_dist = 1
    for i,row in input_data.T.iteritems():
        query = row[0]
        ref = row[1]
        dist = row[2]
        if query == ref:
            continue
        if dist < max_dist:
            max_dist = dist
    print("{}\tSmallest distance is:\t{}".format(args.infile,max_dist))




# call main function
if __name__ == '__main__':
    main()