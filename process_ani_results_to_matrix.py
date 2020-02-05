from Bio import SeqIO
import pandas as pd
from argparse import (ArgumentParser, FileType)
import os, sys
def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(
        description="Filter fasta file with list of contigs")
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Output file')
    parser.add_argument('-i', '--infile', type=str, required=True, help='Input ANI results')
    return parser.parse_args()

def main():
    args = parse_args()
    input_data = args.infile
    outfile = args.outfile
    data = {}

    id_dict = {}

    fh = open(input_data,'r')
    for line in fh:
        row = line.strip().split("\t")
        query =row[0]
        ref = row[1]
        ani =  100 - float(row[2])
        matching_frags = row[3]
        total_frags = row[4]
        if not query in id_dict:
            id_dict[query] = ''
        if not ref in id_dict:
            id_dict[ref] = ''

        if float(matching_frags)/float(total_frags) < 0.5:
            ani = 100

        if not query in data and not ref in data:
            data[query] = {}

        if query in data and ref in data:
            if query in data and ref in data[query]:
                if data[query][ref] > ani:
                    data[query][ref] = ani
            elif ref in data and query in data[ref]:
                if data[ref][query] > ani:
                    data[ref][query] = ani
            continue


        if query in data:
            data[query][ref] = ani
        else:
            data[ref][query] = ani

    fh.close()
    print("#query\t{}".format("\t".join(id_dict.keys())))
    for i1 in id_dict:
        pos = 0
        dists = [100] * len(id_dict)
        for i2 in id_dict:
            if i1 == i2:

                dists[pos] = 0
            else:
                if i1 in data:
                    if i2 in data[i1]:
                        dists[pos] = data[i1][i2]
                    elif i2 in data:
                        if i1 in data[i2]:
                            dists[pos] = data[i2][i1]
                elif i2 in data:
                    if i1 in data[i2]:
                        dists[pos] = data[i2][i1]
                    elif i1 in data:
                        if i2 in data[i1]:
                            dists[pos] = data[i1][i2]
                else:
                    dists[pos] = 100
            pos+=1


        print("{}\t{}".format(i1,"\t".join(str(v) for v in dists)))






# call main function
if __name__ == '__main__':
    main()