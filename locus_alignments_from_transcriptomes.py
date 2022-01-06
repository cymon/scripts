#! /usr/bin/env python3

desc = """
Script published as part of manuscript:
Almeida et al. A low-latitude species pump: Peripheral isolation,
parapatric speciation and mating-system evolution converge in a marine radiation
Molecular Ecology, in submission.

Given a bunch of assembled transcriptomes mapped with similarly identified loci,
build the locus alignments

Assumes all data is in fasta format
Assumes all sequence record labels are unique
Assumes sequence records are labeled:
    >identifier|locus

    - the identifier containing all the necessary labels for down-stream
      analysis
    - the locus being the unique locus number to identify similar loci in other
      transcriptome files

Two files are used to define the input:
    1) a file listing locus identifiers (one per line)
    2) a file listing filenames of the fasta data files (one per line)

Note that this has an awful hack because the locus names can sometimes contain
"/" which obviously causes problems for naming files: here they are converted
to "-" (hypen). Downstream scripts need to take this into account.

Cymon J. Cox - 6-1-2022
"""

import os
import sys
import shutil
import argparse
import textwrap
import subprocess
from Bio import SeqIO

def parse_input_file(input_file):
    """Read the items from the input file"""

    with open(input_file, 'r') as f:
        items = f.read().splitlines()
    if len(items) > len(set(items)):
        print("\tError: not all names in %s are unique")
        print("\tExiting... done.")
        sys.exit()
    else:
        return items

def main(datafilelist, locusnames, outdir="locus_alignments", csv_log=False, quiet=False):
    """Main loop"""

    loci = parse_input_file(locusnames)
    data = parse_input_file(datafilelist)

    if not quiet:
        print("\n")
        print("\tRead %i data files and %s loci" % (
            len(data), len(loci)))
        print("\n")

    for f in data:
        if not os.path.exists(f):
            print("\n")
            print("\tError: cannot find specified data file %s" % f)
            print("\tExiting... done.")
            sys.exit()

    if os.path.isdir(outdir):
        print("\n")
        print("\tError: output directory \"%s\" already present" % outdir)
        print("\tExiting... done.")
        sys.exit()
    else:
        if not quiet:
            print("\tMaking output directory %s..." % outdir,)
        os.makedirs(outdir)
        if not quiet:
            print("done")

    #INDEX
    if os.path.exists("locus.idx"):
        os.remove("locus.idx")
    if not quiet:
        print("\tMaking index file... ", end='')
    lind = SeqIO.index_db("locus.idx", data, "fasta")
    if not quiet:
        print("done")

    if csv_log:
        clog = ["locus,#alleles,max len,min len\n"]

    #Loop over each locus and write seqs to alignment file
    if not quiet:
        print("\tWriting locus files (number, and min and max lengths of alleles):")
    for locus in loci:
        if not quiet:
            print("\t\t%-8s -> " % locus, end='')
        locus_recs = []
        found = False
        for rec in lind:
            rec_locus = rec.split("|")[1]
            if rec_locus == locus:
                locus_recs.append(lind[rec])
                found = True
        if not found:
            print("\t Unable to find locus %s in data file %s" % (locus, rec))

        #Stats
        rec_lens = []
        for r in locus_recs:
            rec_lens.append(len(r))
        if not quiet:
            print("\t%-2i alleles (%-4i - %-4i)" % (len(locus_recs), min(rec_lens),
                    max(rec_lens)))
        if csv_log:
            clog.append("%s,%i,%i,%i\n" % (locus, len(locus_recs),
                min(rec_lens), max(rec_lens)))

        ## BEWARE HACK!!! Get rid of chars causing problems with file names
        if "/" in locus:
            locus = locus.replace("/", "-")

        SeqIO.write(locus_recs, "%s/%s.fasta" % (outdir,locus), "fasta")
        sys.stdout.flush()

    if csv_log:
        fh = open("locus_log.csv", "w")
        for l in clog:
            fh.write(l)
        fh.close()
    if not quiet:
        print("\n\tDone.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(desc),
            )
    parser.add_argument(dest="data_files",
                        help="File listing input data files (one per line)"
                        )
    parser.add_argument(dest="locus_names",
                        help="File listing locus names (one per line)"
                        )
    parser.add_argument("-d", "--output_dir",
                        dest="outdir",
                        help="Output directory name Default: locus_alignments",
                        default="locus_alignments"
                        )
    parser.add_argument("-i", "--reuse_index",
                        dest="reuse_index",
                        help="Reuse an old index file (named locus.idx)",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-o", "--overwrite",
                        dest="overwrite_index",
                        help="Overwrite index file",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-l", "--write_csv",
                        dest="csv_log",
                        help="Write a csv formatted log file",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-q", "--quiet",
                        dest="quiet",
                        help="Quiet - do not write progress. Default: False",
                        default=False,
                        action='store_true'
                        )
    args = parser.parse_args()
    main(args.data_files, args.locus_names, args.outdir, args.reuse_index,
            args.overwrite, args.csv_log, args.quiet)
