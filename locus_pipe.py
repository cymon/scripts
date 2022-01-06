#! /usr/bin/env python3

desc = """
Script published as part of manuscript:
Almeida et al. A low-latitude species pump: Peripheral isolation,
parapatric speciation and mating-system evolution converge in a marine radiation
Molecular Ecology, in submission.

A pipeline for aligning and preparing the loci from transcriptome data for
phylogenetic analysis. Multiple transcriptome experiments, one fasta file per
taxon, are sorted into individual locus data files, aligned with TranslatorX
(inc. GBlocks), and are cleaned.

1. Write locus data files (option -w)
2. Align the locus data files with TranslocatorX (incl Gblocks)
3. Clean to remove any further ambiguities, gaps or stop codons at any codon
sites from the aligned data

 this leaves only solid blocks of in-frame codon aligned data with no gaps,
 ambigs, or stop codons anywhere

data_filenames is a list of the fasta files for each taxon experiment
locus_names is a list of locus names in each taxon experiment - these are the
string after the first "|" character in the fasta entry description ">":

The taxon experiment Mang_BCO-1_0 has the following sequence entry:
e.g.

>Mang_BCO-1_0|NODE_10010_length_3497_cov_20323.525232_g4438_i0_2394_3497_+

the locus name is "NODE_10010_length_3497_cov_20323.525232_g4438_i0_2394_3497_+"
(without the quotes)

Note that this has an awful hack because the locus names can sometimes contain
"/"s which causes problems for naming files: here they are converted to "-"
(hypen). Downstream scripts need to take this into account. It is best to remove
these characters from filenames before running this script.

Cymon J. Cox - 6-1-2022
"""
import argparse
import sys
import subprocess
import os
from Bio import SeqIO
import textwrap
import p4
p4.var.doCheckForBlankSequences = False
p4.var.doCheckForDuplicateSequences = False

from locus_alignments_from_transcriptomes import main as write_locus_alignments
from TranslatorX_functions import run_TranslatorX
from remove_gapped_and_ambiguous_codons import main as clean_loci

##############################################
SUPERFLUOUS_TRANSX_FILES = [
"%s_transX.aa_ali.fasta",
"%s_transX.aa_ali.fasta-gb.txts",
"%s_transX.aa_based_codon_coloured-gb.html",
"%s_transX.aa_based_codon_coloured.html",
"%s_transX.aa_cleanali.fasta",
"%s_transX.aaseqs",
"%s_transX.aaseqs.fasta",
"%s_transX.html",
#"%s_transX.mafft.log", - signal for already having attempted locus
"%s_transX.nt12_ali.fasta",
"%s_transX.nt12_cleanali.fasta",
"%s_transX.nt1_ali.fasta",
"%s_transX.nt1_cleanali.fasta",
"%s_transX.nt2_ali.fasta",
"%s_transX.nt2_cleanali.fasta",
"%s_transX.nt3_ali.fasta",
"%s_transX.nt3_cleanali.fasta"]
##############################################

def main(data_filenames, locus_names, out_dir, write_loci=False):

    ## STEP 1: Write the individual locus files from the sample fasta files
    if write_loci:
        write_locus_alignments(data_filenames, locus_names, outdir=out_dir, csv_log=True, quiet=False)
    else:
        print(f"\nNot writing locus data files, expecting them to be present in the", end="")
        print(f" directory '{out_dir}'\n")

    ## STEP 2: Align with TranslatorX
    f = open(locus_names, 'r')
    loci = f.read().splitlines()
    f.close()

    os.chdir(out_dir)
    for locus in loci:
        sys.stdout.flush()

        ## BEWARE HACK!!! Get rid of chars causing problems with file names
        if "/" in locus:
            locus = locus.replace("/", "-")

        #See if we have already done this locus
        if os.path.exists("%s_transX.mafft.log" % locus):
            print("\tSkipping %s - already done" % locus)
            continue

        #Check file exists and has content
        if not os.path.exists("%s.fasta" % locus):
            print(f"Cannot find locus data file: {locus}.fasta. Exiting")
            sys.exit()
        if os.path.getsize(f"{locus}.fasta") == 0:
            print(f"Locus file {locus}.fasta has 0 bytes content. Exiting")
            sys.exit()

        print("\tAligning %s with TranslatorX " % locus, end="")
        returncode, stderr, stdout = run_TranslatorX(f"{locus}.fasta", code=1, gbl="-b5=n")
        #if r[0] == "":
        #    print("\tTranslatorX failed - perhaps cannot find input file...")
        #    print("\n")
        #    continue
        #else:
        #(returncode, stderr, stdout) = r

        if returncode != 0:
            if b"Gblocks alignment:  0 positions" in stderr:
                print("\tGblocks alignment:  0 positions")
            else:
                print("\tError running TranslatorX - see transX_stdout.text...")
                print("\treturncode = %s" % returncode)
                print("\tSee transX_stdout.text. Continuing...")
                of = open("transX_stdout.text", "a")
                of.write(stdout)
                of.write(stderr)
                of.close()
            #Remove these extra file that are not in SUPERFLUOUS_TRANSX_FILES:
            extra_files = ["%s_transX.aa_ali.fasta-gb.htm",
                           "%s_transX.nt_ali.fasta",
                           "%s_transX.nt_cleanali.fasta"]
            for f in extra_files:
                os.remove(f % locus)
        else:
            target_file = "%s_transX.nt_cleanali.fasta" % locus
            if not os.path.exists(target_file):
                print("\tCannot find %s. Continuing..." % target_file)
                continue
            # Clean out all ambigous, gapped, or stop codons if present - even when
            # not constant - by removing the codon site (3 nucs at which it occurs)
            clean_loci(target_file, quiet=True)

        #Finally clean up these file!
        for f in SUPERFLUOUS_TRANSX_FILES:
            os.remove(f % locus)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(desc),
            )
    parser.add_argument(dest="data_filenames",
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
    parser.add_argument("-w", "--write_loci",
                        dest="write_loci",
                        help="Write the locus data files. Default: False",
                        default=False,
                        action='store_true'
                        )
    args = parser.parse_args()
    main(args.data_filenames, args.locus_names, args.outdir, args.write_loci)
