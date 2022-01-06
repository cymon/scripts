#! /usr/bin/env python3

desc = """
Script published as part of manuscript:
Almeida et al. A low-latitude species pump: Peripheral isolation,
parapatric speciation and mating-system evolution converge in a marine radiation
Molecular Ecology, in submission.


A CODON site is the in-frame 3 nucleotides specifiying 1,2,and 3rd codon
positions that are homologous across all sequences. This script identifies and
can removes CODON sites (sets of 3 nucleotides).

By default, CODON sites which contain a gap, ambiguity, or stop codon of the
specified genetic code are removed (default universal code). By default, this
removes the entire CODON site (all three nucs) if only 1 site in 1 sequence is
gapped, ambiguous or a stop codon. To remove a CODON site only if all sequences
have a a gap, ambiguity, or stop codon at the same site in ALL sequences, use
the -c(onstant) option.

Options can be used to remove only, or combinations, of -g(apped), -a(mbiguous),
and -s(top codons); if none of -g,-a,-s are specified, default is to remove any
site with any of the three.

Cymon J. Cox 6-1-2022
"""

import os
import sys
import time
import argparse
import textwrap
from p4 import *
from p4 import geneticcode
var.doCheckForDuplicateSequences = False
var.doCheckForBlankSequences = False

def contains_ambig(codon):
    for i in codon:
        if i not in ["a", "c", "g", "t"]:
            return True
    return False

def contains_gap(codon):
    for c in codon:
        if "-" in c:
            return True
    return False

def is_stop(codon, stop_codons):
    if codon in stop_codons:
        return True
    return False

def main(matrix, code=1, ambig=False, gaps=False,
        stops=False, constant=False, quiet=False):

    c = geneticcode.GeneticCode(transl_table=code)

    stop_codons = c.codonsForAA["*"]

    read(matrix)
    a = var.alignments[0]
    if a.nChar == 0 or (a.nChar%3) != 0:
        print("Error: Matrix length %s is not exactly divisable by 3." %
                a.nChar)
        print("The assumption is made that the alignment is in-frame codons")
        #Clean up P4 vars in case it's being looped
        var.alignments = []
        var.nexusSets       = None
        var.trees           = []
        var.sequenceLists   = []
        return
    else:
        no_codons = a.nChar/3

    if not quiet:
        print("\tAlignment has %i characters and %i CODONS..." % (a.nChar,
                no_codons))

    codons_by_taxon =  [[seq.sequence[i:i+3] for i in
            range(0, len(seq.sequence), 3)] for seq in a.sequences]
    codon_sites = zip(*codons_by_taxon)

    valid_codon_sites = []
    count = 0
    count_ambigs = 0
    count_gaps   = 0
    count_stops  = 0

    #Default
    #If not one of -g, -a, or -s, remove all gaps ambigs and stops
    if not ambig and not gaps and not stops:
        ambig = True
        gaps  = True
        stops = True

    if not quiet:
        if constant:
            print("\tOnly removing codon sites if all sequences contain target...")
        else:
            print("\tRemoving codon sites if any sequence has target...")
        if gaps:
            print("\tRemoving gaps...")
        else:
            print("\tNOT removing gaps...")
        if ambig:
            print("\tRemoving ambiguous...")
        else:
            print("\tNOT removing ambiguous...")
        if stops:
            print("\tRemoving stop codons...")
        else:
            print("\tNOT remvoving stop codons...")
        print("\tProcessing codon sites: ", end="")

    #Constant
    if not constant:
        #remove either or gaps, ambigs, stops
        for codon_site in codon_sites:
            for codon in codon_site:
                if ambig:
                    if contains_ambig(codon):
                        count_ambigs += 1
                        break
                if gaps:
                    if contains_gap(codon):
                        count_gaps += 1
                        break
                if stops:
                    if is_stop(codon, stop_codons):
                        count_stops += 1
                        break
            else:
                valid_codon_sites.append(codon_site)
                if not quiet:
                    print(".", end="")
                count += 1
    #Per seq
    else:
        for codon_site in codon_sites:
            if ambig:
                if [contains_ambig(codon) for codon in
                        codon_site].count(True) == len(codon_site):
                    count_ambigs += 1
                    continue
            if gaps:
                if [contains_gap(codon) for codon in
                        codon_site].count(True) == len(codon_site):
                    count_gaps += 1
                    continue
            if stops:
                if [is_stop(codon, stop_codons) for codon in
                        codon_site].count(True) == len(codon_site):
                    count_stops += 1
                    continue
            valid_codon_sites.append(codon_site)
            if not quiet:
                print(".",)
            count += 1

    if not quiet:
        print
        if ambig:
            print("\tRemoved %-3i sites containing ambiguities..." %
                    count_ambigs)
        if gaps:
            print("\tRemoved %-3i sites containing gaps..." % count_gaps)
        if stops:
            print("\tRemoved %-3i sites containing stop codons..." %
                    count_stops)

    counted_codons = count+count_gaps+count_ambigs+count_stops
    assert counted_codons == no_codons, "Error: not all codons counted (%i != %i)" % (counted_codons, no_codons)
    valid_codons_by_taxon = zip(*valid_codon_sites)
    valid_sequences = ["".join(seq) for seq in valid_codons_by_taxon]

    if valid_sequences == []:
        print("\tAll codon sites removed; no remaining sites - no output.")
    else:
        for i, seq in enumerate(valid_sequences):
            a.sequences[i].sequence = seq
        a.checkLengthsAndTypes()
        (filename, ext) = os.path.splitext(a.fName)
        fn = filename + "_cleaned.nex"
        a.writeNexus(fn)
        print("\tWritten cleaned alignment (%i codons) to: %s" % \
            (len(valid_codon_sites), fn))

    #Clean up P4 vars in case it's being looped
    var.alignments = []
    var.nexusSets       = None
    var.trees           = []
    var.sequenceLists   = []
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(desc),
            )
    parser.add_argument("alignment",
                        help="Formatted codon alignment (Nexus, PHYLIP, FASTA)"
                        )
    parser.add_argument("-t", "--translation_code",
                        dest="code",
                        help="The number of the genetic code translation " + \
                             "table http://www.ncbi.nlm.nih.gov/Taxonomy/" + \
                             "Utils/wprintgc.cgi?mode=t#SG11 Default: 1 " +\
                             "(standard)",
                        type=int,
                        default=1
                        )
    parser.add_argument("-a", "--ambiguous",
                        dest="ambig",
                        help="Remove (only) codon sites that have amiguous sites. Default: False",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-g", "--gaps",
                        dest="gaps",
                        help="Remove (only) gapped codons. Default: False",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-s", "--stops",
                        dest="stops",
                        help="Remove (only) codon sites that have stop codons. Default: False",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-c", "--constant",
                        dest="constant",
                        help="Only remove if all sequences have the target. Default: False",
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
    main(args.alignment, args.code, args.ambig, args.gaps, args.stops, args.constant, args.quiet)
    print("\n\tDone\n")

