#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import itertools
import argparse
from Bio import SeqIO

def main():
    ####################
    # read and check fasta file
    ####################
    
    # read fasta file => create fasta dict : key = sequence name, value = sequence string
    fasta = fasta2hash(args.fasta)

    # check all sequences have same size
    length = 0
    for seqid in fasta:
        if length == 0:
            length = len(fasta[seqid])
        else:
            if len(fasta[seqid]) != length:
                print("ERROR: fasta input sequences do not have identical lengths. seq {} {}bp # {}bp".format(seqid, len(seqid), length))
                exit(1)
                

    ####################
    # iter on positions to count singletons, multinuc and to collect nucleotids per seqid
    ####################
    
    # sitesHoH[pos][seqid] = nucl
    sitesHoH = dict()
    # multinuc[pos] = 1
    multinuc = dict()
    # singleton[pos] = 1
    singletons = dict()

    for pos in range(length):
        # pos start is 0-index, print_pos is 1-index
        print_pos = pos + 1
        
 	# count nucleotides at that position for all seq => dict nuc_hash : key = nuc [ATCG], value = number of seqids with this nuc
        # nuc_hash[NUC] = sum
        nuc_hash = dict()
        for seqid in fasta: #seqid = ind's names
            nuc = fasta[seqid][pos].upper() #nuc = the nucleotides [ATCGN] - pos = the positions in the seq
            if nuc not in nuc_hash:
                nuc_hash[nuc] = 0
            nuc_hash[nuc] += 1

        # if polymorphe position with exactly 2 alleles
        #   skip singleton, after tracing theses pos in singleton dict
        #   keep others by storing NUCLEOTID in sitesHoH[POS][SEQID]
        # else
        #  this is an error

 	# skip monomorphe positions
        if len(nuc_hash) == 1:
            if args.verbose: 
                print("skip pos {} : monomorph, nuc_hash={}".format(print_pos, nuc_hash))
            continue

        # skip polymorphe positions with more than 3 possible alleles, after tracing these pos in multinuc dict
        if len(nuc_hash) > 3 :
            multinuc[print_pos] = 1
            if args.verbose: 
                print("skip pos {}, multinuc, nuc_hash={}".format(print_pos, nuc_hash))
            continue

        # skip polym positions with 3 alleles if none of them is a 'N'
        if len(nuc_hash) == 3 :
            #if one the allele = 'N' --> keep the position
            if 'N' in nuc_hash.keys() :
                #check for singleton : there is at least one nucleotid found in only one seq
                if min(nuc_hash.values()) == 1:
                    min_key = min(nuc_hash, key=nuc_hash.get)
                    #if the singleton is not a 'N' -->  skip the position
                    if min_key != 'N' :
                        # trace this pos to singletons
                        singletons[print_pos]= 1
                        if args.verbose: 
                            print("skip pos {}, singleton, nuc_hash={}".format(print_pos, nuc_hash))
                    else :
                        continue
                #if the singleton is a 'N' --> keep the pos & continue
                else:
                    for seqid in fasta:
                        nuc = fasta[seqid][pos].upper()
                        if print_pos not in sitesHoH:
                            sitesHoH[print_pos] = dict()
                        sitesHoH[print_pos][seqid] = nuc
                    if args.verbose: 
                        print("keep pos {}, nuc_hash={}, sitesHoH={}".format(print_pos, nuc_hash, sitesHoH[print_pos]))
            #if no 'N' in the 3 alleles, skip the position because 3 'real' alleles
            else :
                multinuc[print_pos] = 1
                if args.verbose: 
                    print("skip pos {}, multinuc, nuc_hash={}".format(print_pos, nuc_hash))
                continue
            continue
            
        # if exactly two alleles    
        if len(nuc_hash) == 2 :
            #if 1 of the 2 alleles = 'N' --> skip the pos
            if 'N' in nuc_hash.keys() :
                multinuc[print_pos] = 1
                if args.verbose: 
                    print("skip pos {}, multinuc, nuc_hash={}".format(print_pos, nuc_hash))
                continue
            else : #if 2 real alleles
                # search for singletons
                if min(nuc_hash.values()) == 1:
                    singletons[print_pos]= 1
                    if args.verbose: 
                        print("skip pos {}, singleton, nuc_hash={}".format(print_pos, nuc_hash))
                else:
 		    # store nucleotid for that position and that seqid
                    # create hash of hash : key1=pos, key2=seqid, val=nuc
                    for seqid in fasta:
                        nuc = fasta[seqid][pos].upper()
                        if print_pos not in sitesHoH:
                            sitesHoH[print_pos] = dict()
                        sitesHoH[print_pos][seqid] = nuc
                    if args.verbose: 
                        print("keep pos {}, nuc_hash={}, sitesHoH={}".format(print_pos, nuc_hash, sitesHoH[print_pos]))
        else:
            print("error pos {}, nuc_hash={}".format(print_pos, nuc_hash))
            print("\nERROR: nucleotide counting failed at site pos {}".format(print_pos))
            exit(2)

    ####################
    # print some global counts
    ####################
    singleton_count = len(singletons)
    multinuc_count = len(multinuc)
    parsimony_count = len(sitesHoH)
    polymorphic_count = singleton_count + multinuc_count + parsimony_count
    print("#Input File: {}".format(args.fasta))
    print("#Alignment Length: {}".format(length))
    print("#Polymorphic Sites: {}".format(polymorphic_count))
    print("#Singleton sites: {}".format(singleton_count))
    print("#Excluded sites with >3 alleles: {}".format(multinuc_count))
    print("#Sites analyzed for FGT: {}".format(parsimony_count))

    ####################
    # print details results for each pair of positions
    ####################
    
    # sorted list of sequence identifiers
    seqids = sorted(fasta)

    # print header line
    columns = ["Site1", "Site2"] + seqids + ["AlleleCount", "AlleleCount_noN"]
    print("\t".join(columns))

    # print other lines, one line per pair of positions
    # sorted list of positions with sitesHoH values
    pos_array = sorted(sitesHoH)
    
    # for each (i, j) pair of positions
    # https://docs.python.org/3/library/itertools.html
    # permutations(['A', 'B', 'C', 'D'], 2) => AB AC AD BA BC BD CA CB CD DA DB DC
    # combinations(['A', 'B', 'C', 'D'], 2) => AB AC AD BC BD CD
    for i, j in itertools.combinations(pos_array, 2):
    #for i, j in itertools.permutations(pos_array, 2):
        # first 2 columns = i and j positions
        columns = [str(i), str(j)]

        # total number of sequences per dinuc
        dinuc_hash = dict()
        # total number of sequences per dinuc with no N
        dinuc_hash_noN = dict()

        # for each seqid
        for seqid in seqids:
            # collect dinuc for that seqid
            nuc1 = sitesHoH[i][seqid]
            nuc2 = sitesHoH[j][seqid]
            dinuc = "{}{}".format(nuc1, nuc2).upper()
            # add this dinuc to the seqid column 
            columns.append(dinuc)
            # per dinuc, count the number of seqids
            if dinuc not in dinuc_hash:
                dinuc_hash[dinuc] = 0
            dinuc_hash[dinuc] += 1
            # for each dinuc with no N, count the number of seqids
            if 'N' not in dinuc:
                if dinuc not in dinuc_hash_noN:
                    dinuc_hash_noN[dinuc] = 0
                dinuc_hash_noN[dinuc] += 1

        # last 2 columns : number of dinuc and dinuc_noN
        allele_count = len(dinuc_hash)
        allele_count_noN = len(dinuc_hash_noN)
        columns.append(str(allele_count))
        columns.append(str(allele_count_noN))

        # print all columns for this line corresponding to (i, j) pair of positions
        print("\t".join(columns))



def fasta2hash(in_fasta):
    """
    read fasta file and create a dict : key=seqid (sequence identifier), value = sequence (as a string)
    """
    fasta = dict()
    with open(in_fasta) as f_in:
        for record in SeqIO.parse(f_in, "fasta") :
            fasta[record.id] = str(record.seq)
    print("{} seq in {}".format(len(fasta), in_fasta))
    return fasta
    
if __name__ == '__main__':
    description = """ 
    Traduction of WHO ?? perl script to python script
    Modif : 
      - add count for alleles ignoring 'N' nucleotide
      - skip lines 
          for pairs for identical positions
          for pairs of positions (j, i), redondant with lines for pairs of positions (i, j)
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--verbose", action="store_const", const="-v", default="")
    parser.add_argument("fasta", help="input fasta_file")
    args = parser.parse_args()
    main() 
