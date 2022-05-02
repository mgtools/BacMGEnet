import sys, os
from pathlib import Path
import pandas as pd
import argparse
import Bio.SeqIO as IO

def to_dict_remove_dups(sequences):
    ''' 
    resolve dups in fasta dictionary 
    '''
    return {record.id: record for record in sequences}

def get_MGE_headers(gff):
    
    MGE_set = set()

    with open(gff) as infile:
        for line in infile:
            line = line.rstrip('\n')
            sline = line.split('\t')
            if len(sline) < 2:
                continue
            MGE_set.add(sline[0])

    return(list(MGE_set))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--phagedb", help = "phage db fasta file", default = '/home/team/phage_databases/all_phage_seq.fna')
    parser.add_argument("--plasmiddb", help = "plasmid db fasta file", default ='/home/team/plasmid_databases/all_plasmid_seq.fna')
    parser.add_argument("--gff", help = "input gff file")
    parser.add_argument("--outdir", help = 'output directory')
    args = parser.parse_args()

    # build dictionary of fasta file
    phage_fasta = to_dict_remove_dups(IO.parse(args.phagedb, "fasta"))
    plasmid_fasta = to_dict_remove_dups(IO.parse(args.plasmiddb, "fasta"))
    
    # merge dictionaries
    db_fasta = {**phage_fasta, **plasmid_fasta}

    # get MGE headers from gff file
    MGE_list = get_MGE_headers(args.gff)

    # make output directory if doesnt exist
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    for mge in MGE_list:
        mge_name = mge.split('|')[0]
        with open(f'{args.outdir}/{mge_name}.fna','w') as outfile:
            print(f'exporing {mge} to fasta')
            outfile.write(f'>{mge_name}\n')
            outfile.write(f'{db_fasta[mge].seq}\n')
    #print(db_fasta[MGE_list[0]].seq)
