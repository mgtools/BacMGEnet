import sys, os
import pandas as pd
import argparse
import Bio.SeqIO as IO

def to_dict_remove_dups(sequences):
    '''
    resolve dups in fasta dictionary 
    '''
    return {record.id: record for record in sequences}

def get_protospacer_dic(blast_infile):
    '''
    blast_dic[query_id] = [subject_id]
    protospacer_dic[subject_id] = [sstart:send]
    source_dic[protospacer_hit_id][start:end] = [source_id] 
    '''
    
    blast_dic = dict()
    protospacer_dic = dict()
    source_dic = dict()

    with open(blast_infile) as infile:
        for line in infile:
            line = line.rstrip('\n')
            sline = line.split('\t')
            qseqid = sline[0]
            sseqid = sline[1]
            pident = sline[2]
            qstart = sline[7]
            qend   = sline[8]
            sstart = sline[9]
            send   = sline[10]

            if qseqid not in blast_dic:
                blast_dic[qseqid] = list()
            blast_dic[qseqid].append(sseqid)

            if sseqid not in protospacer_dic:
                protospacer_dic[sseqid] = list()
            protospacer_dic[sseqid].append(f'{sstart}:{send}')
            protospacer_dic[sseqid] = list(set(protospacer_dic[sseqid]))

            if sseqid not in source_dic:
                source_dic[sseqid] = dict()
            if f'{sstart}:{send}' not in source_dic[sseqid]:
                source_dic[sseqid][f'{sstart}:{send}'] = list()
            source_dic[sseqid][f'{sstart}:{send}'].append(qseqid)
            
    return(blast_dic, protospacer_dic, source_dic)

def build_gff(phage_fasta, plasmid_fasta, phage_blast_dic, phage_protospacer_dic, plasmid_blast_dic, plasmid_protospacer_dic, phage_source_dic, plasmid_source_dic, outfile): 

    '''
    outputs gff file of CRISPR targets with 2 or more protospacers
    '''

    # get list of phage protospacer hits
    phage_list = list()
    for qseqid in phage_blast_dic.keys():
        phage_list += phage_blast_dic[qseqid]
    phage_list = set(phage_list)

    # get list of plasmid protospacer hits
    plasmid_list = list()
    for qseqid in plasmid_blast_dic.keys():
        plasmid_list += plasmid_blast_dic[qseqid]
    plasmid_list = set(plasmid_list)

    # output header
    outfile.write('##gff-version 3\n')

    # itterate over each phage with protospacer
    for phage in phage_list:
        
        # check if phage has multiple protospacer hits 
        #if len(phage_protospacer_dic[phage]) <= 1:
        #    continue

        outfile.write(f'{phage}\t.\tcontig\t1\t{len(phage_fasta[phage])}\t.\t+\t.\tID={phage};Type=phage;n_Protospacers={len(phage_protospacer_dic[phage])}\n')

        # print protospacer information
        for protospacer in phage_protospacer_dic[phage]:
            # get position
            protospacer_start, protospacer_end = protospacer.split(':')

            # determine strandadness
            strand = '+'
            if int(protospacer_start) > int(protospacer_end):
                strand = '-'
                protospacer_start, protospacer_end = protospacer_end, protospacer_start

            outfile.write(f'{phage}\t.\tprotospacer\t{protospacer_start}\t{protospacer_end}\t.\t{strand}\t.\tID={phage};Type=protospacer;SpacerSource={phage_source_dic[phage][protospacer]}\n')

    # itterate over each plasmid with protospacer
    for plasmid in plasmid_list:

        # check if phage has multiple protospacer hits 
        #if len(plasmid_protospacer_dic[plasmid]) <= 1:
        #    continue

        outfile.write(f'{plasmid}\t.\tcontig\t1\t{len(plasmid_fasta[plasmid])}\t.\t+\t.\tID={plasmid};Type=plasmid;n_Protospacers={len(plasmid_protospacer_dic[plasmid])}\n')

        # print protospacer information
        for protospacer in plasmid_protospacer_dic[plasmid]:
            # get position
            protospacer_start, protospacer_end = protospacer.split(':')

            # determine strandadness
            strand = '+' 
            if int(protospacer_start) > int(protospacer_end):
                strand = '-' 
                protospacer_start, protospacer_end = protospacer_end, protospacer_start

            outfile.write(f'{plasmid}\t.\tprotospacer\t{protospacer_start}\t{protospacer_end}\t.\t{strand}\t.\tID={plasmid};Type=protospacer;SpacerSource={plasmid_source_dic[plasmid][protospacer]}\n') 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--phagedb", help = "phage db fasta file", default = '/home/team/phage_databases/all_phage_seq.fna')
    parser.add_argument("--plasmiddb", help = "plasmid db fasta file", default ='/home/team/plasmid_databases/all_plasmid_seq.fna')
    parser.add_argument("--plasmid", help = "blast output for plasmids")
    parser.add_argument("--phage", help = "blast output for phage")
    parser.add_argument("--outdir", help = 'output directory')
    parser.add_argument("-l", "-t", "--tag", help="L# (e.g. L47, L35 and L29, and all)", default="")
    args = parser.parse_args()

    # build dictionary of fasta file
    phage_fasta = to_dict_remove_dups(IO.parse(args.phagedb, "fasta"))
    plasmid_fasta = to_dict_remove_dups(IO.parse(args.plasmiddb, "fasta"))
    
    # merge dictionaries
    db_fasta = {**phage_fasta, **plasmid_fasta}

    # get blast dictionary
    phage_blast_dic, phage_protospacer_dic, phage_source_dic = get_protospacer_dic(args.phage)
    plasmid_blast_dic, plasmid_protospacer_dic, plasmid_source_dic = get_protospacer_dic(args.plasmid)

    # get list of phage protospacer hits
    phage_list = list()
    for qseqid in phage_blast_dic.keys():
        phage_list += phage_blast_dic[qseqid]
    phage_list = set(phage_list)

    # get list of plasmid protospacer hits
    plasmid_list = list()
    for qseqid in plasmid_blast_dic.keys():
        plasmid_list += plasmid_blast_dic[qseqid]
    plasmid_list = set(plasmid_list)

    # build gff 
    with open(f'{args.outdir}/{args.tag}protospacer_w_multihits.gff', 'w') as outfile:
        build_gff(phage_fasta, plasmid_fasta, phage_blast_dic, phage_protospacer_dic, plasmid_blast_dic, plasmid_protospacer_dic, phage_source_dic, plasmid_source_dic, outfile)
