import sys, os
import argparse

def get_crt_dic(crispralign, crt_dic):
    '''
    builds a dictionary of crispralign 
    '''
    with open(crispralign) as crt:
        for line in crt:
            line = line.rstrip('\n')
            # get header 
            if line.startswith('SEQ'):
                header = (line.split()[1])
            # get crispr range
            if line.startswith('CRISPR'):
                sline = line.split()
                crispr_start = sline[3]
                crispr_end = sline[5]

                # add to dictionary
                if header not in crt_dic:
                    crt_dic[header] = list()
                
                crt_dic[header].append(f'{crispr_start}:{crispr_end}')
                # remove duplicate entries
                crt_dic[header] = list(set(crt_dic[header]))

    return(crt_dic)

def filter_blast(blast_output, crt_dic, outfile):

    with open(blast_output) as m8:
        for line in m8:

            # true hit by default
            hit = True

            nline = line.rstrip('\n')
            sline = nline.split()
            qseqid = sline[0]
            sseqid = sline[1]
            pident = sline[2]
            qstart = sline[7]
            qend   = sline[8]
            protospacer_start = sline[9]
            protospacer_end   = sline[10]

            # check if blast hit matches to a crispr found within phage or plasmid db 
            if sseqid in crt_dic:

                # check if the blast hit matches the crispr region, check all crisprs predicted in this contig
                for crispr in crt_dic[sseqid]:
                    crispr_start, crispr_end = crispr.split(':')
        
                    # if blast hit falls within crispr region, mark hit as False positive
                    if crispr_start <= protospacer_start <= crispr_end and crispr_start <= protospacer_end <= crispr_end:
                        hit = False


            # subject sequence has no crispr predicted from crispralign
            if hit == True:
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--blast", help="blast results")
    parser.add_argument("--phage_crt",  help = "CRISPRalign results crt file - phage db")
    parser.add_argument("--plasmid_crt", help = "CRISPRalign results crt file - plasmid db")
    args = parser.parse_args()

    # get crt crispr range dictionary
    crt_dic = dict()
    crt_dic = get_crt_dic(args.plasmid_crt, crt_dic)
    crt_dic = get_crt_dic(args.phage_crt, crt_dic)

    # filter blast and save to output file
    with open(os.path.splitext(args.blast)[0]+".filtered.m8", 'w') as outfile:
        filter_blast(args.blast, crt_dic, outfile)

