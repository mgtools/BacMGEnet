import sys, os
import argparse

class Spacer:
    '''
    build spacer metadata information
    '''

    def __init__(self, spacer):
        self.spacer = spacer
        self.individual = None
        self.cluster = None
        self.mge_hits = list()

    def add_individual(self, individual):
        self.individual = individual

    def add_cluster(self, cluster):
        self.cluster = cluster

    def add_hits(self, mge_header):
        self.mge_hits.append(mge_header)

    def which_spacer(self):
        return(self.spacer)
    
    def which_individual(self):
        return(self.individual)

    def which_cluster(self):
        return(self.cluster)

    def what_hits(self):
        return(self.mge_hits)

def get_blast_dict(blast_input, spacer_list):
    ''' 
    returns blast dictionaries
    dict[query_header] = list of subject
    dict[subject_header] = list of query
    '''
    
    with open(blast_input) as blast:
        for line in blast:
            line = line.rstrip('\n')
            sline = line.split('\t')
            qseqid = sline[0] # spacer ID
            sseqid = sline[1] # MGE id
            pident = sline[2]
            qstart = sline[7]
            qend   = sline[8]
            protospacer_start = sline[9]
            protospacer_end   = sline[10]

            sseqid = sseqid.split('|')[0]

            S_number = qseqid.split(':')[1].split('-')[0]

            spacer_list.append(Spacer(qseqid))            

    return(spacer_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", help = "output GML file path")
    parser.add_argument("--phage", help="phage blast results")
    parser.add_argument("--plasmid", help = "plasmid blast results")
    parser.add_argument("-c", "--cdhit", help="cdhit output")
    parser.add_argument("--ori", help = "module ori file, crispr spacer graph")
    args = parser.parse_args()

    spacer_list = list()

    spacer_list = get_blast_dict(args.phage, spacer_list)

    print(spacer_list)

    print([ i.which_spacer() for i in spacer_list ])
    print(SRR9713593:S06-0001-1:NODE_15_length_137879_cov_13.360278:c1:p106396.which_spacer())
