import sys, os
import argparse

def build_derep_list(Wdb):
    '''
    returns a list of derep genomes; representatives
    '''
    
    rep_list = list()

    with open(Wdb) as infile:
        for i, line in enumerate(infile):
            if i == 0:
                continue
            sline = line.split('.fna')
            rep_list.append(sline[0])

    return(rep_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--blast", help="blast results")
    parser.add_argument("--derep", help = "derep output Wdb file")
    args = parser.parse_args()

    rep_list = build_derep_list(args.derep)

    with open(args.blast) as blast_out:
        with open(os.path.splitext(args.blast)[0]+'.derep.m8', 'w') as outfile:
            for line in blast_out:
                nline = line.rstrip('\n')
                sline = nline.split('\t')
                header = (sline[1].split('|')[0])    

                if header in rep_list:
                    outfile.write(line)
