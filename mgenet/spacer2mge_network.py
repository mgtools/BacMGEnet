import sys, os
import argparse
import networkx as nx

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", help = "output GML file path")
    parser.add_argument("--gff", help="gff file")
    args = parser.parse_args()

    G = nx.Graph()

    spacer_dic = {}

    with open(args.gff) as infile:
        for line in infile:
            if line[0] == '#':
                continue
            line = line.rstrip('\n')
            sline = line.split('\t')
            sline2 = sline[-1].split(";")
            sseqid = sline[0]
            spacerid = sline2[-1][15:-2]
            if sline[2] == 'contig':
                if "phage" in sline2[1]:
                    G.add_node(sseqid , vtype='phage')
                elif "plasmid" in sline2[1]:
                    G.add_node(sseqid , vtype='plasmid')
            else:
                if spacerid not in spacer_dic:
                    G.add_node(spacerid, vtype='spacer')
                    spacer_dic[spacerid] = 1
                G.add_edge(spacerid, sseqid)

    nx.write_gml(G, args.outfile)
