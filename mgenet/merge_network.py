import sys, os
import argparse
import networkx as nx

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", help = "output GML file path")
    parser.add_argument("--gml", help="gml files (at least two) seperated by a space", required=True)
    args = parser.parse_args()

    gmls = args.gml.split()
    if len(gmls) < 2:
        sys.exit("need at least gml files")
    G = nx.read_gml(gmls[0])

    for gml in gmls[1:]:
        H = nx.read_gml(gml)
        G = nx.compose(G,H)

    nx.write_gml(G, args.outfile)
