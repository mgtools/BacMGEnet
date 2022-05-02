import sys, os
import argparse
import networkx as nx
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", help = "output GML file path")
    parser.add_argument("--gml", help="gml files (at least two) seperated by a space", required=True)
    parser.add_argument("--bac", help="annotation of the host", required=False)
    parser.add_argument("--phage", help="annotation of the host", required=False)
    args = parser.parse_args()

    G = nx.read_gml(args.gml)
    host_nodes = [n for n,v in G.nodes(data=True) if v['vtype'] == 'host']
    phage_nodes = [n for n,v in G.nodes(data=True) if v['vtype'] == 'phage']
    plasmid_nodes = [n for n,v in G.nodes(data=True) if v['vtype'] == 'plasmid']
    spacer_nodes = [n for n,v in G.nodes(data=True) if v['vtype'] == 'spacer']
    print(f"total nodes: {G.number_of_nodes()}")
    print(f"host-nodes: {len(host_nodes)}")
    print(f"spacer-nodes: {len(spacer_nodes)}")
    print(f"phage-nodes: {len(phage_nodes)}")
    print(f"plasmid-nodes: {len(plasmid_nodes)}")
    print(f"  phage and plasmids: {len(phage_nodes) + len(plasmid_nodes)}")

    comp = sorted(nx.connected_components(G), key=len, reverse=True)
    compsize = [len(c) for c in comp]
    idx = 0
    for one in comp:
        idx += 1
        print(f"comp {idx} size {compsize[idx - 1]} {one}")

    if args.bac:
        inf = open(args.bac, "r")
        taxon = {}
        for aline in inf:
            sub = aline.strip().split("\t")
            tmp = sub[-1]
            sub2 = tmp.split(";")
            tlist = []
            for one in sub2:
                tlist.append(one[3:])
            if len(tlist) < 2:
                print(f"wrong format {aline}")
            taxon[sub[0]] = tlist 
        inf.close()

        idx = 0 
        for c in comp:
            idx += 1
            print(f"Component {idx} size {compsize[idx - 1]}")
            names = list(c)
            phylum, spe = {}, {}
            bac = 0
            for one in names:
                if one not in taxon:
                    continue
                bac += 1
                p, s = taxon[one][1], taxon[one][-1]
                if p not in phylum:
                    phylum[p] = 1
                else:
                    phylum[p] = phylum[p] + 1
                if s:
                    if s not in spe:
                        spe[s] = 1
                    else:
                        spe[s] = spe[s] + 1
            print(f"  total bac {bac} phylum {len(phylum.keys())}")
            for p in phylum.keys():
                print(f"     phylum {p}  {phylum[p]}")
            print(f"  total species {len(spe.keys())}")
            for s in spe.keys():
                print(f"     species {s}  {spe[s]}")

