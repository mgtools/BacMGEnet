import sys, os
import argparse
import networkx as nx

def build_cluster2header_dict(cdhit_infile):

    cluster2header_dict = dict()
    header2cluster_dict = dict()

    with open(cdhit_infile) as infile:
        for line in infile:
            line = line.rstrip('\n')
            # get cluster numer
            if line.startswith('>Cluster'):
                cluster_num = int(line.lstrip('>Cluster '))+1
            else:
                header = (line.split()[2].strip('...').strip('>')) 
                if cluster_num not in cluster2header_dict:
                    cluster2header_dict[cluster_num] = list()
                cluster2header_dict[cluster_num].append(header)
                header2cluster_dict[header]=cluster_num

    return(cluster2header_dict, header2cluster_dict)

def get_cluster2node(argsori, header2cluster_dict):
    ''' 
    returns a clsuter2node dict
    '''
    cluster2node_dict = dict()

    node2header_dict = dict()
    header2node_dict = dict()

    with open(argsori) as ori:
        for line in ori:
            line = line.strip('\n')
            sline = line.split()

            # get cluster information
            nodes = sline[2:-1]
            for n in nodes:
                node, spacer_header = n.strip(']').split('[')
    
                # add cluster number if not available
                if node not in node2header_dict:
                    node2header_dict[node] = list()

                node2header_dict[node].append(spacer_header)
                header2node_dict[spacer_header] = node

    for header in header2cluster_dict:
        cluster = str(header2cluster_dict[header])
        node = str(header2node_dict[header])

        if cluster not in cluster2node_dict:
            cluster2node_dict[cluster] = node
        else:
            if node == cluster2node_dict[cluster]:
                continue
            else:
                print('node with different cluster annotations detected, please debug me')

    return(cluster2node_dict, header2node_dict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", help = "output GML file path")
    parser.add_argument("--phage", help="phage blast results")
    parser.add_argument("--plasmid", help = "plasmid blast results")
    parser.add_argument("-c", "--cdhit", help="cdhit output")
    parser.add_argument("--ori", help = "module ori file, crispr spacer graph")
    args = parser.parse_args()

    G = nx.Graph()

    #cluster2header_dict, header2cluster_dict = build_cluster2header_dict(args.cdhit) 
    #cluster2node_dict, header2node_dict = get_cluster2node(args.ori, header2cluster_dict)

    #for cluster in cluster2header_dict:
    #    G.add_node(f'Node {cluster}', vtype='spacer')

    with open(args.phage) as phage:
        for line in phage:
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

            G.add_node(S_number, vtype='Host')
            G.add_node(sseqid, vtype='phage')
            G.add_edge(S_number, sseqid)

    with open(args.plasmid) as plasmid:
        for line in plasmid:
            line = line.rstrip('\n')
            sline = line.split('\t')
            qseqid = sline[0]
            sseqid = sline[1]
            pident = sline[2]
            qstart = sline[7]
            qend   = sline[8]
            protospacer_start = sline[9]
            protospacer_end   = sline[10] 

            sseqid = sseqid.split('|')[0]

            S_number = qseqid.split(':')[1].split('-')[0]

            G.add_node(S_number, vtype='Host')
            G.add_node(sseqid, vtype='plasmid')
            G.add_edge(S_number, sseqid)


    nx.write_gml(G, args.outfile)
