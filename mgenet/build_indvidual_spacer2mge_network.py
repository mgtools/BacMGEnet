import sys, os
import argparse
import networkx as nx

def individual2spacerclusters(reference_arrays):

    individual2spacercluster_dict = dict()

    with open(reference_arrays) as reference:
        for line in reference:
            line = line.rstrip('\n').split()
            S_number = line[0].split(':')[1].split('-')[0]

            #check for key
            if S_number not in individual2spacercluster_dict:
                individual2spacercluster_dict[S_number] = list()

            # add cluster numbers
            individual2spacercluster_dict[S_number] += line[2:-1]
            individual2spacercluster_dict[S_number] = list(set(individual2spacercluster_dict[S_number]))

    return(individual2spacercluster_dict)

def build_cluster2header_dict(cdhit_infile):

    cluster2header_dict = dict()
    header2cluster_dict = dict()

    with open(cdhit_infile) as infile:
        for line in infile:
            line = line.rstrip('\n')
            # get cluster numer
            if line.startswith('>Cluster'):
                cluster_num = int(line.lstrip('>Cluster '))+1
                cluster_num = str(cluster_num)
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

def get_blast_dict(blast_input, header2cluster_dict, cluster2mge):
    '''
    updates blast dict to include results
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

            if qseqid not in header2cluster_dict:
                print(f'{qseqid} not in header2cluster_dict')

            if header2cluster_dict[qseqid] not in cluster2mge:
                cluster2mge[header2cluster_dict[qseqid]] = list()
            cluster2mge[header2cluster_dict[qseqid]].append(sseqid)
            cluster2mge[header2cluster_dict[qseqid]] = list(set(cluster2mge[header2cluster_dict[qseqid]]))

    return(cluster2mge)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference_arrays", help = "reference array for each individual")
    parser.add_argument("-o", "--outdir", help = "output GML file path")
    parser.add_argument("--phage", help="phage blast results")
    parser.add_argument("--plasmid", help = "plasmid blast results")
    parser.add_argument("-c", "--cdhit", help="cdhit output")
    parser.add_argument("--ori", help = "module ori file, crispr spacer graph")
    args = parser.parse_args()

    individual2spacercluster_dict = individual2spacerclusters(args.reference_arrays)
    cluster2header_dict, header2cluster_dict = build_cluster2header_dict(args.cdhit) 
    cluster2node_dict, header2node_dict = get_cluster2node(args.ori, header2cluster_dict)


    # initiate dictionary
    phage_cluster2mge = dict()
    plasmid_cluster2mge = dict()

    phage_cluster2mge = get_blast_dict(args.phage, header2cluster_dict, phage_cluster2mge)
    plasmid_cluster2mge = get_blast_dict(args.plasmid, header2cluster_dict, plasmid_cluster2mge)

    #print(phage_cluster2mge)
    #print(plasmid_cluster2mge)

    for S in individual2spacercluster_dict:

        # start network for individual
        G = nx.Graph()

        for cluster in individual2spacercluster_dict[S]:
            G.add_node(f"Spacer {cluster2node_dict[cluster]}", vtype='spacer') # add spacer node

            if cluster in phage_cluster2mge:
                for phage in phage_cluster2mge[cluster]:
                    G.add_node(phage, vtype='phage') 
                    G.add_edge(f"Spacer {cluster2node_dict[cluster]}", phage)
                    
    
            if cluster in plasmid_cluster2mge:
                for plasmid in plasmid_cluster2mge[cluster]:
                    G.add_node(plasmid, vtype='plasmid')
                    G.add_edge(f"Spacer {cluster2node_dict[cluster]}", plasmid)

        nx.write_gml(G, f'{args.outdir}/{S}_L47_reference_spacer2mge.gml')
