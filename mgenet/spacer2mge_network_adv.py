import sys, os
import argparse
import networkx as nx
import numpy as np

cluster2spacer_dict = dict()
spacer2cluster_dict = dict()
rep2cluster_dict = dict()
cluster2rep_dict = dict()
spacer2host_dict = dict()
spacer2ctype_dict = dict()
spacer2leader_dict = dict()
spacer2trailer_dict = dict()
host2spacer_dict = dict()
host2cluster_dict = dict()
host2sample_dict = dict()
spacer2mge_dict = dict()
multitype = 0
host_select = []

def get_spacercluster(cdhit_infile):
    '''
    get clustering of the spacers information from cdhit clsr file
    '''
    global cluster2spacer_dict
    global spacer2cluster_dict
    global rep2cluster_dict

    with open(cdhit_infile, "r") as infile:
        for line in infile:
            line = line.rstrip('\n')
            # get cluster numer
            if line.startswith('>Cluster'):
                cluster_num = int(line.lstrip('>Cluster '))+1
            else:
                header = (line.split()[2].strip('...').strip('>')) 
                if cluster_num not in cluster2spacer_dict:
                    cluster2spacer_dict[cluster_num] = list()
                cluster2spacer_dict[cluster_num].append(header)
                spacer2cluster_dict[header]=cluster_num
                if line[-1] == '*':
                    rep2cluster_dict[header] = cluster_num 
                    cluster2rep_dict[cluster_num] = header 

def get_crisprarray(argsori):
    ''' 
    get host information from the crispr array files (.ori files)
    argsori: "ori-file:type ori-file:type.."
    '''
    global spacer2host_dict 
    global spacer2ctype_dict
    global spacer2leader_dict, spacer2trailer_dict
    global host2spacer_dict
    global host2cluster_dict
    global host2sample_dict

    subs = argsori.split()
    for one in subs:
        tmps = one.split(":")
        if len(tmps) > 1:
            ctype, orifile = tmps[0], tmps[1]
        else:
            ctype, orifile = "all", one
        with open(orifile) as ori:
            for line in ori:
                line = line.strip('\n')
                sline = line.split()
                sline2 = sline[0].split(":")
                hostid = sline2[0]
                if len(sline2) > 1:
                    sample = sline2[1]
                else:
                    sample = "unk"
                host2sample_dict[hostid] = sample
                # get cluster information
                #thislist, thisclust = [], []
                if hostid in host2spacer_dict:
                    thislist, thisclust = host2spacer_dict[hostid], host2cluster_dict[hostid]
                else:
                    thislist, thisclust = [], []
                nodes = sline[2:-1]
                add = 0
                for n in nodes:
                    add += 1
                    node, spacer_header = n.strip(']').split('[')
                    spacer2host_dict[spacer_header] = hostid 
                    spacer2ctype_dict[spacer_header] = ctype
                    if sline[1] == 'end':
                        spacer2leader_dict[spacer_header] = add
                    else:
                        spacer2leader_dict[spacer_header] = -1
                    if sline[-1] == 'end':
                        spacer2trailer_dict[spacer_header] = len(nodes) - add + 1
                    else:
                        spacer2trailer_dict[spacer_header] = -1
                    if spacer_header not in thislist:
                        thislist.append(spacer_header)
                        clust = spacer2cluster_dict[spacer_header]
                        if clust not in thisclust:
                            thisclust.append(clust)
                host2spacer_dict[hostid] = thislist
                host2cluster_dict[hostid] = thisclust

def sel_host(contain):
    '''
    select the minimum number of hosts to cover all spacers (with or without protospacers in MGEs)
    using greedy algorithm
    '''
    global host2cluster_dict, rep2cluster_dict, host2sample_dict
    host = list(host2cluster_dict.keys())
    cluster = list(rep2cluster_dict.keys())
    tocover = len(cluster)
    info = []
    for idx in range(len(host)):
        #print(f'samples {host2sample_dict[host[idx]]}')
        if (not contain) or (contain in host2sample_dict[host[idx]]): #if restrict sample/host, set to 0 
            info.append([idx, len(host2cluster_dict[host[idx]])])
        else:
            info.append([idx, 0])
    global host_select
    cluster_explained = [0] * (len(cluster) + 1)
    tocover0 = tocover
    print(f"tocover {tocover}")
    while tocover > 0:
        infosrt = sorted(info, key=lambda item: item[1], reverse=True)
        thishost = infosrt[0][0]
        thisnum = 0
        #once a host is selected, all of its unexplained spacer clusters are consided to be explained
        for cidx in host2cluster_dict[host[thishost]]:
            if cluster_explained[cidx] == 0:
                thisnum += 1
                #the same spacer cluster in all hosts are also considered to be explained
                for allhost in range(len(host)):
                    if cidx in host2cluster_dict[host[allhost]]:
                        info[allhost][1] -= 1
                cluster_explained[cidx] = 1
        #print(f" host {host[thishost]} spacer {len(host2cluster_dict[host[thishost]])} explained {thisnum}, tocover {tocover}")
        if thisnum == 0: 
            if not contain:
                print(f"Warning: cannot find representive host {info}")
                sys.exit()
            else:
                print(f"Warning: not all spacers are considered as hosts are constrained")
                break
        tocover -= thisnum
        host_select.append(host[thishost])
    print(f"total host {len(host)} selected {len(host_select)} total spacers {tocover0} explained {tocover0 - tocover}")

    return (cluster_explained)

def get_ctype(spacer_rep):
    cluster_idx = rep2cluster_dict[spacer_rep]
    spacers = cluster2spacer_dict[cluster_idx]
    ctypes = []
    global multitype
    for spacer in spacers:
        ctype = spacer2ctype_dict[spacer]
        if ctype not in ctypes:
            ctypes.append(ctype)
    if not ctypes:
        print(f"Warning no type found for {spacer} (cluster_idx {cluster_idx} spacers {spacers}")
    elif len(ctypes) == 1:
        return ctypes[0]
    else:
        print(f'ctypes {ctypes}')
        multitype += 1
        ctypes.sort()
        return " ".join(ctypes)
    
def process_spacer2mge0(args):
    G = nx.Graph()
    spacer_dict = {}
    global multitype
    global spacer2mge_dict
    nodelist = []
    edgelist = []
    mgetotal = 0
    spacertotal = 0
    with open(args.gff) as infile:
        for line in infile:
            if line[0] == '#':
                continue
            line = line.rstrip('\n')
            sline = line.split('\t')
            sline2 = sline[-1].split(";")
            sseqid = sline[0]
            if sline[2] == 'contig':
                if "phage" in sline2[1]:
                    G.add_node(sseqid , vtype='phage', ctype='phage')
                    nodelist.append([sseqid, "phage", "phage"])
                elif "plasmid" in sline2[1]:
                    G.add_node(sseqid , vtype='plasmid', ctype='plasmid')
                    nodelist.append([sseqid, "plasmid", "plasmid"])
                mgetotal += 1
            else:
                spacerid0 = sline2[-1][15:-2]
                sline3 = spacerid0.split("', '")
                for spacerid in sline3:
                    if spacerid not in spacer_dict:
                        ctype = get_ctype(spacerid)
                        G.add_node(spacerid, vtype='spacer', ctype=ctype)
                        nodelist.append([spacerid, "spacer", ctype])
                        spacer_dict[spacerid] = 1
                        spacer2mge_dict[spacerid] = sseqid
                        spacertotal += 1
                    G.add_edge(spacerid, sseqid)
                    edgelist.append([spacerid, sseqid])
        infile.close()

    print(f"spacer2mge saved to {args.spacer2mge} mge {mgetotal} spacer {spacertotal} edge {G.number_of_edges()} {multitype} spacers of multiple CRISPR types")
    nx.write_gml(G, args.spacer2mge)

    return (nodelist, edgelist)

def isvalid_spacer(spacer):
    '''
    as long as one of the spacer in the same cluster is from a selected host, the input spacer is valid
    '''
    global host_select, spacer2cluster_dict, cluster2spacer_dict, spacer2host_dict
    c = spacer2cluster_dict[spacer]
    for s in cluster2spacer_dict[c]:
        h = spacer2host_dict[s]
        if h in host_select:
            return True
    return False

def leader_trailer(spacer):
    c = spacer2cluster_dict[spacer]
    leader, trailer = [], []
    for s in cluster2spacer_dict[c]:
        if spacer2leader_dict[s] != -1:
            leader.append(spacer2leader_dict[s])
        if spacer2trailer_dict[s] != -1:
            trailer.append(spacer2trailer_dict[s])
    leaderdis = trailerdis = 100
    if leader:
        leaderdis = np.median(leader)
    if trailer:
        trailerdis = np.median(trailer)

    if (leaderdis < trailerdis) and (leaderdis <= 2):
        return "leader"
    elif (trailerdis < leaderdis) and (trailerdis <= 2):
        return "trailer"
    return "middle"

def process_spacer2mge(args):
    G = nx.Graph()
    global multitype
    global spacer2mge_dict
    nodelist = []
    edgelist = []
    mgetotal = 0
    spacertotal = 0
    spacer_dict = {}
    nodeincluded = []
    with open(args.gff) as infile:
        for line in infile:
            if line[0] == '#':
                continue
            line = line.rstrip('\n')
            sline = line.split('\t')
            sline2 = sline[-1].split(";")
            sseqid = sline[0]
            if sline[2] == 'contig':
                if "phage" in sline2[1]:
                    nodelist.append([sseqid, "phage", "phage"])
                elif "plasmid" in sline2[1]:
                    nodelist.append([sseqid, "plasmid", "plasmid"])
                mgetotal += 1
            else:
                spacerid0 = sline2[-1][15:-2]
                sline3 = spacerid0.split("', '")
                for spacerid in sline3:
                    if not isvalid_spacer(spacerid):
                        continue
                    if spacerid not in spacer_dict:
                        ctype = get_ctype(spacerid)
                        nodelist.append([spacerid, "spacer", ctype])
                        spacer_dict[spacerid] = 1
                        spacer2mge_dict[spacerid] = sseqid
                        spacertotal += 1
                    edgelist.append([spacerid, sseqid])
                    if sseqid not in nodeincluded:
                        nodeincluded.append(sseqid)
        infile.close()

    nodelist2 = []
    mgetotal2 = 0
    spacertotal2 = 0
    for node in nodelist:
        #some MGEs are excluded if all protospacers were excluded
        if (node[1] == 'spacer') or (node[0] in nodeincluded):
            nodelist2.append(node)
    for node in nodelist2:
        if node[1] == 'spacer':
            G.add_node(node[0], vtype=node[1], ctype=node[2], dtype=leader_trailer(node[0]))
            spacertotal2 += 1
        else:
            G.add_node(node[0], vtype=node[1], ctype=node[2], dtype=node[1])
            mgetotal2 += 1
    for edge in edgelist:
        G.add_edge(edge[0], edge[1])

    print(f"spacer2mge saved to {args.spacer2mge} mge {mgetotal} selected {mgetotal2} spacer {spacertotal} included {spacertotal2} edge {G.number_of_edges()} {multitype} spacers of multiple CRISPR types")
    nx.write_gml(G, args.spacer2mge)

    return (nodelist2, edgelist)

def process_host2mge(args, nodelist, edgelist):
    global host_select
    G = nx.Graph()
    spacerwithmge = {}
    mgetotal = 0
    hosttotal = 0
    spacertotal = 0
    mge2node_dict = {}
    mgeincluded = []
    for idx in range(len(nodelist)):
        node = nodelist[idx]
        if node[1] != 'spacer':
            mgetotal += 1
            mge2node_dict[node[0]] = idx 
            #G.add_node(node[0], vtype=node[1], ctype=node[2])
            #mgeincluded.append(node[0])
        else:
            spacerwithmge[node[0]] = 1
    for host in host_select:
        add = 0
        for cidx in host2cluster_dict[host]:
            for spacer in cluster2spacer_dict[cidx]:
                if spacer in spacerwithmge:
                    mge = spacer2mge_dict[spacer]
                    if mge not in mgeincluded:
                        mgeincluded.append(mge)
                        node = nodelist[mge2node_dict[mge]]
                        G.add_node(node[0], vtype=node[1], ctype=node[2])
                    G.add_edge(spacer2mge_dict[spacer], host) 
                    add += 1
                    spacertotal += 1
        if add > 0:
            if args.contain:
                ctype = args.contain
            else:
                ctype = "host"
            G.add_node(host, vtype="host", ctype=ctype)
            hosttotal += 1
    nx.write_gml(G, args.host2mge)
    print(f"host2mge network saved to {args.host2mge} mge included {len(mgeincluded)} host included {hosttotal} spacer {spacertotal} edge {G.number_of_edges()} spacerwithmge {len(spacerwithmge)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", help="gff file of the MGEs (containing protospacer information)", required=True)
    parser.add_argument("--cdhit", help="clusting of the spacers by cdhit", required=True)
    parser.add_argument("--ori", help = "module ori file, of crispr arrays in the format of tag:file, e.g., L47:BfragL47_module.all.ori; multiple files seperated by space", required=True)
    parser.add_argument("--spacer2mge", help = "output GML file of the spacer-MGE network", required=True)
    parser.add_argument("--host2mge", help = "output GML file of the host-MGE network", required=True)
    parser.add_argument("--contain", help = "restrict host (e.g., S01)")
    args = parser.parse_args()

    get_spacercluster(args.cdhit) 
    get_crisprarray(args.ori)
    cluster_explained = sel_host(args.contain)
    (nodelist, edgelist) = process_spacer2mge(args)
    process_host2mge(args, nodelist, edgelist)

