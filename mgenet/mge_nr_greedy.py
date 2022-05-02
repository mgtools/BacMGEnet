import sys, os
import pandas as pd
import argparse
import Bio.SeqIO as IO

def get_protospacer_dic(gff_infile, gff_outfile):
    '''
    mge_dic[query_id] = [idx]
    '''
    
    mge_dic = dict()
    protospacer_dic = dict()
    mge, spacer = 0, 0
    mge_info, mge_to_protospacer = [], []
    protospacer_info, protospacer_to_mge = [], []
    mge_lines = []
    with open(gff_infile) as infile:
        for line in infile:
            if line[0] == '#':
                continue
            line = line.rstrip('\n')
            sline = line.split('\t')
            qseqid = sline[0]
            sline2 = sline[-1].split(';')
            if sline[2] == 'contig':
                slen, mgetype, nprotospacer = int(sline[4]), sline2[1][5:], int(sline2[2][15:])
                mge_dic[qseqid] = mge
                mge_info.append([qseqid, slen, mgetype, nprotospacer])
                mge_to_protospacer.append([])
                mge_lines.append([line,])
                mge += 1
            else:
                sbeg, send, spacerid = int(sline[3]), int(sline[4]), sline2[2][14:-2]
                if spacerid not in protospacer_dic:
                    protospacer_dic[spacerid] = spacer
                    sidx = spacer
                    protospacer_info.append([spacerid, 1, sbeg, send])
                    protospacer_to_mge.append([mge - 1, ])
                    spacer += 1 
                else:
                    sidx = protospacer_dic[spacerid]
                    if (mge - 1) not in protospacer_to_mge[sidx]:
                        protospacer_to_mge[sidx].append(mge - 1)
                if sidx not in mge_to_protospacer[mge - 1]:
                    mge_to_protospacer[mge - 1].append(sidx)
                mge_lines[mge - 1].append(line)
        mge_to_protospacer.append(spacer)
        infile.close()

    tocover = spacer
    info = []
    for idx in range(mge):
        info.append([idx, len(mge_to_protospacer[idx])])
    mge_select = []
    mge_explained = [0] * mge
    while tocover > 0:
        infosrt = sorted(info, key=lambda item: item[1], reverse=True)
        thismge = infosrt[0][0]
        thisnum = 0
        toexplain = []
        for sidx in mge_to_protospacer[thismge]:
            if protospacer_info[sidx][1] == 1:
                toexplain.append(sidx)
        #print(f"tocover {tocover} best mge {thismge} {mge_info[thismge][0]}, n-protospacer {mge_info[thismge][3]} to-explain {infosrt[0][1]}, all {mge_to_protospacer[thismge]}, to explain {toexplain}")
        if len(toexplain) == 0:
            break 
        for sidx in mge_to_protospacer[thismge]:
            if protospacer_info[sidx][1] == 1:
                thisnum += 1
                for allmge in protospacer_to_mge[sidx]:
                    info[allmge][1] -= 1
                protospacer_info[sidx][1] = 0 
        #print(f" --- explained {thisnum}")
        tocover -= thisnum
        mge_select.append(thismge)
        mge_explained[thismge] = thisnum

    phage, phage_sel, plasmid, plasmid_sel = 0, 0, 0, 0
    for thismge in range(mge):
        if mge_info[thismge][2] == 'phage':
            phage += 1
            if mge_explained[thismge]:
                phage_sel += 1
        elif mge_info[thismge][2] == 'plasmid':
            plasmid += 1
            if mge_explained[thismge]:
                plasmid_sel += 1

    outfile = open(gff_outfile, "w")
    outfile.write("##gff-version 3\n")
    for thismge in mge_select:
        outfile.write(f"{mge_lines[thismge][0]};explained={mge_explained[thismge]}\n")
        for line in mge_lines[thismge][1:]:
            outfile.write(line + "\n")
    outfile.close()

    print(f"input {gff_infile}: spacer {spacer} mge {mge} (phage {phage} plasmid {plasmid}) selected mge {len(mge_select)} (phage {phage_sel} plasmid {plasmid_sel})")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", help = "gff file as input", default = 'BfragLall_protospacer_w_multihits.gff')
    parser.add_argument("--out", help = "new gff file (nr)", default = 'BfragLall_protospacer_w_multihits-nr.gff')
    args = parser.parse_args()

    get_protospacer_dic(args.gff, args.out)
