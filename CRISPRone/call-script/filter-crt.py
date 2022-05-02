#!/usr/bin/env python3

import sys
import re
if len(sys.argv) < 4:
    sys.exit(sys.argv[0] + " old-crt-file gff-file new-crt-file")
oldcrtfile, gfffile, newcrtfile = sys.argv[1:4]
crisprseq, crisprbeg, crisprend, crisprinfo, crispr2seq = [], [], [], [], [] 
seqlist = []
seqlen = []
seqcrispr = [] 
inf = open(gfffile, "r")
for aline in inf:
    if "direct_repeat" in aline:
        subs = aline.strip().split()
        crisprseq.append(subs[0])
        crisprbeg.append(int(subs[3]))
        crisprend.append(int(subs[4]))
        crisprinfo.append("")
        if subs[0] not in seqlist:
            seqlist.append(subs[0])
            seqcrispr.append(1)
            seqlen.append(0)
        else:
            seqcrispr[seqlist.index(subs[0])] += 1
inf.close()

inf = open(oldcrtfile, "r")
sidx = thiscidx = -1
valid = False
for aline in inf:
    aline = aline.strip()
    m1 = re.match("^SEQ: (\S+)", aline)
    m2 = re.match("CRISPR (\d+)\s+Range: (\S+) - (\S+)", aline)
    m3 = re.match("^Bases: (\d+)", aline)
    if m1:
        currseq = m1.groups()[0]
        print(f"currseq found {currseq}")
        if currseq in seqlist:
            sidx = seqlist.index(currseq)
            m1b = re.match("^SEQ: (\S+) length (\d+)", aline)
            if m1b:
                seqlen[seqlist.index(currseq)] = m1b.groups()[1]
            valid = True
        else:
            sidx = -1
            valid = False
    elif m3:
        if currseq in seqlist:
            seqlen[seqlist.index(currseq)] = m3.groups()[0]
    elif m2:
        print(f"range found {m2.groups()[1]} {m2.groups()[2]}")
        thisbeg, thisend = int(m2.groups()[1]), int(m2.groups()[2])
        posmatch = False
        for cidx in range(len(crisprseq)):
            if (crisprseq[cidx] == currseq) and (thisbeg >= crisprbeg[cidx] - 40) and (thisend <= crisprend[cidx] + 40): #not exact coordinate match for partial repeats
                posmatch = True
                thiscidx = cidx
                break
        if (sidx != -1) and posmatch:
            print(f"range matched {m2.groups()[1]} {m2.groups()[2]}")
            valid = True
        else:
            valid = False
    elif valid and aline:
        m4 = re.match("^Bases: ", aline)
        m5 = re.match("^Total ", aline)
        if not (m4 or m5):
            crisprinfo[thiscidx] += aline + "\n"
inf.close()

out = open(newcrtfile, "w")
for sidx in range(len(seqlist)):
    out.write(f"SEQ: {seqlist[sidx]} length {seqlen[sidx]}\n")
    out.write(f"Bases: {seqlen[sidx]}\n\n")
    clustidx = 0
    for cidx in range(len(crisprseq)):
        if crisprseq[cidx] == seqlist[sidx]:
            clustidx += 1
            out.write(f"CRISPR {clustidx} Range: {crisprbeg[cidx]} - {crisprend[cidx]}\n")
            out.write(crisprinfo[cidx] + "\n")
    
out.close()
