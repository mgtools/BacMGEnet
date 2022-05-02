#!/usr/bin/env python

import sys
fontsize = 18

class Graph:
        def __init__(self, nodelist, adjmatrix):
                self.nodelist = nodelist
                self.adjmatrix = adjmatrix
                self.num = len(self.nodelist)
                self.supernode = []

        #collapse graph -- collapse one-in-one-out nodes (spacers)
        def nocollapse(self):
            self.supernode = []
            for idx in range(self.num):
                self.supernode.append([idx,])
            self.supernodegraph()
            print(f"total supernodes {len(self.supernode)}")

        def collapse(self):
                self.supernode = []
                notvisited = list(range(self.num)) #notvisited -- store node indices not actual names 
                print(f"notvisited {notvisited}")
                while(notvisited):
                        tovisit = notvisited.pop(0)
                        self.supernode.append([tovisit,])
                        #forward
                        outnodes = self.getOutNode(tovisit)
                        while(len(outnodes) == 1):
                                if outnodes[0] not in notvisited:
                                        break
                                innodes = self.getInNode(outnodes[0])
                                if len(innodes) > 1: 
                                        break
                                if int(self.nodelist[outnodes[0]]) != int(self.nodelist[self.supernode[-1][-1]]) + 1: #discontinue in numbering
                                        break
                                self.supernode[-1].append(outnodes[0])
                                del notvisited[notvisited.index(outnodes[0])]
                                outnodes = self.getOutNode(outnodes[0])
                                
                        #backward       
                        innodes = self.getInNode(tovisit)
                        while(len(innodes) == 1):
                                if innodes[0] not in notvisited:
                                        break
                                outnodes = self.getOutNode(innodes[0])
                                if len(outnodes) > 1:
                                        break
                                if int(self.nodelist[innodes[0]]) != int(self.nodelist[self.supernode[-1][0]]) - 1: #discontinue in numbering
                                        break
                                self.supernode[-1].insert(0, innodes[0])
                                del notvisited[notvisited.index(innodes[0])]
                                innodes = self.getInNode(innodes[0])
                        print(f"supernode {self.supernode[-1]}")
                self.supernodegraph()
                print(f"total supernodes {len(self.supernode)}")

        def getpath(self, spacerlist):
                print("getpath..")
                spidxlist = [-1] * len(spacerlist)
                for sidx in range(len(spacerlist)):
                        spacer = spacerlist[sidx]
                        spidx = -1
                        for s1idx in range(len(self.supernode)):
                                s1 = self.supernode[s1idx]
                                s1n = self.supernodename(s1idx)
                                tmp = s1n.split("-")
                                if (len(tmp) == 1) and (int(spacer) == int(tmp[0])):
                                        spidx = s1idx
                                        break
                                elif (len(tmp) > 1) and (int(spacer) <= int(tmp[1])) and (int(spacer) >= int(tmp[0])):
                                        spidx = s1idx
                                        break
                        if spidx == -1:
                                print(f"spacer {spacer} not found in spacergraph with supernodes")
                                sys.exit("wrong spacergraph")   
                        spidxlist[sidx] = spidx
                self.onpath = []
                for idx in range(len(self.supernode)):
                        self.onpath.append([0] * len(self.supernode))
                for idx in range(len(spacerlist) - 1):
                        self.onpath[spidxlist[idx]][spidxlist[idx + 1]] = 1

        def supernodelen(self, n):
                i = int(self.nodelist[self.supernode[n][0]])
                j = int(self.nodelist[self.supernode[n][-1]])
                return j - i + 1

        def supernodename(self, n):
                i = self.nodelist[self.supernode[n][0]]
                j = self.nodelist[self.supernode[n][-1]]
                #tmp = "g" + str(n) + "_" + str(len(self.supernode[n])) + "_" + i       
                #tmp = "b" + str(n) + "_" + i   
                tmp = i         
                if j != i:
                        tmp += "-" + j  
                return tmp

        def supernodegraph(self):
                self.supernode_edge = []
                self.supernode_edgeweighted = []
                self.supernode_edgeweighted_all = []
                for s1idx in range(len(self.supernode)):
                        self.supernode_edge.append([0] * len(self.supernode))
                        self.supernode_edgeweighted.append([0] * len(self.supernode))
                        self.supernode_edgeweighted_all.append([0] * len(self.supernode))
                        s1 = self.supernode[s1idx]
                        s1n = self.supernodename(s1idx)
                        for s2idx in range(len(self.supernode)):
                                s2 = self.supernode[s2idx]
                                s2n = self.supernodename(s2idx) 
                                edge = 0 
                                edgeweighted = 0
                                edgeweighted_all = 0
                                for n1 in s1:
                                        for n2 in s2:
                                                if self.adjmatrix[n1][n2] > 0:
                                                        #self edge
                                                        edgeweighted_all += self.adjmatrix[n1][n2]
                                                        if s1idx != s2idx:
                                                                edge += 1
                                                                edgeweighted += self.adjmatrix[n1][n2]
                
                                self.supernode_edge[s1idx][s2idx] = edge
                                self.supernode_edgeweighted[s1idx][s2idx] = edgeweighted
                                self.supernode_edgeweighted_all[s1idx][s2idx] = edgeweighted_all

        def writedot(self, outfile, spacers_detail=[], colorfile = "", colordot=""):
                out = open(outfile, "w")
                out.write("digraph SpacerGraph {\n")
                out.write('rankdir="LR";\n')
                out.write(f'node [fontsize = {fontsize}];\n')
                if colorfile: #color the nodes according to the specifications
                    inf = open(colorfile, "r")
                    colorlist, taglist = [], []
                    for aline in inf:
                        subs = aline.split()
                        colorlist.append(subs[0])
                        taglist.append(subs[1])
                    inf.close()
                    saved = ""
                    print(f"colorlist {colorlist} taglist {taglist}")
                    colorused = [False] * len(colorlist)
                    for s1idx in range(len(self.supernode)):
                        s1 = self.supernode[s1idx]
                        s1n = self.supernodename(s1idx)
                        print(f"supernode {s1idx}, nodes {s1}")
                        thiscolor, thiscoloridx = [], []
                        for nidx in s1: 
                            print(f"spacer {nidx} {spacers_detail[nidx]}")
                            for cidx in range(len(colorlist)):
                                if taglist[cidx] in spacers_detail[nidx]:
                                    if colorlist[cidx] not in thiscolor:
                                        thiscolor.append(colorlist[cidx])
                                        thiscoloridx.append(cidx)
                        print(f"Supernode {s1idx} thiscolor {thiscolor}")
                        if len(thiscolor) == 1:
                            color = thiscolor[0]
                            saved += '"' + s1n + '"' + "[style=filled, fillcolor=" + color + "]\n"
                            colorused[thiscoloridx[0]] = True
                        else:
                            saved += '"' + s1n + '"' + '\n'
                    out.write(" subgraph cluster_0 {\n")
                    out.write("  style=filled;\n")
                    out.write("  color=lightgray;\n")
                    schemetag = "Color scheme"
                    out.write(f'  "a"[style=filled, color=lightgray, fillcolor=lightgray, fontcolor=lightgray, fontsize={fontsize + 10}]\n')
                    out.write(f'  "Color scheme"[style=filled, color=lightgray, fillcolor=lightgray, fontsize={fontsize + 10}]\n')
                    tmp = f'"a" -> "Color scheme"'
                    for cidx in range(len(colorlist)):
                        if colorused[cidx]:
                            out.write(f'  "{taglist[cidx]}"[style=filled, fillcolor={colorlist[cidx]}, fontsize={fontsize + 10}]\n')
                            tmp += f' -> "{taglist[cidx]}"'
                    tmp += f' -> "shared"'
                    out.write(f'  "shared"[style=filled, fillcolor=white, fontsize={fontsize + 10}]\n')
                    out.write(f'{tmp} [color=lightgray]\n')
                    out.write('}\n')
                    out.write(saved)
                else:
                    for s1idx in range(len(self.supernode)):
                        s1 = self.supernode[s1idx]
                        s1n = self.supernodename(s1idx)
                        shownode = True
                        color=""
                        if self.supernodelen(s1idx) > 0:
                                shownode = True
                                innode = self.getInSupernode(s1idx)
                                outnode = self.getOutSupernode(s1idx)
                                if (innode) and (not outnode):
                                        color = "yellow"
                                elif (not innode) and (outnode):
                                        color = "skyblue"
                        if shownode:
                                if color:
                                        out.write('"' + s1n + '"' + "[style=filled, fillcolor=" + color + "]\n")
                                else:
                                        out.write('"' + s1n + '"' + '\n')

                #assign weights to edges
                edgeweightlist = []
                for s1idx in range(len(self.supernode)):
                        for s2idx in range(len(self.supernode)):
                                if self.supernode_edge[s1idx][s2idx] > 0:
                                        edgeweightlist.append(self.supernode_edgeweighted[s1idx][s2idx])
                edgeweightlist.sort()
                totedge = len(edgeweightlist)
                if totedge >= 10:
                        quarter = int(totedge / 4)
                        low, medium, high = edgeweightlist[quarter], edgeweightlist[2 * quarter], edgeweightlist[3 * quarter]
                else:
                        low, medium, high = 1, 1, 1
                                
                for s1idx in range(len(self.supernode)):
                        s1 = self.supernode[s1idx]
                        s1n = self.supernodename(s1idx)
                        for s2idx in range(len(self.supernode)):
                                #allow self edge for supernode with only one spacer
                                s2 = self.supernode[s2idx]
                                s2n = self.supernodename(s2idx) 
                                if (s1idx == s2idx) and (len(s1) > 1):
                                        continue
                                edge = self.supernode_edge[s1idx][s2idx] 
                                edgeweighted = self.supernode_edgeweighted_all[s1idx][s2idx]  #use edgeweighted_all
                                if edgeweighted > 0:
                                        penweight = 1
                                        if totedge >= 10:
                                                if edgeweighted > high:
                                                        penweight = 4   
                                                elif edgeweighted > medium:
                                                        penweight = 3
                                                elif edgeweighted > low:
                                                        penweight = 2
                                        if colorfile:
                                                out.write(f'"{s1n}" -> "{s2n}" [penwidth=1]\n')
                                        elif self.onpath[s1idx][s2idx]:
                                                out.write(f'"{s1n}" -> "{s2n}" [color=red, penwidth={penweight}]\n')
                                        else:
                                                out.write(f'"{s1n}" -> "{s2n}" [penwidth={penweight}]\n')
                out.write("}\n")
                out.close()
                print("graph saved to", outfile)

        def getInNode(self, n):
                a = []
                for r in range(self.num):       
                        if self.adjmatrix[r][n] > 0:
                                a.append(r)             
                return a

        def getInSupernode(self, n):
                a = []
                for r in range(len(self.supernode)):    
                        if self.supernode_edge[r][n] > 0:
                                a.append(r)             
                return a

        def getOutNode(self, n):
                a = []
                for c in range(self.num):       
                        if self.adjmatrix[n][c] > 0:
                                a.append(c)
                return a

        def getOutSupernode(self, n):
                a = []
                for c in range(len(self.supernode)):    
                        if self.supernode_edge[n][c] > 0:
                                a.append(c)
                return a

def getspacer(one0):
       if one0[-1] == ']':
                subs = one0[:-1].split("[")
                one, des = subs[0], subs[1]
       else:
                one, des = one0, ""
       return one, des

def main():
        inpfile, dotfile, colorfile, collapse = "", "", "", True
        global fontsize
        for idx in range(len(sys.argv)):
            if (sys.argv[idx] == '-i') and (len(sys.argv) > idx + 1):
                inpfile = sys.argv[idx + 1]
            elif (sys.argv[idx] == '-o') and (len(sys.argv) > idx + 1):
                dotfile = sys.argv[idx + 1]
            elif (sys.argv[idx] == '-c') and (len(sys.argv) > idx + 1):
                colorfile = sys.argv[idx + 1]
            elif (sys.argv[idx] == '-f') and (len(sys.argv) > idx + 1):
                fontsize = int(sys.argv[idx + 1])
            elif sys.argv[idx] == "--no-collapse":
                collapse = False
        if not (inpfile and dotfile):
                sys.exit(sys.argv[0] + " -i inputfile -o outfile <--no-collapse> <-c color.txt> <-d color-dot-file>")
        inf = open(inpfile, "r")
        lines = inf.readlines()
        inf.close()
        seqs = []
        spacers = []
        spacers_detail = []
        for aline in lines:
                subs = aline.strip().split()    
                if len(subs) < 2:
                        continue
                seqs.append(subs[0])
                for one0 in subs[2:-1]:
                        one, des = getspacer(one0)
                        if one not in spacers:
                                spacers.append(one)
                                spacers_detail.append(des)
                        else:
                                spacers_detail[spacers.index(one)] += f"+{des}"
        tot = len(spacers)

        #build graph -- using adjacency matrix
        adjmatrix = []
        for idx in range(tot):
                tmp = [0] * tot
                adjmatrix.append(tmp)

        print("spacers", spacers)
        for aline in lines:
                subs = aline.strip().split()
                one, des = getspacer(subs[2])
                pre = spacers.index(one)
                for one0 in subs[3:-1]:
                        one, des = getspacer(one0)
                        curr = spacers.index(one)
                        adjmatrix[pre][curr] += 1
                        pre = curr

        spacergraph = Graph(spacers, adjmatrix)
        if collapse:
            spacergraph.collapse()
        else:
            spacergraph.nocollapse()
        tmplist = []
        for one0 in lines[0].split()[2:-1]:
            one, des = getspacer(one0)
            tmplist.append(one)
        spacergraph.getpath(tmplist)
        spacergraph.writedot(dotfile, spacers_detail=spacers_detail, colorfile=colorfile)
main()
