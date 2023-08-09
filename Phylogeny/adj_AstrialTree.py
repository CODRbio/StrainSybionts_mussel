#!/usr/bin/env python3
import sys, getopt, os


def usage():
    print('''
    used for prepare trees from astrial compatible trees (index numbers for leave names)
    -t tree file
    -o output file
    ''')

opts,args=getopt.getopt(sys.argv[1:],'-h-t:-a:-s:-T:-o:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-t"):
        tree = v
    elif k in ("-o"):
        output  = v
    else:
        usage()
        sys.exit()


from ete3 import EvolTree


def adjust_nodename(tree):
    t = EvolTree(tree)
    for node in t.traverse():
        if node.is_leaf():
            if node.name.isdigit():
                node.name = "T_" + str(node.name)
        else:
            node.name = "N_" + str(node.node_id)
    return(t)

adj = adjust_nodename(tree)
adj.write(outfile=output, format=1)


