#!/usr/bin/env python3
import sys, getopt, os
from scipy import stats

def usage():
    print('''
    use to find the sepecific orthologs in each clade
    -t adjusted tree file with node names
    -T targeted clades to find specific orthologs
    -O Orthogroup gene count table
    -p adjusted pvalue
    -s speciesId in the workingdir of orthofinder
    -o output
    ''')
    
padj = False  # assign a default value
opts,args=getopt.getopt(sys.argv[1:],'-h-t:-p-O:-s:-T:-o:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-t"):
        tree = v
    elif k in ("-T"):
        target = v
    elif k in ("-O"):
        Orthogroup = v
    elif k in ("-s"):
        speMap = v
    elif k in ("-o"):
        output  = v
    elif k in ("-p"):
        padj = True
    else:
        usage()
        sys.exit()

if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

from ete3 import EvolTree
import pandas as pd
from statsmodels.stats.multitest import multipletests
def get_leafs(target,tree):    
    t = EvolTree(tree,format=1)
    node = t&target
    TarLeafList = node.get_leaf_names()
    return(TarLeafList)


def ortho_mapper(id, matrix):
    name = matrix.iloc[:,0][id]
    return(name)

def id_mapper(SpeMap):
    mapper = {}
    IN=open(SpeMap)
    for lines in IN.readlines():
        ID = lines.split()[0].replace(":","")
        ID = "T_"+str(ID)
        SeqID = lines.split()[1].replace(".faa","")
        mapper[ID] = SeqID
    return(mapper)


def contigency_TabGen(tarList, bgList , Matrix):
    a = (Matrix[tarList] > 0).sum(1)
    # Cell a: occurrences of the gene of interest in the species of interest
    c = (Matrix[bgList] > 0).sum(1)
    # Cell c: occurrences of the gene of interest in other species
    b = [sum(a.drop(x)) for x in a.index]
    # Cell b: occurrences of other genes in the species of interest
    d = [sum(c.drop(x)) for x in c.index]
    # Cell d: occurrences of other genes in other species
    list2x2 = [ [[a[x], b[x]], [c[x], d[x]]] for x in a.index]
    return list2x2

def FisherTest(contingency_list):
    fisherP = []
    for contingency in contingency_list:
        oddsratio, pvalue = stats.fisher_exact(contingency, alternative="greater")
        fisherP.append(pvalue)
    if padj:
        #fisherP_adj = stats.p_adjust(FloatVector(fisherP), method = 'BH')
        fisherP_adj = multipletests(fisherP, method='fdr_bh')[1]
        return(fisherP_adj)
    else:
        return(fisherP)

def FindSpecific(tarNode,tarList,Matrix):  #Fisher needes the gene frequency in front clade significantly higher than the other clade of interest and all tips other than front taxa, however, the uniq genes only need the front gene frequency is significantly higher than background and the genes occur in front clade only
    pMatrix = []
    otherNode = list(set(tarList) - set([tarNode])) 
    ftList = list(map(lambda x: SpeMapper[x], get_leafs(tarNode, tree)))
    if len(otherNode)>0:
        for bgNode in otherNode:
            bgList = list(map(lambda x: SpeMapper[x], get_leafs(bgNode,tree)))
            contigencyTab = contigency_TabGen(ftList, bgList , Matrix)
            pMatrix.append(FisherTest(contigencyTab))
    bgList = list(set(Matrix) - set(ftList))
    contigencyTab = contigency_TabGen(ftList, bgList, Matrix)
    pMatrix.append(FisherTest(contigencyTab))
    df = pd.DataFrame(list(map(list, zip(*pMatrix))))
    df.to_csv(f"fisherHigh_{tarNode}.tsv",sep="\t")
    df_b = df <= 0.05
    fisher_filtered = [ortho_mapper(id, matrixDF) for id in df_b[df_b.all(axis=1)].index]
    #df_b = df_b.iloc[:,-1]
    df_b["unique"] = [(x[0][0] > 0 and x[1][0] == 0) for x in contigencyTab]
    df_b.to_csv(f"uniq_{tarNode}.csv",sep="\t")
    all_filtered = [ortho_mapper(id, matrixDF) for id in df_b[df_b.iloc[:,[-2,-1]].all(axis=1)].index]
    return(fisher_filtered, all_filtered)

if "," in target:
    tarList = target.split(',')
else:
    tarList = [target]
matrixDF = pd.read_table(Orthogroup,sep="\t",header=0)
nSpecies = matrixDF.shape[1] - 1
Matrix=matrixDF.iloc[:,1:nSpecies]
SpeMapper = id_mapper(speMap)

for tarNode in tarList:
    #leaves = get_leafs(tarNode,tree)
    FisherRes,SpecificRes = FindSpecific(tarNode, tarList, Matrix)
    fo = open(output+"_"+tarNode+"_FisherRes.tsv", "w")
    fo.write("\n".join(FisherRes))
    fo.close()
    fo = open(output+"_"+tarNode+"_UnicRes.tsv", "w")
    fo.write("\n".join(SpecificRes))
    fo.close()

