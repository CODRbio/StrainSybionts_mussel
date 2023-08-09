#!/usr/bin/env python3
import sys, getopt, os
from collections import Counter, defaultdict
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from ete3 import EvolTree
from pathos.pools import ProcessPool 
from statsmodels.sandbox.stats.multicomp import multipletests


surfix=("fas","fa","faa","pep","fasta","aa","phy")

def usage():
    print('''
    -h help
    -r reanalysis using pre-computed tmp-results
    -a folders containing aligned cds in fasta format
    -t target species
    -s SpeciesID map table
    -T species tree
    -o output prefix
    ''')


def id_mapper(SpeMap):
    mapper = {}
    IN=open(SpeMap)
    for lines in IN.readlines():
        ID = lines.split()[0].replace(":","")
        SeqID = lines.split()[1].replace(".faa","")
        mapper[SeqID] = ID
    return(mapper)
#speRep  the representative seq id for the sepcies in the ortholog
#seqDic  the sequences of each id
#speOrt  the seq ids for each species in one ortholog
def get_repreSeq(ortho):  #获得ortho中每个物种的代表性序列 
    align=AlignIO.read(ortho, "fasta")
    seqDic = {}
    speOrt = defaultdict(list)
    lenDic = {}
    speRep = {}
    for record in align:
        speID = record.id.split("_")[0]
        OrtID = record.id
        seqDic[OrtID] = str(record.seq)
        lenDic[OrtID] = len(str(record.seq).replace("-",""))
        speOrt[speID].append(OrtID)
    for key in speOrt:
        if len(speOrt[key]) > 1:
            OrtList = sorted(speOrt[key], reverse=True, key=lambda x: lenDic[x])
            speRep[key] = OrtList[0]
        else:
            speRep[key] = speOrt[key][0]
    return(speRep,seqDic)

def get_leafs(target,tree):
    t = EvolTree(tree,format=1)
    node = t&target
    TarLeafList = node.get_leaf_names()
    return(TarLeafList)

def interList(lists, AllList):
    if type(lists) == list:
        pass
    else:
        lists = [ lists ]
    return(list(set(AllList) & set(lists)))



def Align_tre_gen(speRep, seqDic, targetSpe, ortho, tree):
    OrtIDconver = {}
    tarLeaves = get_leafs(targetSpe,tree)
    species = ["T_"+x for x in list(speRep.keys())] 
    tarInter = interList(species, tarLeaves)
# 1, 3 , 5
    if len(tarInter) >= 1 and (len(species) - len(tarInter)) >=2 and len(species) >=4:
        #print("target:%d, all:%d" %(len(tarInter), len(species)))
        rec_out = []
        for Rec in speRep:
            rec_new = SeqRecord(Seq(seqDic[speRep[Rec]]),id="".join(["T_",Rec]))
            rec_out.append(rec_new)
        basename=os.path.split(ortho)[1]
        basename=os.path.splitext(basename)[0]
        OrtOut=basename + ".phy"
        SeqIO.write(rec_out, OrtOut, "phylip-sequential")
        t = EvolTree(tree, format=1)
        t.prune(species, preserve_branch_length=False)
        treefile=basename+".tre"
        t.write(outfile=treefile)
        return(basename)

def re_test(align):
    basename=os.path.split(align)[1]
    basename=os.path.splitext(basename)[0]
    t = EvolTree(basename+".tre")
    t.link_to_alignment(basename+".phy")
    t.link_to_evol_model("/tmp/ete3-tmp/bsA."+basename+"/out", 'bsA')
    t.link_to_evol_model("/tmp/ete3-tmp/bsA1."+basename+"/out", 'bsA1')
    t.link_to_evol_model("/tmp/ete3-tmp/M0."+basename+"/out", 'M0')
    ps = t.get_most_likely ('bsA', 'bsA1')
    rx = t.get_most_likely ('bsA1', 'M0')
    if ps<0.05:
        status = "positive selection"
    elif rx<0.05 and ps>=0.05:
        status = "relaxation"
    else:
        status = "Fail"
    print(basename+"\t"+str(ps)+"\t"+str(rx)+"\t"+status)
    return([basename,ps,rx,status])

def bs_test(tree, align, target, speTree):
    try:
        t= EvolTree(tree)
    except:
        print("The tree %s cannot be opened" %tree)
    basename=os.path.split(tree)[1]
    basename=os.path.splitext(basename)[0]
    t.link_to_alignment(align)
    tarLeaves = get_leafs(target,speTree)
    align = AlignIO.read(align, "phylip")
    species = []
    for record in align:
        OrtID = record.id
        species.append(OrtID)
    try:
        tarInter = interList(species, tarLeaves)
        ancestor = t.get_common_ancestor(tarInter)
        markList = [x.node_id for x in ancestor.get_descendants()] + [ancestor.node_id] 
        marks = ['#1'] * len(markList)
    except:
        return([basename,1,1,"Failed"])
    t.mark_tree (markList, marks = marks)
    t.run_model ('M0.'+basename)
    t.run_model ('bsA.'+basename)
    t.run_model ('bsA1.'+basename)
    ps = t.get_most_likely ('bsA.'+basename, 'bsA1.'+basename) 
    rx = t.get_most_likely ('bsA1.'+basename, 'M0.'+basename)
    
    if ps<0.05:
        status = "positive selection"
    elif rx<0.05 and ps>=0.05:
        status = "relaxation"
    else:
        status = "Fail"
    print(basename+"\t"+str(ps)+"\t"+str(rx)+"\t"+status)
    return([basename,ps,rx,status])

def get_align(path):
    name = []
    for files in os.listdir(path):
        if files.endswith(surfix):
            name.append(os.path.join(path,files))
    return name

def get_tree(path):
    name = []
    for files in os.listdir(path):
        if files.endswith([".tre",".tree","newick"]):
            name.append(os.path.join(path,files))
    return name

redo = False
opts,args=getopt.getopt(sys.argv[1:],'-h-t:-a:-s:-T:-o:-r',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-o"):
        output = v
    elif k in ("-t"):
        tarSpe = v
    elif k in ("-a"):
        alignment = v
    elif k in ("-s"):
        SpeMap = v
    elif k in ("-T"):
        tree = v 
    elif k in ("-r"):
        redo = True
    else:
        usage()
        sys.exit()

speTree = tree
mapper = id_mapper(SpeMap)
aligns = get_align(alignment)
#tarList = tarSpe.split(',')  # leave for the future works with multiple spcies as forground
if not redo:
    print("Preparing necessary files for selection Test")
    OrthoSelected = []
    for orthologs in aligns:
        basename=os.path.split(orthologs)[1]
        basename=os.path.splitext(basename)[0]
        if not os.path.exists(basename+".phy"):
            FamilyInfo = get_repreSeq(ortho=orthologs)
            task = Align_tre_gen(FamilyInfo[0], FamilyInfo[1], tarSpe, orthologs, tree) 
            if task is not None:
                OrthoSelected.append(task)
        else:
            OrthoSelected.append(basename)
    print("Running Selection Test")
    results = []
    tree=sorted(list(map(lambda x:x+".tre", OrthoSelected )))
    align=sorted([ x+".phy" for x in OrthoSelected ])

#tNode = ["T_"+str(tarSpe)]
    target= [tarSpe] * len(tree)
    speTreeL = [speTree] * len(tree)

    pool = ProcessPool()
#print(target)

    results = pool.map(bs_test,tree, align, target, speTreeL)
    df = pd.DataFrame(results, columns = [ "orthoID", "p_ps", "p_rx", "status" ])
    p_ps = df['p_ps'].tolist()
    p_rx = df['p_rx'].tolist()
    ps_adjusted = multipletests(p_ps, alpha=.05, method='bonferroni')[1]
    rx_adjusted = multipletests(p_rx, alpha=.05, method='bonferroni')[1]
    df['padj_ps'] = ps_adjusted
    df['padj_rx'] = rx_adjusted
    
    df_ps = df[df['padj_ps']<0.05]
    df_rx = df[(df['padj_ps']>0.05) & (df['padj_rx']<0.05)]
    df_ps.to_csv(output+"PS_selection.tsv",sep='\t', index=0)
    df_rx.to_csv(output+"RX_selection.tsv",sep='\t', index=0)
            
else:
    results = []
    aligns = get_align("./")
    for orthologs in aligns:
        results.append(re_test(orthologs))
    df = pd.DataFrame(results, columns = [ "orthoID", "p_ps", "p_rx", "status" ])
    p_ps = df['p_ps'].tolist()
    p_rx = df['p_rx'].tolist()
    ps_adjusted = multipletests(p_ps, alpha=.05, method='bonferroni')[1]
    rx_adjusted = multipletests(p_rx, alpha=.05, method='bonferroni')[1]
    df['padj_ps'] = ps_adjusted
    df['padj_rx'] = rx_adjusted

    df_ps = df[df['padj_ps']<0.05]
    df_rx = df[(df['padj_ps']>0.05) & (df['padj_rx']<0.05)]
    df_ps.to_csv(output+"PS_selection.tsv",sep='\t', index=0)
    df_rx.to_csv(output+"RX_selection.tsv",sep='\t', index=0)
