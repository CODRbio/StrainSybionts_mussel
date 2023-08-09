#!/usr/bin/env python3
import sys, getopt, os
from collections import Counter, defaultdict
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
from ete3 import EvolTree
#from pathos.pools import ProcessPool 
from multiprocessing import Pool
from statsmodels.sandbox.stats.multicomp import multipletests
import random
import statistics
import shutil

surfix = ".fa.aln_adj_cds.fa"

def usage():
    print('''
    -h help
    -r analysis after discarding all previous results
    -a folders containing aligned cds in fasta format
    -t target species, seperated by comma
    -s SpeciesID map table
    -S surfix of aligned cds files, default .fa.aln_adj_cds.fa
    -T species tree
    -f dnds_res for the clades to filter unsutiable genes
    -n nodes to label without tips
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

def targetINort(speRep,target,ortho,tree):
    tarlist = target.split(",")
    resList = []
    species = ["T_"+x for x in list(speRep.keys())]
    for item in tarlist:
        tarLeaves = get_leafs(item,tree)
        #species = ["T_"+x for x in list(speRep.keys())]
        tarInter = interList(species, tarLeaves)
        resList.append(tarInter)
    return(resList)

def Align_tre_gen(speRep, seqDic, refTree, ortho): 
    species = ["T_"+x for x in list(speRep.keys())]
    AllLeaves = refTree.get_leaf_names()
    if(set(AllLeaves) <= set(species)):
        rec_out = []
        for Rec in speRep:
            rec_new = SeqRecord(Seq(seqDic[speRep[Rec]]),id="".join(["T_",Rec]))
            rec_out.append(rec_new)
        basename=os.path.split(ortho)[1]
        basename=basename.replace(surfix,"")
        basename = os.path.splitext(basename)[0]
        OrtOut=basename + ".phy"
        SeqIO.write(rec_out, OrtOut, "phylip-sequential")
        return(basename)

def permulation_gen(seqList, num_sampling, ids):
    count = len(seqList)
    sel = random.sample(range(count), num_sampling)
    sample = [seqList[x]+".phy" for x in sel]
    sequences_con = defaultdict(str)
    for seqs in sample:
        align=AlignIO.read(seqs, "phylip")
        for record in align:
            ID = record.id
            sequences_con[ID] += str(record.seq)
    rec_out = []
    for key in sequences_con:
        rec = SeqRecord(Seq(sequences_con[key]),id=key)
        rec_out.append(rec)
    AlignOut = "permulation_{:04d}.phy".format(ids)
    SeqIO.write(rec_out, AlignOut, "phylip-sequential")
    return(AlignOut)


def mark_tree(t,target,speTree,label):
    tarLeaves = get_leafs(target,speTree)
    species = t.get_leaf_names()
    tarInter = interList(species, tarLeaves)
    ancestor = t.get_common_ancestor(tarInter)
    markList = [x.node_id for x in ancestor.get_descendants()] + [ancestor.node_id]
    marks = [label] * len(markList)
    t.mark_tree (markList, marks = marks)
    return(t)

def mark_node(t,target,speTree,label):
    tarLeaves = get_leafs(target,speTree)
    species = t.get_leaf_names()
    tarInter = interList(species, tarLeaves)
    ancestor = t.get_common_ancestor(tarInter)
    markList = [ancestor.node_id]
    marks = [label] 
    t.mark_tree (markList, marks = marks)
    return(t)

def res_parser(res,basename,target,tarNode,model="b_free"):
    lnL = res.lnL
    resBuffer = {}
    NdN = {}
    SdS = {}
    dS = {}
    omega = {}
    if model == "b_free":
        for key in res.branches:
            if(int(res.branches[key]['mark'].split('#')[1]) <20 and int(res.branches[key]['mark'].split('#')[1])>0):
                nodes = target.split(",")[(int(res.branches[key]['mark'].split('#')[1])-1)]
            elif(int(res.branches[key]['mark'].split('#')[1])>=20):
                nodes = tarNode.split(",")[(int(res.branches[key]['mark'].split('#')[1])-21)]
            else:
                nodes = "root"
                continue
            mark = res.branches[key]['mark']
            if nodes not in NdN:
                NdN[nodes] = []
            NdN[nodes].append(res.branches[key]['dN']*res.branches[key]['N'])
            if nodes not in SdS:
                SdS[nodes] = []
            SdS[nodes].append(res.branches[key]['dS']*res.branches[key]['S'])
            if nodes not in dS:
                dS[nodes] = []
            dS[nodes].append(res.branches[key]['dS'])
            omega.update({nodes:res.branches[key]['w']})
        for nodes in set(NdN.keys()) - {"root"}:
            resBuffer[nodes]=[nodes, mark, basename, max(NdN[nodes]), max(SdS[nodes]), statistics.median(dS[nodes]), omega[nodes], lnL]
    else:
        nodes = "M0All"
        for key in res.branches:
            if 'dN' not in res.branches[key]:
                continue
            if nodes not in NdN:
                NdN[nodes] = []
            NdN[nodes].append(res.branches[key]['dN'] * res.branches[key]['N'])
            if nodes not in SdS:
                SdS[nodes] = []
            SdS[nodes].append(res.branches[key]['dS'] * res.branches[key]['S'])
            if nodes not in dS:
                dS[nodes] = []
            dS[nodes].append(res.branches[key]['dS'])
            omega.update({nodes:res.branches[key]['w']})
        resBuffer[nodes]=[nodes, "No_mark", basename, max(NdN[nodes]), max(SdS[nodes]), statistics.median(dS[nodes]), omega[nodes], lnL]
    resLis = [x for x in resBuffer.values()]
    return(resLis)
def dnds_cal(tree, align, target, tarNode):
    basename=os.path.split(align)[1]
    basename = basename.replace(surfix, "")
    basename=os.path.splitext(basename)[0]
    tree.link_to_alignment(align)
    try:
        tree.link_to_evol_model(f"/tmp/ete3-tmp/b_free.{basename}/out", 'b_free.' + basename)
    except:
        tree.run_model ('b_free.'+basename, **{'Small_Diff':1e-5})
    res = tree.get_evol_model('b_free.'+basename)
    resUP = res_parser(res,basename,target,tarNode,model="b_free")
    return(resUP)

def re_test(tree, align, target, tarNode):
    basename=os.path.split(align)[1]
    basename=basename.replace(surfix,"")
    basename = os.path.splitext(basename)[0]
    tree.link_to_alignment(align)
    try:
        shutil.rmtree(f"/tmp/ete3-tmp/b_free.{basename}/")
    except:
        print(f"No {basename} available")
    tree.run_model ('b_free.'+basename, **{'Small_Diff':1e-5})
    res = tree.get_evol_model('b_free.'+basename)
    resUP = res_parser(res,basename,target,tarNode,model="b_free")
    return(resUP)

def get_align(path):
    name = []
    for files in os.listdir(path):
        if files.endswith(surfix):
            name.append(os.path.join(path,files))
    return name

def get_permulation(path):
    name = []
    for files in os.listdir(path):
        if files.startswith("permulation_") and files.endswith(".phy"):
            name.append(files)
    return name

def get_tree(path):
    name = []
    for files in os.listdir(path):
        if files.endswith([".tre",".tree","newick"]):
            name.append(os.path.join(path,files))
    return name

redo = False
tarNode = False
output = "RES"
opts,args=getopt.getopt(sys.argv[1:],'-h-f:-n:-t:-a:-s:-T:-o:-S:-r',['help'])
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
    elif k in ("-f"):
        filterF = v
    elif k in ("-n"):
        tarNode = v
    elif k in ("-S"):
        surfix = v
    else:
        usage()
        sys.exit()
filterDF = pd.read_table(filterF,sep="\t",header=0)

filterDF["passed"] = (filterDF['w']<=5) | ((filterDF['w']>5) & (filterDF["N*dN"] >0.5) & (filterDF["S*dS"] >0.5) & (filterDF["dS"]<2))
tmplis = [not(x) for x in filterDF["passed"]]
tmpDF = filterDF[tmplis]
filterDic = {}
for index, row in tmpDF.iterrows():
    filterDic[row["Gene"]] = 1

speTree = tree
mapper = id_mapper(SpeMap)
aligns = get_align(alignment)
#tarList = tarSpe.split(',')  # leave for the future works with multiple spcies as forground
refTree = EvolTree(speTree, format=1)
count = 1 
for tarItem in tarSpe.split(","):
    label = "#"+str(count)
    count = count + 1
    refTree= mark_tree(refTree,tarItem,speTree,label)
if(tarNode):
    count = 21
    for tarItem in tarNode.split(","):
        label = "#"+str(count)
        count = count + 1
        refTree= mark_node(refTree,tarItem,speTree,label)

print("Preparing necessary files for selection Test")
files = os.listdir(".")
OrthoSelected = []
for orthologs in aligns:
    basename=os.path.split(orthologs)[1]
    basename=basename.replace(surfix,"")
    if(basename in filterDic):
        continue
    if(not os.path.exists(basename+".phy")):
        FamilyInfo = get_repreSeq(ortho=orthologs)
        task = Align_tre_gen(FamilyInfo[0], FamilyInfo[1], refTree, orthologs) 
        if task is not None:
            OrthoSelected.append(task)
    else:
        if(os.path.exists(f"{basename}.phy")):
            OrthoSelected.append(basename)

PermuSeqs = []
if(redo):
    for num in range(1000):
        PermuSeqs.append(permulation_gen(OrthoSelected, 25, num))

print("Running Selection Test")
results = []
PerAlign = get_permulation("./")
align=sorted(PerAlign)
target= [tarSpe] * len(align)
tarNodeL = [tarNode] * len(align)
refTreeL = [refTree] * len(align)
pool = Pool() 
resLis = []
zip_args = list(zip(refTreeL, align, target, tarNodeL))
if(redo):
    results = pool.starmap(re_test,zip_args,chunksize = 1)
    pool.close()
    pool.join()
else:
    results = pool.starmap(dnds_cal,zip_args,chunksize = 1)
    pool.close()
    pool.join()
for level1 in results:
    for level2 in level1:
        resLis.append(level2)
df = pd.DataFrame(resLis, columns = [ "Nodes", "Mark", "Gene", "N*dN", "S*dS", "dS", "w", "lnL" ])
df.to_csv(output+"permulation_pnps_each_clade.tsv",sep='\t', index=0)
