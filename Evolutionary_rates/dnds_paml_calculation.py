#!/usr/bin/env python3
import sys, getopt, os
from collections import Counter, defaultdict
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
from ete3 import EvolTree
from multiprocessing import Pool
import statistics
from statsmodels.sandbox.stats.multicomp import multipletests



def usage():
    print('''
    -h help
    -r reanalyze the data even when the pre-computed tmp-results are available
    -a folders containing aligned cds in fasta format
    -t target species, seperated by comma
    -s SpeciesID map table
    -S surfix, default .fa.aln_adj_cds.fa
    -T species tree
    -n nodes to label without tips (optional, used to know whether some ancstor node evolve in uniq rate)
    -o output prefix
    -f dnds_res for the clades to filter unsutiable genes (some clades may evolve in special rate, can be filtered, can be gained using M0 model first)
    -m model to use "M0 or b_free, default b_free"
    ''')

#use generate a dict to transfer the species to species_ID
def id_mapper(SpeMap):
    mapper = {}
    IN=open(SpeMap)
    for lines in IN.readlines():
        ID = lines.split()[0].replace(":","")
        SeqID = lines.split()[1].replace(".faa","")
        mapper[SeqID] = ID
    return(mapper)

#speOrt  the seq ids for each species in one ortholog
#to parse the longest seq in each ortholog, the representative id were stored in speRep
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

# get the species ID for each node of interest
def get_leafs(target,tree):
    t = EvolTree(tree,format=1)
    node = t&target
    TarLeafList = node.get_leaf_names()
    return(TarLeafList)

# get the intersection of two lists
def interList(lists, AllList):
    if type(lists) == list:
        pass
    else:
        lists = [ lists ]
    return(list(set(AllList) & set(lists)))

# give the interaction of the representative seq IDs of the target nodes in the species tree
def targetINort(speRep,target,tree):
    tarlist = target.split(",")
    resList = []
    species = ["T_"+x for x in list(speRep.keys())]
    for item in tarlist:
        tarLeaves = get_leafs(item,tree)
        #species = ["T_"+x for x in list(speRep.keys())]
        tarInter = interList(species, tarLeaves)
        resList.append(tarInter)
    return(resList)

# At least 2 species are expected in each node
def Align_tre_gen(speRep, seqDic, targetSpe, ortho, tree): 
    species = ["T_"+x for x in list(speRep.keys())]
    OrtLinTar = targetINort(speRep,targetSpe,tree) # the representative seq IDs of the target nodes in the species tree
    Cnt = map(lambda x:len(x),OrtLinTar)
    CntInOrt = map(lambda x:len(x)>1,OrtLinTar)
    if(all(CntInOrt) and (len(species)-sum(Cnt)>=0)):
        #print("target:%d, all:%d" %(len(tarInter), len(species)))
        rec_out = []
        for Rec in speRep:
            rec_new = SeqRecord(Seq(seqDic[speRep[Rec]]),id="".join(["T_",Rec]))
            rec_out.append(rec_new)
        basename=os.path.split(ortho)[1]
        basename=basename.replace(surfix,"")
        OrtOut=basename + ".phy"
        SeqIO.write(rec_out, OrtOut, "phylip-sequential")
        t = EvolTree(tree, format=1)
        t.prune(species, preserve_branch_length=False)
        treefile=basename+".tre"
        t.write(outfile=treefile)
        return(basename)

def mark_tree(t,target,speTree,label):
    tarLeaves = get_leafs(target,speTree)
    species = t.get_leaf_names()
    tarInter = interList(species, tarLeaves)
    ancestor = t.get_common_ancestor(tarInter)
    markList = [x.node_id for x in ancestor.get_descendants()] + [ancestor.node_id]
    marks = [label] * len(markList)
    t.mark_tree (markList, marks = marks)
    return(t)

def mark_node(t,target,speTree,label): # mark only some nodes, useful when try to find some uniq feature
    tarLeaves = get_leafs(target,speTree)
    species = t.get_leaf_names()
    tarInter = interList(species, tarLeaves)
    ancestor = t.get_common_ancestor(tarInter)
    markList = [ancestor.node_id]
    marks = [label] 
    t.mark_tree (markList, marks = marks)
    return(t)

def res_parser(res,basename,target,tarNode,model):
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


def dnds_cal(tree, align, target, tarNode, speTree, model):
    try:
        t= EvolTree(tree)
    except:
        print("The tree %s cannot be opened" %tree)
        sys.exit()
    basename=os.path.basename(tree)
    basename=os.path.splitext(basename)[0]
    print(f"Running {basename}")
    t.link_to_alignment(align)
    count = 1
    for tarItem in target.split(","):
        label = "#"+str(count)
        count = count + 1
        t= mark_tree(t,tarItem,speTree,label)
    if(tarNode):
        count = 21
        for tarItem in tarNode.split(","):
            label = "#"+str(count)
            count = count + 1
            t= mark_node(t,tarItem,speTree,label)
    try:
        t.link_to_evol_model(f"/tmp/ete3-tmp/{model}.{basename}/out", model + '.' + basename)
        res = t.get_evol_model(model + '.'+basename)
        resUP = res_parser(res,basename,target,tarNode,model)
    except:
        print(f"No previous records found, running {basename}!")
        try:
            t.run_model(model + '.'+basename, **{'Small_Diff':1e-4})
            res = t.get_evol_model(model + '.'+basename)
            resUP = res_parser(res,basename,target,tarNode,model)
        except:
            resUP=[["NA", "NA", basename, 0, 0, 0, 999, 99999]] 
    return(resUP)

def re_test(tree, align, target, tarNode, speTree, model):
    try:
        t= EvolTree(tree)
    except:
        print("The tree %s cannot be opened" %tree)
    basename=os.path.split(tree)[1]
    basename=os.path.splitext(basename)[0]
    t.link_to_alignment(align)
    count = 1
    for tarItem in target.split(","):
        label = "#"+str(count)
        count = count + 1
        t= mark_tree(t,tarItem,speTree,label)
    if(tarNode):
        count = 21
        for tarItem in tarNode.split(","):
            label = "#"+str(count)
            count = count + 1
            t= mark_node(t,tarItem,speTree,label)
    try:
        shutil.rmtree(f"/tmp/ete3-tmp/{model}.{basename}/")
    except:
        print(f"NO {basename} available")
    try:
        t.run_model(model + '.' + basename, **{'Small_Diff': 1e-6})
        res = t.get_evol_model(model + '.' + basename)
        resUP = res_parser(res, basename, target, tarNode, model)
        return (resUP)
    except:
        print(f"There are something wromg during the test for {basename}")
        return ([["NA", "NA", basename, 0, 0, 0, 999, 9999]])

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

def concanate_gen(seqList):
    sample = [x+".phy" for x in seqList]
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
    AlignOut = "concanated.phy"
    SeqIO.write(rec_out, AlignOut, "phylip-sequential")

#initiation of the variables
redo = False
tarNode = False
output = "RES"
model = "b_free"
surfix = ".fa.aln_adj_cds.fa"
filterF = False

opts,args=getopt.getopt(sys.argv[1:],'-h-f:-n:-t:-a:-s:-T:-o:-m:-S:-r',['help'])
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
    elif k in ("-n"):
        tarNode = v
    elif k in ("-m"):
        model = v
    elif k in ("-f"):
        filterF = v
    else:
        usage()
        sys.exit()

speTree = tree
mapper = id_mapper(SpeMap)
aligns = get_align(alignment) #get the cds alignment files
#if fiter option is selected, filter the data
if(filterF):
    filterDF = pd.read_table(filterF,sep="\t",header=0)
    filterDF["passed"] = (filterDF['w']<=5) | ((filterDF['w']>5) & (filterDF["N*dN"] >0.5) & (filterDF["S*dS"] >0.5) & (filterDF["dS"]<2))
    FailLis = [not(x) for x in filterDF["passed"]]
    failDic = {}
    failDF = filterDF[FailLis]
    passDic = {}
    for index, row in failDF.iterrows():
        failDic[row["Gene"]] = 1

print("Preparing necessary files for selection Test")

t = EvolTree(speTree,format=1)
cnt = len(t.get_leaf_names())
to_con = []
OrthoSelected = []
for orthologs in aligns:
    basename=os.path.split(orthologs)[1]
    basename=basename.replace(surfix,"")
    try:
        if(basename in failDic): 
            continue 
    except:
        print(f"{basename} passed")
    if(not os.path.exists(basename+".phy")):
        FamilyInfo = get_repreSeq(ortho=orthologs)
        task = Align_tre_gen(FamilyInfo[0], FamilyInfo[1], tarSpe, orthologs, tree)
        if (len(FamilyInfo[0]) == cnt):
            to_con.append(basename)
        if task is not None:
            OrthoSelected.append(task)
    else: 
        if(os.path.exists(f"{basename}.phy")):
            OrthoSelected.append(basename)

if(not os.path.exists("concanated.phy")):
    concanate_gen(to_con)
    t = EvolTree(speTree, format=1)  
    treefile="concanated.tre"
    t.write(outfile=treefile)

if ("concanated" not in OrthoSelected):
    OrthoSelected.append("concanated")


print("Running Selection Test")
results = []

tree=sorted(list(map(lambda x:x+".tre", OrthoSelected )))
align=sorted([ x+".phy" for x in OrthoSelected ])


#tNode = ["T_"+str(tarSpe)]
target= [tarSpe] * len(align)
tarNodeL = [tarNode] * len(align)
speTreeL = [speTree] * len(align)
modelL = [model] * len(align)
pool = Pool()
resLis = []
zip_args = list(zip(tree, align, target, tarNodeL, speTreeL, modelL))
if(redo):
    results = pool.starmap(re_test, zip_args, chunksize=1)
    pool.close()
    pool.join()
else:
    results = pool.starmap(dnds_cal, zip_args, chunksize=1)
    pool.close()
    pool.join()
for level1 in results:
    for level2 in level1:
        resLis.append(level2)
df = pd.DataFrame(resLis, columns = [ "Nodes", "Mark", "Gene", "N*dN", "S*dS", "dS", "w", "lnL" ])
df.to_csv(output+"_pnps_each_clade.tsv",sep='\t', index=0)
