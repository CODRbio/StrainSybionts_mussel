#!/nfs_genome/anaconda/envs/rnaseq/bin/python
import sys, getopt, os
import pandas as pd
import numpy as np
from collections import Counter, defaultdict

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

surfix=("fas","fa","faa","pep","fasta","aa","out","trimmed")

def usage():
    print('''
    -h help
    -a folders containing aligned seqs in fasta format
    -O Orthogroup_count.tsv file
    -s SpeciesID map table
    -S surfix of the aligned file, the part following OG00XXXX, .fa.aln_adj_cds.fa for cds
    -o output prefix
    -p prefix of the Ortholog
    ''')

def GetMulticopyCutoff(nSpecies, factor = 0.25, pMax = 0.25):
    """ Get the maximum number of species with multicopy genes for each number of non-single copy genes. A species is non-single
    copy if it is absent or if it is multicopy
    Args:
        nSpecies - number of species in the analysis
        factor - probability that a single copy gene is actually non-orthologous relative to the proportion of multi-copy species
        pMax - maximum allowed probability that one of the genes is non-orthologous
    Returns:
        maxMulticopy - array, i-th element is maximum number of allowed multicopy species if ispecies are non-single copy
    """        
    allowed_multicopy = []
    for iNonSingle in range(nSpecies):
        nSingle = nSpecies - iNonSingle
        qFoundMax = False
        for nMultiple in range(iNonSingle+1):
            pFalse = factor*nMultiple/float(nSpecies)
            pAnyFalse = 1.-(1-pFalse)**nSingle
            if pAnyFalse > pMax:
                allowed_multicopy.append(nMultiple - 1)
                qFoundMax = True
                break
        if not qFoundMax:
            allowed_multicopy.append(iNonSingle)
    return allowed_multicopy

def SingleCopy_WithProbabilityTest(fraction, ogMatrix):
    """ X fraction of species have a single copy of the gene, all other species can have whatever number of copies including 0"""
    nSpecies = ogMatrix.shape[1]
    allowed_multicopy = GetMulticopyCutoff(nSpecies)
    multi = (ogMatrix > 1 * np.ones((1, nSpecies))).sum(1) # multicopy
    excluded = nSpecies - (ogMatrix == 1 * np.ones((1, nSpecies))).sum(1) # single-copy
    temp2 = (ogMatrix == np.ones((1, nSpecies))).sum(1) / float(nSpecies)   # proportion of species with gene present
    c2 = temp2 > fraction
    passFracOGs = c2.values.nonzero()[0]
    ok_ogs = [iog for iog in passFracOGs if multi[iog] <= allowed_multicopy[excluded[iog]]]
    return ok_ogs


def GetOrthogroupOccupancyInfo(m):
    """
    Args:
        m - orthogroup matrix
    """
    N = m.shape[1]
    f = 1./N
    fractions = []
    nOrtho = []
    for n in range(N):
        F = 1.-n*f
        fractions.append(F)
        nOrtho.append(len(SingleCopy_WithProbabilityTest(F-1e-5, m)))
    nOrtho = list(map(float, nOrtho))
    return fractions, nOrtho

def DetermineOrthogroupsForSpeciesTree(m, nOGsMin=100, nSufficient=1000, increase_required=2.):
    """Orthogroups can be used if at least a fraction f of the species in the orthogroup are single copy, f is determined as described 
    below. Species that aren't single copy are allowed in the orthogroup as long as there aren't too many (criterion is quite strict). 
    The presence of species with multiple copies suggests that the single copy species might actually have arisen from duplication 
    and loss. Therefore, we use a probability model to determine how many of the excluded species can be multicopy before it is likely
    that the single copy genes are actually hidden paralogues.
    Args:
        m - the orthogroup matrix shape=(nOGs, nSpecies), each entry is number of genes from that species
        nOGsMin - the minimum number of orthogroups. Need to have at least this many before considering the next criteria
        p - proportionalIncreaseMin: Each decrease in the proportion of single copy species should bring about a relative increase 
        in the number of orthogroups that will be used of 'p', otherwise it is not worth lowering the bar
    """
    fractions, nOrtho = GetOrthogroupOccupancyInfo(m)
    if nOrtho[0] > nSufficient:
        ogsToUse = SingleCopy_WithProbabilityTest(1.0-1e-5, m)
        return ogsToUse, 1.0
    nSpecies = m.shape[1]
    for i in range(1, nSpecies):
        if nOrtho[i-1] < nOGsMin: continue
        if nOrtho[i-1] > nSufficient: break
        p = ((nOrtho[i] - nOrtho[i-1])/nOrtho[i-1]) /  (-(fractions[i] - fractions[i-1])/fractions[i-1])
        if fractions[i] > 0.5:
            if p < increase_required: break
        else:
            # if fewer than half the species are in any orthogroup then set highr bar for reducing fraction further
            if p < 2*increase_required: break
    f = fractions[i-1]
    ogsToUse = SingleCopy_WithProbabilityTest(f-1e-5, m)
    return ogsToUse, f   

def id_mapper(SpeMap):
    mapper = {}
    IN=open(SpeMap)
    for lines in IN.readlines():
        ID = lines.split()[0].replace(":","")
        SeqID = lines.split()[1].replace(".faa","")
        mapper[SeqID] = ID
    return(mapper)


def ortho_mapper(id, matrix):
    name = matrix.iloc[:,0][id]
    try:
        name = prefix + name
        return(name)
    except:
        return(name)

def GetSingleOrthoSeqs(orthoID, m,alignment,mapper,surfix):
    seqDic = {}
    sogID = m.loc[orthoID][m.loc[orthoID,]==1].keys()
    print("There are %d SOGs" %len(sogID))
    sogID_tran =[]
    for Sid in sogID: 
        sogID_tran.append(mapper[Sid])
    sogDis_tran = []
    sogDis = list(set(m.loc[orthoID][m.loc[orthoID,]!=1].keys()))
    print("There are %d SEGs" %(len(sogID)+len(sogDis)))
    for Sid in sogDis:
        sogDis_tran.append(mapper[Sid])
    orthoName = ortho_mapper(orthoID,matrix)
    print(orthoName)
    aliFile = os.path.join(alignment,orthoName)+surfix
    try:
        align=AlignIO.read(aliFile, "fasta")
    except:
        print("The alignment %s cannot be opened" %aliFile) #s.path.join(alignment,orthoName,aliFile))
        return(0)
    length=align.get_alignment_length()
    tmpDic = {}
    for record in align:
        speID = record.id.split("_")[0]
        tmpDic[speID] = str(record.seq)

    for speID in sogID_tran:
        #print("seq%s:%s" % (speID,record.seq))
        try:
            seqDic[speID] = tmpDic[speID]
        except:
            seq = "-" * length
            seqDic[speID] = str(seq)
    for Sid in sogDis_tran:
        seq = "-" * length
        #print("seq%s:%s" %(Sid, seq))
        seqDic[Sid] = str(seq)
    return(seqDic)
prefix=""
opts,args=getopt.getopt(sys.argv[1:],'-h-p:-a:-s:-S:-o:-O:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-o"):
        output = v
    elif k in ("-p"):
        prefix = v
    elif k in ("-a"):
        alignment = v
    elif k in ("-s"):
        SpeMap = v
    elif k in ("-S"):
        surfix = v
    elif k in ("-O"):
        Orthogroup = v
    else:
        usage()
        sys.exit()

mapper = id_mapper(SpeMap)
matrix = pd.read_table(Orthogroup,sep="\t",header=0)
nSpecies = matrix.shape[1] - 1
m=matrix.iloc[0:,1:nSpecies]

targetOrgNum = 1500
while DetermineOrthogroupsForSpeciesTree(m,nSufficient=targetOrgNum)[1] <=0.75:
    targetOrgNum = targetOrgNum - 50
OrgSel = DetermineOrthogroupsForSpeciesTree(m,nSufficient=targetOrgNum)[0]
concanated = {}
print(len(OrgSel))
for ID in OrgSel:  
    ortSel = GetSingleOrthoSeqs(ID,m,alignment,mapper,surfix)
    if ortSel:
        for rec in ortSel:
            try:
                concanated[rec] = concanated[rec] + ortSel[rec]
            except:
                concanated[rec] = ortSel[rec]

conseq = open("_".join([output,"concanated.fa"]), "w")
for key in sorted(concanated):
    conseq.write(">"+key+"\n")
    conseq.write(concanated[key]+"\n")
