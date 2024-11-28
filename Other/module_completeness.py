#!/usr/bin/env python
import re

import sys, getopt, os


def usage():
    print('''
    use to find the sepecific orthologs in each clade
    -d directory containing *.Module_chart.csv
    -m module defination file, default /nfs_genome/BioDB/kegg/202212/ftp.kegg.net/kegg/module/module
    -o output
    -s suffix, default _Module_chart.csv
    ''')
suffix = "_Module_chart.csv"
opts,args=getopt.getopt(sys.argv[1:],'-h-s:-m:-d:-M:-o:',['help'])
module = "/nfs_genome/BioDB/kegg/202212/ftp.kegg.net/kegg/module/module"
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-M"):
        Module_ko = v
    elif k in ("-m"):
        module = v
    elif k in ("-o"):
        output = v
    elif k in ("-d"):
        directory = v
    elif k in ("-s"):
        suffix = v
    else:
        usage()
        sys.exit()
if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

def ko_hierachy_split(string):
    string = string.strip()
    status = 0
    result=[]
    for ch in string:
        if ch =="(":
            status += 1
        if ch == ")":
            status -= 1
        if status > 0:
            if re.search("\s+",ch):
                ch = "_"
        result.append(ch)
        transformed = "".join(result)
    return(transformed.split(" "))

def get_MO(mo): #get definition of Module form a list
    Mid = ""
    DEF = {}
    line = 0
    Name_Module = {}
    
    while line < len(mo):
        mo[line] = mo[line].strip()
        if re.match("ENTRY\s+(M\d{5})\s+.*",mo[line]):
            Mid = re.match("ENTRY\s+(M\d{5})\s+.*",mo[line]).group(1)
            DEF[Mid] = []
        elif re.match("ENTRY\s+[^M]\w.*",mo[line]):
            while not re.search("///",mo[line]):  
                line += 1
        elif re.match("NAME\s+(\S+.*)", mo[line]):
            Name = re.match("NAME\s+(\S+.*)", mo[line]).group(1)
            Name_Module[Mid] = Mid + " " + Name
        elif re.match("DEFINITION\s+(\S+.*)",mo[line]):
            while not re.search("ORTHOLOGY", mo[line]):
                if re.match("DEFINITION\s+(\S+.*)",mo[line]):
                    DEF[Mid].append(re.match("DEFINITION\s+(\S+.*)",mo[line]).group(1))
                else:
                    DEF[Mid].append(mo[line])
                    
                line += 1
        line += 1    
    return(DEF, Name_Module)

import pandas as pd

def get_dict(gene2module):
    KOdic = {}
    df = pd.read_table(gene2module,sep="\t", header=None)
    rows = df.shape[0]
    for cnt in range(0, rows):
        KOdic[df.iloc[cnt,1]] = 1
    return(KOdic)

import numpy as np

def de_plus(string):
    plus = re.findall("\d(?:\+\d)+", string)
    for plus_item in plus:             # simplify strucutre(+)
        if eval(plus_item) / ((len(plus_item)+1)/2) >= 0.75: #多个亚基的类型，需要出来的比例
            value = 1
        else:
            value = 0
        string = string.replace(plus_item,str(value),1)
    return(string)

def de_under(item):
    under = re.findall("\d(?:_\d)+", item)
    for under_item in under:
        under_list = under_item.split("_")
        under_list = list(map(int, under_list))
       
        under_array = np.array(under_list)
        if under_array.all() > 0:
            value = 1
        else:
            value = 0
        item = item.replace(under_item,str(value),1)
    return(item)

def de_comma(item):
    comma = re.findall("\d(?:,\d)+", item)
    for comma_item in comma:         # simplify alternative step(,)
        domma_list = comma_item.split(",")
        domma_list = list(map(int, domma_list))
        domma_array = np.array(domma_list)
        if domma_array.any() > 0:
            value = 1
        else:
            value = 0
        item = item.replace(comma_item, str(value),1)
    return(item)

def de_bracket(lists,KOdict):
    matrix = []
    total_progress = len(lists)
    for cnt in range(0, total_progress):   #replace with 0,1 for the absence and prsence of each ko
        #lists[cnt] = re.sub(r'^\(',"",lists[cnt])
        #lists[cnt] = re.sub(r'\)$',"",lists[cnt])
        koLists = re.findall("K\d{5}",lists[cnt])

        for ko in koLists:
            if ko in KOdict:
                value = 1
            else:
                value = 0
            lists[cnt] = lists[cnt].replace(ko,str(value))
        
        lists[cnt] = lists[cnt].replace("--","1")
        while re.search("\([^(]+\)",lists[cnt]): # get the elements bracketed
            enbraced = re.findall("\([^(]+?\)",lists[cnt])
            for item in enbraced:
                oriItem = item
                #print(f"1:{item}")
                # deal with accessory situation
                item = item.replace("-1","")
                item = item.replace("-0","")
                #print(f"2:{item}")
                # deal with plus situation
                item = de_plus(item)
                #print(f"3:{item}")
                item = de_under(item)
                #print(f"4:{item}")
                item = de_comma(item)
                #print(f"5:{item}")
                # deal with the alternative situation
                item = item.strip("(")
                item = item.strip(")") 
                lists[cnt] = lists[cnt].replace(oriItem, str(item), 1)
        matrix.append(lists[cnt])
    return(matrix)

def inter(a,b):
    return(list(set(a)&set(b)))

def simplify(matrix,completeness): 
    result = []
    second_list = []
    for item in matrix: 
        if re.search("M\d{5}", item):
            for Module in re.findall("M\d{5}", item): 
                if Module in completeness:
                    if completeness[Module] == 1:
                        value = "1"
                    else:
                        value = "0"
                    item = item.replace(Module, value)
                else:
                    return(0)
                    break
   
        while inter(["+","-","_",","], item):
            item = item.replace("-1","")
            item = item.replace("-0","")
            item = de_plus(item)
            item = de_under(item)
            item = de_comma(item)
        result.append(item)
    return(result)

def completeness_cal(lists):
    lists = list(map(lambda x:'1' if x=="" else x, lists))
    lists = list(map(int, lists))
    return(sum(lists)/len(lists))

with open(module,'r') as MD:
    moduleL = MD.readlines()
Mo_list = get_MO(moduleL)[0]
ModuleName = get_MO(moduleL)[1]


def complete_MAG(Mo_list,KO_dict):
    complete_res = {}
    for module in Mo_list:
        values = []
        for mo_item in Mo_list[module]:          
            first = ko_hierachy_split(mo_item)
            #print(first)
            debraket = de_bracket(first,KO_dict)
            #print(debraket)
            trans = simplify(debraket,complete_res)
            #print(trans)
            if trans:
                completeness = completeness_cal(trans)
                values.append(completeness)
        if len(values):
            res = max(values)
        else:
            res=""
        complete_res[module] = res
    return(complete_res)
df = pd.DataFrame()
import glob
pathpattern = f"{directory}/*{suffix}"
source = glob.glob(pathpattern)
for item in source:
    sample = re.search(f"\S+/(\w\S+?){suffix}",item).group(1)
    KO_dict = get_dict(item)
    completeness = complete_MAG(Mo_list,KO_dict)
    #print(sample)
    #print(completeness)
    df[sample] = pd.Series(completeness)
df.index = df.index.map(lambda x:ModuleName[x])
df.to_csv(f"{output}.tsv",sep="\t")
df.drop(df[df.iloc[:,0] == ""].index, inplace = True)
heatDF = df[df.sum(1)>=1]

import seaborn as sns
import matplotlib.pyplot as plt
fontsize = 9
plt.rcParams['ytick.labelsize'] = fontsize

fontsize_pt = plt.rcParams['ytick.labelsize']

dpi = 72.27
matrix_width_pt = fontsize_pt * (heatDF.shape[1]+80)
matrix_width_in = float(matrix_width_pt) / dpi
left_margin = 0.05
right_margin = 0.15
figure_width = matrix_width_in / (1 - left_margin - right_margin)
matrix_height_pt = fontsize_pt * heatDF.shape[0]
matrix_height_in = float(matrix_height_pt) / dpi
top_margin = 0.06  # in percentage of the figure height
bottom_margin = 0.06 # in percentage of the figure height
figure_height = matrix_height_in / (1 - top_margin - bottom_margin)


maps = sns.clustermap(data=heatDF,
               row_cluster=True, #行方向不聚类
               col_cluster=True, #列方向聚类
               cmap = "vlag",
      
               xticklabels = True,
               yticklabels = True,
               figsize = (figure_width,(figure_height)),
              )

maps.ax_heatmap.set_xticklabels(maps.ax_heatmap.get_xmajorticklabels(), fontsize = fontsize)
maps.ax_heatmap.set_yticklabels(maps.ax_heatmap.get_ymajorticklabels(), fontsize = fontsize)
maps.savefig(f"{output}.pdf")


