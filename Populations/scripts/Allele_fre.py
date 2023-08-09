#!/usr/bin/env python
import pandas as pd
import sys
import numpy as np
refLen = float(sys.argv[1])

data = pd.read_table(sys.argv[2].strip(),sep="\t")
col = data.shape[1]
row = data.shape[0]

def count_above_zero(AD):
    info = AD.split(",")
    if np.sum(list(map(lambda x:float(x)>0, info))) > 1:
        return(True)
    else:
        return(False)

den = {}
newDataFrame = pd.DataFrame()
title = data["CHROM"] + data["POS"].astype(str)
for item in range(3,col,2):
    ref = pd.to_numeric(data.iloc[:,item].apply(lambda x:x.split(",")[0])) #get the count value of the ref
    sample = list(data)[item].strip('.AD')             #get the sample name
    snp = len(data.iloc[:,item][data.iloc[:,item].apply(count_above_zero)])
    density = snp / refLen
    den[sample] = [density]
    newDataFrame[sample] = ref / data.iloc[:,item+1].astype(float)
print(den)
denDF = pd.DataFrame.from_dict(den)
newDataFrame.insert(0, "POS", title, True)
newDataFrame.to_csv(sys.argv[3], index=None)
denDF.to_csv(sys.argv[4])
