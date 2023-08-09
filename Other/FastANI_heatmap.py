#!/nfs_genome/anaconda/envs/rnaseq/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import sys
import getopt

def main(argv):
    input_file = ''
    output_file = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=", "ofile="])
    except getopt.GetoptError:
        print('heatmap.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('heatmap.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_file = arg
        elif opt in ("-o", "--ofile"):
            output_file = arg
    return input_file, output_file

file_name, output_file = main(sys.argv[1:])

df = pd.read_csv(file_name, sep="\t", header=None, names=['Query', 'Reference', 'ANI', 'Mapped','Total'])

reciprocal_df = df.copy()
reciprocal_df['Query'], reciprocal_df['Reference'] = df['Reference'], df['Query']
df = pd.concat([df, reciprocal_df])

df = df.drop_duplicates(subset=['Query', 'Reference'])
pivot_df = df.pivot(index='Query', columns='Reference', values='ANI')

pivot_df = pivot_df.fillna(0)

# Calculate the distances and linkage
distances = pdist(pivot_df.values)
Z = linkage(distances, method='average')

# Create a new figure
fig = plt.figure(figsize=(10,10))

# Create a subplot for the dendrogram
ax1 = fig.add_axes([0.05,0.1,0.2,0.6])
dendrogram(Z, orientation='right', ax=ax1)
ax1.axis('off')

# Mirror the dendrogram
ax1.invert_xaxis()
ax1.invert_yaxis()

# Reorder the DataFrame according to the dendrogram
idx = dendrogram(Z, no_plot=True)['leaves']
pivot_df = pivot_df.iloc[idx, idx]

# Create a subplot for the heatmap
ax2 = fig.add_axes([0.3,0.1,0.6,0.6])
sns.heatmap(pivot_df, cmap='RdBu_r', ax=ax2)

# Move the labels to the right
ax2.yaxis.tick_right()
plt.yticks(rotation=0)

# Save the plot to a PDF file
plt.savefig(output_file)

