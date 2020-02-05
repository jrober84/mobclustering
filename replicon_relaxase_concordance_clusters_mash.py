from sklearn.metrics.cluster import adjusted_rand_score
from skbio.diversity.alpha import shannon
from scipy.stats import entropy
from statistics import mean
from statistics import median
from statistics import stdev
from math import e
from math import log
import pandas as pd
import sys

input_data = pd.read_csv('/Users/jrobertson/PycharmProjects/__mobcluster/__2019-12-results/2019-12-10-mash-clustering-average.integrated.txt',header=0,sep="\t")

thresholds = [
'clust_0.1',
'clust_0.09',
'clust_0.08',
'clust_0.07',
'clust_0.06',
'clust_0.05',
'clust_0.04',
'clust_0.03',
'clust_0.02',
'clust_0.01',
'clust_0'
]

cdf = input_data[input_data['rep_type(s)'] != '-']


rel_codes = cdf['relaxase_type(s)'].value_counts().index.to_list()
rep_codes = cdf['rep_type(s)'].value_counts().index.to_list()

print(">>>>>>>Replicons<<<<<<<<<<<<")
overall_entropy = {}
overal_cluster_size = {}
for thresh in thresholds:
    cluster_codes = cdf[thresh].value_counts().index.to_list()



    # Relaxase first
    rel_shannon_list = []
    cluster_sizes = cdf[thresh].value_counts().to_list()
    cluster_members = []
    for code in cluster_codes:
        sdf = cdf[cdf[thresh] == code]


        if len(sdf) == 0:
            continue
        vc = sdf['rep_type(s)'].value_counts()
        #print(vc)
        vcs = vc.sum()
        vc = vc.apply(lambda x: float(x)/vcs)
        #print(vc)
        max = log(len(vc))
        #print(f"num_cat:{len(vc)} max:{max}")
        cluster_members.append(len(vc))
        if (len(vc) == 1):
            rel_shannon_list.append(0.0)
        else:
            rel_shannon_list.append(entropy(vc.to_list())/max)
        # rel_shannon_list.append(shannon(vc.to_list(),e))


    overall_entropy[thresh] = rel_shannon_list
    overal_cluster_size[thresh] = cluster_members
    if len(rel_shannon_list) > 1:
        new_rel_shannon_list = []
        new_cluster_sizes = []
        new_cluster_members = []
        for i in range(0,len(rel_shannon_list)):
            for k in range(0,cluster_sizes[i]):
                new_rel_shannon_list.append(rel_shannon_list[i])
                new_cluster_sizes.append(cluster_sizes[i])
                new_cluster_members.append(cluster_members[i])



        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(thresh,mean(new_rel_shannon_list),stdev(new_rel_shannon_list),mean(new_cluster_sizes),stdev(new_cluster_sizes),mean(new_cluster_members),stdev(new_cluster_members)))
    else:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(thresh, rel_shannon_list[0], 0,cluster_sizes[0],0,cluster_members[0],0))





print(">>>>>>>Relaxases<<<<<<<<<<<<")

cdf = input_data[input_data['relaxase_type(s)'] != '-']
for thresh in thresholds:
    cluster_codes = cdf[thresh].value_counts().index.to_list()

    # Relaxase first
    rel_shannon_list = []
    cluster_sizes = cdf[thresh].value_counts().to_list()
    cluster_members = []
    for code in cluster_codes:
        sdf = cdf[cdf[thresh] == code]


        if len(sdf) == 0:
            continue
        vc = sdf['relaxase_type(s)'].value_counts()

        vcs = vc.sum()
        vc = vc.apply(lambda x: float(x)/vcs)

        max = log(len(vc))

        cluster_members.append(len(vc))
        if (len(vc) == 1):
            rel_shannon_list.append(0.0)
        else:
            rel_shannon_list.append(entropy(vc.to_list())/max)

    if len(rel_shannon_list) > 1:
        new_rel_shannon_list = []
        new_cluster_sizes = []
        new_cluster_members = []
        for i in range(0,len(rel_shannon_list)):
            for k in range(0,cluster_sizes[i]):
                new_rel_shannon_list.append(rel_shannon_list[i])
                new_cluster_sizes.append(cluster_sizes[i])
                new_cluster_members.append(cluster_members[i])

        print(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(thresh, mean(new_rel_shannon_list), stdev(new_rel_shannon_list),
                                                        mean(new_cluster_sizes), stdev(new_cluster_sizes),
                                                        mean(new_cluster_members), stdev(new_cluster_members)))
    else:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(thresh, rel_shannon_list[0], 0,cluster_sizes[0],0,cluster_members[0],0))
