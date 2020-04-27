"""
Given clustering parameters, this file conducts Agglomerative Clustering and selects different clusters of contacts to be removed for the prediciton of new structures

Run as:
python <codename> <input_file> 

Returns:
stage1: a folder containing the best 3D model for each cluster 
"""
import os, sys
from lib import *

# Read input file
f_log=open("generateCandidates.log","w")
if len(sys.argv)==2:
    input_file = sys.argv[1]
    f=open(input_file,"r")
    for line in f:
        exec (line)
    f.close()
else:
    print ("This script needs exactly "+str(len(sys.argv)-1)+" arguments, aborting")
    f_log.write("This script needs exactly "+str(len(sys.argv)-1)+" arguments, aborting\n")
    sys.exit()

id=os.path.basename(contact_file).split(".")[0]
dir_default=output_path+"/stage0" 
contact_default=dir_default+"/"+id+"-default.rr"
structure_default=dir_default+"/"+id+"-model1.pdb"
dir_out=output_path+"/stage1"

# Read the optimal contact set (-default.rr)
f_log.write("Start reading the optimal contact set: ")
contact_default_info=[]
optimal_num=0
f_default=open(contact_default, "r")
for line in f_default:
	if optimal_num !=0:
		columns=line.strip().split(" ")
		contact_default_info.append([int(columns[0]), int(columns[1]), float(columns[-1])])
	optimal_num+=1
optimal_num-=1
f_log.write(str(optimal_num)+" contact pairs\n")
f_default.close()
contact_default_info=pd.DataFrame(contact_default_info, columns = ["res_i", "res_j", "p"])

# Calculate the maximum contact satisfication: the sum of contact weights
Q_max=np.sum(contact_default_info["p"].values)
f_log.write("Q_max: "+str(Q_max)+"\n")

# Calculate the center of each residue pair
contact_center_list=[]
for i in range(len(contact_default_info)):
	contact_pair=[int(contact_default_info.iloc[i]["res_i"]), int(contact_default_info.iloc[i]["res_j"])]
	center=list(contact_center(structure_default, contact_pair))
	contact_center_list.append(center)
contact_center_list=np.array(contact_center_list)

# Agglomerative Clustering
## Cluster list is 0-indexed
n_clusters_list=np.arange(n_clusters_range[0], n_clusters_range[1]+n_clusters_range[2], n_clusters_range[2])
cluster_list=[]
for n_clusters in n_clusters_list:
	clustering = AgglomerativeClustering(n_clusters=n_clusters).fit(contact_center_list)
	labels=clustering.labels_
	for i in range(n_clusters):
		group=[]
		for j in range(len(labels)):
			if labels[j] == i:
				group.append(j)
		group=sorted(group)
		## remove redudant clusters
		if bound[0]<=float(len(group))/float(len(contact_default_info))<=bound[1] and group not in cluster_list:
			cluster_list.append(group)
f_log.write("Stage 1, total clusters: "+str(len(cluster_list))+"\n")
# Save cluster_list
np.save(dir_out+"/cluster_list_stage1", cluster_list)

# Predict new structures for each cluster
for i in range(len(cluster_list)):
	cluster=cluster_list[i]
	contact_new=dir_out+"/"+id+"-cluster"+str(i+1)+".rr"
	contact_num=update_contact_file(contact_default, cluster, contact_new)
	output_path_new=dir_out+"/cluster"+str(i+1)
	fold(confoldPath, contact_new, secondary_structure_file, output_path_new)
f_log.write("Submit all jobs\n")
f_log.close()

















