"""
This file selects candidates from stage 1 based on zscore_reward and generates a new set of clusters based on all possible combinations among candidates and the predefined cluster size upper bound 

Run as:
python <codename> <input_file> <cluster_list_file_1> <reward_info_file_1>

Returns:
stage2: a folder containing the best 3D model for each cluster and candidates.csv (selected candidates for growing based on zscore_reward)
"""
import os, sys
from lib import *

# Read input file
f_log=open("growCandidates.log","w")
if len(sys.argv)==4:
    input_file = sys.argv[1]
    cluster_list_file = sys.argv[2]
    reward_info_file = sys.argv[3]
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
dir_out=output_path+"/stage2"

# Read the optimal contact set (-default.rr)
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

# Read cluster_list file
cluster_list=list(np.load(cluster_list_file))

# Read reward_info file
reward_info=pd.read_csv(reward_info_file)

# Grow clusters by self-joining 
cluster_list_new=[]
candidates=reward_info[reward_info["zscore_reward"]>=zscore_thres]
candidates.to_csv(dir_out+"/candidates.csv", index=None, header=True)
candidates_index=candidates["index"].values

for k in range(2, len(candidates_index)+1):
	for combo in combinations(candidates_index, k):
		cluster_new=[]
		for i in range (len(combo)):
			cluster_new=list(set(cluster_new+cluster_list[combo[i]-1]))
		cluster_new=sorted(cluster_new)
		if cluster_new not in cluster_list and cluster_new not in cluster_list_new:
			#print (k, len(cluster_new), cluster_new)
			if float(len(cluster_new))/float(len(contact_default_info))<=bound[2]:
				cluster_list_new.append(cluster_new)
				#print (cluster_new)		
f_log.write("Stage 2, total clusters: "+str(len(cluster_list_new))+"\n")
# Save cluster_list
np.save(dir_out+"/cluster_list_stage2", cluster_list_new)

# Predict new structures for each cluster
for i in range(len(cluster_list_new)):
	cluster=cluster_list_new[i]
	contact_new=dir_out+"/"+id+"-cluster"+str(i+1)+".rr"
	contact_num=update_contact_file(contact_default, cluster, contact_new)
	output_path_new=dir_out+"/cluster"+str(i+1)
	fold(confoldPath, contact_new, secondary_structure_file, output_path_new)
f_log.write("Submit all jobs\n")
f_log.close()
