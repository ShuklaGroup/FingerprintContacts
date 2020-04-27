"""
This file selects high scoring structures from stage 1 and stage 2 based on reward

Run as:
python <codename> <input_file> <cluster_list_file_1> <cluster_list_file_2> <reward_info_file_1> <reward_info_file_2>

Returns:
stage3: a folder containing high scoring structures and corresponding information
"""
import os, sys 
from lib import *

# Read input file
f_log=open("final_results.outcfg","w")
if len(sys.argv)==6:
    input_file = sys.argv[1]
    cluster_list_file_1 = sys.argv[2]
    cluster_list_file_2 = sys.argv[3]
    reward_info_file_1 = sys.argv[4]
    reward_info_file_2 = sys.argv[5]
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
dir_out=output_path+"/stage3"

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
f_default.close()
contact_default_info=pd.DataFrame(contact_default_info, columns = ["res_i", "res_j", "p"])

# Read cluster_list file
cluster_list_1=np.load(cluster_list_file_1)
cluster_list_2=np.load(cluster_list_file_2)

# Read reward_info file
reward_info_1=pd.read_csv(reward_info_file_1)
reward_info_2=pd.read_csv(reward_info_file_2)

# Select high scoring clusters and structures
cluster_list_final=[]
cluster_final_1=reward_info_1[reward_info_1["reward"]>=reward_thres]
cluster_final_2=reward_info_2[reward_info_2["reward"]>=reward_thres]
reward_info_final=pd.concat([cluster_final_1, cluster_final_2])
reward_info_final["index"]=range(1, len(reward_info_final)+1)
reward_info_final=reward_info_final.drop("zscore_reward", axis=1)
if len(reward_info_final) !=0:
	reward_info_final.to_csv(dir_out+"/reward_info_final.csv", index=None, header=True)

count=0
for i in cluster_final_1["index"].values:
	count+=1
	cluster=cluster_list_1[i-1]
	cluster_list_final.append(cluster)
	cmd="cp -r "+output_path+"/stage1/cluster"+str(i)+" "+dir_out+"/cluster"+str(count)
	os.system(cmd)

for i in cluster_final_2["index"].values:
	count+=1
	cluster=cluster_list_2[i-1]
	cluster_list_final.append(cluster)
	cmd="cp -r "+output_path+"/stage2/cluster"+str(i)+" "+dir_out+"/cluster"+str(count)
	os.system(cmd)

# Save final cluster_list
if len(cluster_list_final) !=0:
	np.save(dir_out+"/cluster_list_stage3", cluster_list_final)

# Save all the residue pair info in the final cluster_list
for i in range(len(cluster_list_final)):
	contact_final_info=pd.DataFrame(columns = ["res_i", "res_j", "p"])
	cluster=cluster_list_final[i]
	for index in cluster:
		contact_final_info=contact_final_info.append(contact_default_info.iloc[index])
	contact_final_info.insert(loc=0, column="index", value=contact_final_info.index.values+1)
	contact_final_info.to_csv(dir_out+"/contact_final_info_"+str(i+1)+".csv", index=None, header=True)

f_log.write("Results Summary\n")
f_log.write("Stage 0: Default structure information... \n")
f_contact=open(contact_file, "r")
seq=f_contact.readline().strip().split(" ")[0]
f_contact.close()
L=len(seq)
f=open(dir_default+"/clustering/centroids.txt","r")
line=f.readline()
f.close()
best_model=os.path.basename(line.strip().split(" ")[-1])
optimal_num=int(round(float(best_model.split("-")[1][:-1])*L))
f_log.write("Sequence_length: "+str(L)+"\n")
f_log.write("Best_model: "+best_model+"\n")
f_log.write("Optimal_contact_number: "+str(optimal_num)+"\n")
f_log.write("Optimal_contact_set: "+dir_default+"/"+id+"-default.rr\n")
f_log.write("\n")
f_log.write("Stage 1, total clusters: "+str(len(cluster_list_1))+"\n")
f_log.write("Stage 2, total clusters: "+str(len(cluster_list_2))+"\n")
f_log.write("\n")
f_log.write("Stage 1, high scoring clusters: "+str(len(cluster_final_1))+"\n")
f_log.write("Stage 2, high scoring clusters: "+str(len(cluster_final_2))+"\n")
f_log.write("\n")
f_log.write("Total high scoring clusters: "+str(len(cluster_final_1)+len(cluster_final_2))+"\n")
if len(reward_info_final) != 0:
	f_log.write("Maximum reward: "+str(max(reward_info_final["reward"].values))+"\n")
	f_log.write("Minimum tmscore: "+str(min(reward_info_final["tmscore"].values))+"\n")
f_log.close()

print("Final results:")
if len(reward_info_final)!=0:
	print("Total high scoring clusters: "+str(len(cluster_final_1)+len(cluster_final_2)))
	print("Maximum reward: "+str(max(reward_info_final["reward"].values)))
	print("Minimum tmscore: "+str(min(reward_info_final["tmscore"].values)))
else:
	print ("No high scoring clusters found")

