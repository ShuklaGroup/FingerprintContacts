"""
This file generates reward information file for each cluster 

Run as:
python <codename> <input_file> <cluster_list_file> <stage>

Returns:
reward_info_stage*.csv: the .csv file with columns: "index", "contact_num", "reward", "Q_s", "tmscore", "zscore_reward"
"""
import os, sys
from lib import *

# Read input file
if len(sys.argv)==4:
    input_file = sys.argv[1]
    cluster_list_file = sys.argv[2]
    stage = sys.argv[3]
    f=open(input_file,"r")
    for line in f:
        exec (line)
    f.close()
else:
    print ("This script needs exactly "+str(len(sys.argv)-1)+" arguments, aborting")
    sys.exit()

id=os.path.basename(contact_file).split(".")[0]
dir_default=output_path+"/stage0" 
contact_default=dir_default+"/"+id+"-default.rr"
structure_default=dir_default+"/"+id+"-model1.pdb"
dir_out=output_path+"/stage"+stage

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

# Calculate the maximum contact satisfication: the sum of contact weights
Q_max=np.sum(contact_default_info["p"].values)

# Read cluster_list file
cluster_list=np.load(cluster_list_file)

# Calculate reward info
reward_info=pd.DataFrame(index=range(1, len(cluster_list)+1), columns=["index", "contact_num", "reward", "Q_s", "tmscore"])
for i in range(len(cluster_list)):
	output_path_new=dir_out+"/cluster"+str(i+1)
	structure_new=output_path_new+"/"+id+"-cluster"+str(i+1)+"_model1.pdb"
	if os.path.isfile(structure_new):
		cluster=cluster_list[i]
		contact_num=len(cluster)
		reward, Q_s, tmscore=cal_reward(confoldPath, structure_new, structure_default, contact_default_info, Q_max, optimal_num, tm_thres, Q_s_thres)
		reward_info.loc[i+1]=[i+1, contact_num, reward, Q_s, tmscore]

reward_info=reward_info.dropna()
reward_list=reward_info["reward"].values
if len(reward_list)>1:
	z_list=stats.zscore(reward_list)
	reward_info["zscore_reward"]=z_list
else:
	reward_info["zscore_reward"]=np.nan

reward_info.to_csv(dir_out+"/reward_info_stage"+stage+".csv", index=None, header=True)


