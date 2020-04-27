"""
Given two reference structures (structure_A and structure_B), this file compares the predicted high scoring structures with reference structures 

Run as:
python <codename> <input_file> 

Returns:
1. structure_default_comparison.csv: the .csv file comparing default structure and reference structures with columns: "Q_s", "tmscore_A", "tmscore_B", "rmsd_A", "rmsd_B"
2. reward_info_final_comparison.csv: the .csv file with columns: "index", "contact_num", "reward", "Q_s", "tmscore", "tmscore_A", "tmscore_B"
"""

import os, sys
from lib import *
import glob

# Read input file
if len(sys.argv)==2:
    input_file = sys.argv[1]
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
Q_max=np.sum(contact_default_info["p"].values)
contacts_default=contact_default_info[["res_i", "res_j"]].values

cluster_list=output_path+"/stage3/cluster_list_stage3.npy"
cluster_list=np.load(cluster_list)

# Calcluate TM by comparing structure_default and reference structures
reward_default=pd.DataFrame(index=[1], columns=["Q_s", "tmscore_A", "tmscore_B", "rmsd_A", "rmsd_B"])
reward, Q_s, tmscore=cal_reward(confoldPath, structure_default, structure_default, contact_default_info, Q_max, optimal_num, tm_thres, Q_s_thres)
tmscore_A, rmsd_A=tm_score(confoldPath, structure_default, structure_A)
tmscore_B, rmsd_B=tm_score(confoldPath, structure_default, structure_B)
reward_default.loc[1]=[Q_s, tmscore_A, tmscore_B, rmsd_A, rmsd_B]
reward_default.to_csv(dir_out+"/structure_default_comparison.csv", index=None, header=True)

# Calcluate TM by comparing structure_new and reference structures
reward_info=pd.DataFrame(index=range(1, len(cluster_list)+1), columns=["index", "contact_num", "reward", "Q_s", "tmscore", "tmscore_A", "tmscore_B"])
for structure_new in sorted(glob.glob(dir_out+"/cluster*/*pdb")):
	index=int(structure_new.split("/")[-2][7:])
	cluster=cluster_list[index-1]
	contact_num=len(cluster)
	reward, Q_s, tmscore=cal_reward(confoldPath, structure_new, structure_default, contact_default_info, Q_max, optimal_num, tm_thres, Q_s_thres)
	tmscore_A, rmsd_A=tm_score(confoldPath, structure_new, structure_A)
	tmscore_B, rmsd_B=tm_score(confoldPath, structure_new, structure_B)
	reward_info.loc[index]=[index, contact_num, reward, Q_s, tmscore, tmscore_A, tmscore_B]
	
reward_info.to_csv(dir_out+"/reward_info_final_comparison.csv", index=None, header=True)







