"""
This file defines useful functions
"""

import os,sys
import numpy as np
import pandas as pd 
from sklearn.cluster import AgglomerativeClustering
from scipy import stats
from itertools import combinations

def generateConfoldscript(confoldPath, contact_new, secondary_structure_file, output_path_new):
	""" 
	Generate bash script for CONFOLD2

    Parameters
    ----------
    confoldPath: str
        the path of confold2 software
    contact_new: str
        the updated contact filename
    secondary_structure_file: str
        the filename of the secondary structure provided in the input file (*.inp)
    output_path_new: str
        the output path for CONFOLD2 results

    Returns
    ----------
    sge_file: str
        the submission script for CONFOLD2
	"""

	jobname=os.path.basename(output_path_new)
	sge_file=os.path.dirname(output_path_new)+"/sge.confold."+jobname
	f_sge=open(sge_file,"w")
	f_sge.write("#$ -S /bin/bash"+"\n"
+"#$ -q analysis.q"+"\n"
+"#$ -pe orte 1"+"\n"
+"#$ -o "+os.path.dirname(output_path_new)+"/confold_"+jobname+".log"+"\n"
+"#$ -e "+os.path.dirname(output_path_new)+"/confold_"+jobname+".err"+"\n"
+"#$ -N confold_"+jobname+"\n"
+""+"\n"
+"cd "+os.path.dirname(output_path_new)+"\n"
+""+"\n"
+confoldPath+"/core.pl -mcount 20 -rr "+contact_new+" -ss "+secondary_structure_file+" -selectrr all -o "+output_path_new+"\n"
+""+"\n"
+"cd "+output_path_new+"\n"
+"cp "+output_path_new+"/stage2/*model1*pdb "+output_path_new+"/\n"
+"tar -czf stage1n2.tar.gz stage1 stage2 input"+"\n"
+"rm -r stage1 stage2 input"+"\n"
+"mv "+sge_file+" "+output_path_new+"/\n"
)
	f_sge.close()
	return sge_file

def submitConfoldscript(sge_file):
	""" 
	Submit CONFOLD2 job

    Parameters
    ----------
    sge_file: str
        the submission script for CONFOLD2
	"""

	cmd="chmod +x "+sge_file
	os.system(cmd)
	cmd="qsub "+sge_file
	os.system(cmd)

def update_contact_file(contact_default, cluster, contact_new):
	"""
	Compute the new contact file based on the state vector

    Parameters
    ----------
    contact_default: str
        the filename of the default contact set genearted by generateDefault.py
    cluster: numpy.ndarray (int)
        the cluster vector with the shape N*2, where N is the number of to be removed contacts
    contact_new: str
        the output filename for the new contact set

    Returns
    ----------
    contact_new: str
        the updated contact file serves as input for CONFOLD2   
    contact_num: int
        the updated number of contacts
    """

	f_new=open(contact_new, "w")
	f_default=open(contact_default, "r")
	indx=0
	contact_num=0
	for line in f_default:
		if indx==0:
			# copy the first line: sequence info
			f_new.write(line)
			
		elif indx-1 not in cluster:
			# copy the active contact info
			f_new.write(line)
			contact_num+=1
		indx+=1
	f_default.close()
	f_new.close()
	return contact_num

def fold(confoldPath, contact_new, secondary_structure_file, output_path_new):
	""" 
	Compute the best 3D structure using CONFOLD2

	Parameters
	----------
    confoldPath: str
	    the path of confold2 software
    contact_new: str
	    the updated contact filename
    secondary_structure_file: str
	    the filename of the secondary structure provided in the input file (*.inp)
    output_path_new: str
	    the output path for CONFOLD2 results

	Returns
    ----------
	contact_new: str
	    the updated contact file serves as input for CONFOLD2   
	contact_num: int
	    the updated number of contacts
	"""

	sge_file=generateConfoldscript(confoldPath, contact_new, secondary_structure_file, output_path_new)
	submitConfoldscript(sge_file)

def tm_score(confoldPath, structure_new, structure_default):
	""" 
	Compute the TM score and RMSD between the new and default structure

	Parameters
	----------
    confoldPath: str
	    the path of confold2 software
    structure_new: str
	    the predicted structure based on the updated contact set
    structure_default: str
	    the predicted structure based on the default contact set

	Returns
    ----------
	tmscore: float (range: (0, 1]) 
	    the TM score between the new and default structure  
	rmsd: float
	    the RMSD between the new and default structure
	"""

	output=structure_new.split(".")[0]+"-evaluation.txt"
	cmd=confoldPath+"/third-party-programs/TMscore "+structure_new+" "+structure_default+" &> "+output
	os.system(cmd)
	# read TM output file
	f=open(output,"r")
	for line in f:
		columns=line.strip().split()
		if len(columns)>3 and columns[0]=="RMSD" and columns[3]=="common":
			rmsd=float(columns[-1])
		if len(columns)>1 and columns[0]=="TM-score" and columns[1]=="=":
			tmscore=float(columns[2])
	f.close()
	return tmscore, rmsd

def get_atom_position(structure_new, residue_id):
	""" 
	Compute the xyz position of CB atom (CA atom for GLY) for the given residue 

	Parameters
	----------
    structure_new: str
	    the predicted structure based on the updated contact set
    residue_id: int
	    the 1-indexed residue index

	Returns
    ----------
	xyz: numpy.ndarray (float)
	    the atom position vector with the shape (3, )
	"""

	f=open(structure_new, "r")
	for line in f:
		columns=line.strip().split()
		if len(columns)==11 and columns[0]=="ATOM" and columns[4]==str(residue_id):
			if columns[3]=="GLY" and columns[2]=="CA":
				# For GLY, return the position vector of CA atom
				xyz=[float(columns[5]), float(columns[6]), float(columns[7])]
			elif columns[2]=="CB":
				xyz=[float(columns[5]), float(columns[6]), float(columns[7])]

		if len(columns)==12 and columns[0]=="ATOM" and columns[4]=="A" and columns[5]==str(residue_id):
			if columns[3]=="GLY" and columns[2]=="CA":
				# For GLY, return the position vector of CA atom
				xyz=[float(columns[6]), float(columns[7]), float(columns[8])]
			elif columns[2]=="CB":
				xyz=[float(columns[6]), float(columns[7]), float(columns[8])]
	f.close()
	xyz=np.array(xyz)
	return xyz

def contact_center(structure_new, contact_pair):
	""" 
	Compute the center of two vectors defining xyz positions

	Parameters
	----------
    structure_new: str
	    the predicted structure based on the updated contact set
    contact_pair: list (int)
        the residue index vector with the shape (2, )

	Returns
    ----------
	center: numpy.ndarray (float)
	    the center of two atom position vectors with the shape (3, )
	"""
	r_1=get_atom_position(structure_new, contact_pair[0])
	r_2=get_atom_position(structure_new, contact_pair[1])
	center=0.5*(r_1+r_2)
	return center 

def contact_filter(structure_new, contact_pair):
	""" 
	Check whether a contact pair is satified in the given structure
	Definition of a satisfied contact pair: CB-CB distance <= 8 Angstorms (CA for Glycine)

	Parameters
	----------
    structure_new: str
	    the predicted structure based on the updated contact set
    contact_pair: list (int)
        the residue index vector with the shape (2, )

	Returns
    ----------
	contact_fil: bool 
	    1: contact pair is satisfied in the given structure
	    0: contact pair is not satisfied in the given structure
	"""

	r_1=get_atom_position(structure_new, contact_pair[0])
	r_2=get_atom_position(structure_new, contact_pair[1])
	dist=np.sqrt(np.sum((r_1-r_2)**2))
	if dist <= 8:
		contact_fil=1
	else:
		contact_fil=0
	return contact_fil 


def contact_satisfaction(structure_new, contact_default_info, Q_max, optimal_num):
	""" 
	Compute the normalized contact satisfaction score for the given structure
	Definition of the normalized contact satisfaction score: Q_s=Q/Q_max
	Q: the sum of the weights of satisfied contacts in the given structure
	Q_max (the maximum contact satisfication): the sum of all contact weights

	Parameters
	----------
    structure_new: str
	    the predicted structure based on the updated contact set
    contact_default_info: pandas.DataFrame (three columns: res_i, res_j, p)
        the information of the default contacts containing the resiude index and weight of each contact
    Q_max: float
        the maximum contact satisfication (the sum of all contact weights)
    optimal_num: int
        the total number of the default contacts

	Returns
    ----------
	Q_s: float (range: (0, 1]) 
	    the normalized contact satisfaction score for the given structure
	"""

	Q=0
	for i in range(optimal_num):
		contact_pair=[int(contact_default_info.iloc[i]["res_i"]), int(contact_default_info.iloc[i]["res_j"])]
		if contact_filter(structure_new, contact_pair)==1:
			Q+=float(contact_default_info.iloc[i]["p"])
	Q_s=Q/Q_max
	return Q_s

def cal_reward(confoldPath, structure_new, structure_default, contact_default_info, Q_max, optimal_num, tm_thres, Q_s_thres):
	""" 
	Compute the reward function
	Reward=1/tmscore-1
	Q_s: the normalized contact satisfaction score for the new structure
	tmscore: the TM score between the new and default structure  

	Parameters
	----------
	confoldPath: str
	    the path of confold2 software
    structure_new: str
	    the predicted structure based on the updated contact set
    structure_default: str
	    the predicted structure based on the default contact set
    contact_default_info: pandas.DataFrame (three columns: res_i, res_j, p)
        the information of the default contacts containing the resiude index and weight of each contact
    Q_max: float
        the maximum contact satisfication (the sum of all contact weights)
    optimal_num: int
        the total number of the default contacts
    tm_thres: float
    	the threshold value of tmscore 
    Q_s_thres: float
    	the threshold value of Q_s (the normalized contact satisfaction score for the given structure)

	Returns
    ----------
	reward: float
	    the reward of the new structure
	"""

	tmscore, rmsd=tm_score(confoldPath, structure_new, structure_default)
	Q_s=contact_satisfaction(structure_new, contact_default_info, Q_max, optimal_num)
	if tmscore>=tm_thres and Q_s>=Q_s_thres: 
		reward=1.0/tmscore-1
	else:
		reward=0
	return reward, Q_s, tmscore 

