"""
Given contact file and secondary structure file, this file generates the best 3D model using CONFOLD2 and outputs the optimal set of top contacts 

Run as:
python <codename> <input_file> 

Returns:
1. -model1.pdb: the best 3D model based on the input contact set
2. -default.rr: the optimal set of contacts that generates the best 3D model 
"""

import os, sys

# Read input file
f_log=open("generateDefault.log","w")
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
dir_out=output_path+"/stage0" 

# Read sequence info from input contact file
f_contact=open(contact_file, "r")
seq=f_contact.readline().strip().split(" ")[0]
f_contact.close()
L=len(seq)

f_log.write("Start CONFOLD jobs...\n")
cmd=confoldPath+"/confold2-main.pl -rr "+contact_file+" -ss "+secondary_structure_file+" -out "+dir_out
os.system(cmd)
f_log.write("\n")

# Read CONFOLD result
f_log.write("Read CONFOLD results...\n")
f=open(dir_out+"/clustering/centroids.txt","r")
line=f.readline()
f.close()
best_model=os.path.basename(line.strip().split(" ")[-1])
optimal_num=int(round(float(best_model.split("-")[1][:-1])*L))
f_log.write("Sequence_length: "+str(L)+"\n")
f_log.write("Best_model: "+best_model+"\n")
f_log.write("Optimal_contact_number: "+str(optimal_num)+"\n")
f_log.write("\n")

# Output the best model (.pdb) and the optimal contact set file (-default.rr)
f_log.write("Extract the optimal contact set: "+dir_out+"/"+id+"-default.rr\n")

f_default=open(dir_out+"/"+id+"-default.rr", "w")
f_contact=open(contact_file, "r")
indx = 0
for line in f_contact:
	if indx < optimal_num+1:
		f_default.write(line)
		indx += 1
f_contact.close()
f_default.close()

f_log.close()




