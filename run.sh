# Predicting protein alternative conformations from coevolution

# Run as:
# <codename> <input_file> 

#!/bin/bash
begin=$(date +"%s")


input_file=$1

while read line;
do
       	if [[ $line == *"sourceCodePath"* ]]
       	then
       		eval "$line";
       	fi

        if [[ $line == *"confoldPath"* ]]
        then
            eval "$line";
        fi
       	
        if [[ $line == *"output_path"* ]]
        then
            eval "$line";
            break;
       	fi
done < ${input_file}


echo "Stage 0: Generate the best default 3D structure and the optimal number of contacts"
#python $sourceCodePath/generateDefault.py ${input_file}

echo "Stage 1: Generate a list of candidate clusters"
folder=stage1
mkdir ${folder}
python $sourceCodePath/generateCandidates.py ${input_file}

sge=$(ls ${output_path}/stage1/sge.confold.* 2> /dev/null | wc -l)
while [ "$sge" != "0" ]
do
	sleep 2m
	sge=$(ls ${output_path}/stage1/sge.confold.* 2> /dev/null | wc -l)
done

echo "Evaluate stage 1 structures"
python $sourceCodePath/evaluation.py ${input_file} ${output_path}/stage1/cluster_list_stage1.npy 1

echo "Stage 2: Grow candidate clusters"
folder=stage2
mkdir ${folder}
python $sourceCodePath/growCandidates.py ${input_file} ${output_path}/stage1/cluster_list_stage1.npy ${output_path}/stage1/reward_info_stage1.csv

sge=$(ls ${output_path}/stage2/sge.confold.* 2> /dev/null | wc -l)
while [ "$sge" != "0" ]
do
	sleep 2m
	sge=$(ls ${output_path}/stage2/sge.confold.* 2> /dev/null | wc -l)
done

echo "Evaluate stage 2 structures"
python $sourceCodePath/evaluation.py ${input_file} ${output_path}/stage2/cluster_list_stage2.npy 2

echo "Stage 3: Select high scoring clusters and predicted structures"
folder=stage3
mkdir ${folder}
python $sourceCodePath/selection.py ${input_file} ${output_path}/stage1/cluster_list_stage1.npy ${output_path}/stage2/cluster_list_stage2.npy ${output_path}/stage1/reward_info_stage1.csv ${output_path}/stage2/reward_info_stage2.csv

if [[ -f $output_path/*log ]]; then
  rm $output_path/*log;
fi

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "Runtime: $(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds" > $output_path/runtime

if [ -f $structure_A ] && [ -f $structure_B ]; then
  if [ -f $output_path/stage3/reward_info_final.csv ]; then
    echo "Compare high scoring structures with reference structures";
    python $sourceCodePath/comparison.py ${input_file};
    echo "All jobs are finished";
  else
    echo "No high scoring structures to be compared with, all jobs are finished";
  fi
else
  echo "No reference structures to be compared with, all jobs are finished";
fi


