#!/bin/bash
#$ -cwd -V
#$ -N csat2
#$ -l h_rt=72:00:00
##$ -q all.q@comp702.bullx

ulimit -s unlimited

export RASPA_DIR=/home/lchen/EXEs/raspa-2.0-cf-v1.9.5_may

FILE_LIST=LIST_POCs
NUM_FILES=42
START_FILE=1

for ((INDEX = $START_FILE; INDEX <= $NUM_FILES; INDEX++))
do

SOURCE=../../GetCoords/MakeMolFile_P1

### Get the crystal filename
ARG1=$(sed -n "$INDEX p" $FILE_LIST | awk '{print $1}')

### Get mol file
cp ${SOURCE}/mol_files/${ARG1}.mol ./

### Get cell info
cell_a=$(tail -4 ${ARG1}.mol | head -1 | awk '{print $1}')
cell_b=$(tail -4 ${ARG1}.mol | head -1 | awk '{print $2}')
cell_c=$(tail -4 ${ARG1}.mol | head -1 | awk '{print $3}')

cell_alpha=$(tail -3 ${ARG1}.mol | head -1 | awk '{print $1}')
cell_beta=$(tail -3 ${ARG1}.mol | head -1 | awk '{print $2}')
cell_gamma=$(tail -3 ${ARG1}.mol | head -1 | awk '{print $3}')

# Density in kg/m3
density=$(cat ${SOURCE}/Output/System_0/output_${ARG1}_*.data | grep "Framework Density" | tail -1 | awk '{print $3}')
cell_density=$(printf "%0.6f\n" $(bc <<< "$density * 0.001"))

echo $ARG1"   "${cell_density} >> results_density.dat

### Make input file
cp INPUT INPUT_$ARG1

sed -i "s/job_coords/${ARG1}.mol/g" ./INPUT_$ARG1

sed -i "s/cell_a/${cell_a}/g" ./INPUT_$ARG1
sed -i "s/cell_b/${cell_b}/g" ./INPUT_$ARG1
sed -i "s/cell_c/${cell_c}/g" ./INPUT_$ARG1

sed -i "s/CELL_ALPHA/${cell_alpha}/g" ./INPUT_$ARG1
sed -i "s/CELL_BETA/${cell_beta}/g" ./INPUT_$ARG1
sed -i "s/CELL_GAMMA/${cell_gamma}/g" ./INPUT_$ARG1

sed -i "s/Cell_Density/${cell_density}/g" ./INPUT_$ARG1

### Run

../build/nonorthoSA.exe < INPUT_$ARG1 > OUTPUT_$ARG1

done
