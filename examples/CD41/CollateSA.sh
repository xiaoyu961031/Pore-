#!/bin/bash
##$ -cwd -V
##$ -N collate
##$ -l h_rt=72:00:00

ulimit -s unlimited

FILE_LIST=LIST_POCs
NUM_FILES=42
START_FILE=1

for ((INDEX = $START_FILE; INDEX <= $NUM_FILES; INDEX++))
do

### Get the crystal filename ###################
ARG1=$(sed -n "$INDEX p" $FILE_LIST | awk '{print $1}')

### Get & write results
SA_H_v=$(cat OUTPUT_$ARG1 | grep "(m^2/cm^3)" | head -1 | tail -1 | awk '{print $9}')
SA_C_v=$(cat OUTPUT_$ARG1 | grep "(m^2/cm^3)" | head -2 | tail -1 | awk '{print $9}')
SA_N_v=$(cat OUTPUT_$ARG1 | grep "(m^2/cm^3)" | head -3 | tail -1 | awk '{print $9}')
SA_O_v=$(cat OUTPUT_$ARG1 | grep "(m^2/cm^3)" | head -4 | tail -1 | awk '{print $9}')
SA_F_v=$(cat OUTPUT_$ARG1 | grep "(m^2/cm^3)" | head -5 | tail -1 | awk '{print $9}')
SA_S_v=$(cat OUTPUT_$ARG1 | grep "(m^2/cm^3)" | head -6 | tail -1 | awk '{print $9}')
SA_Br_v=$(cat OUTPUT_$ARG1 | grep "(m^2/cm^3)" | head -7 | tail -1 | awk '{print $9}')

SA_H_g=$(cat OUTPUT_$ARG1 | grep "(m^2/g)" | head -1 | tail -1 | awk '{print $9}')
SA_C_g=$(cat OUTPUT_$ARG1 | grep "(m^2/g)" | head -2 | tail -1 | awk '{print $9}')
SA_N_g=$(cat OUTPUT_$ARG1 | grep "(m^2/g)" | head -3 | tail -1 | awk '{print $9}')
SA_O_g=$(cat OUTPUT_$ARG1 | grep "(m^2/g)" | head -4 | tail -1 | awk '{print $9}')
SA_F_g=$(cat OUTPUT_$ARG1 | grep "(m^2/g)" | head -5 | tail -1 | awk '{print $9}')
SA_S_g=$(cat OUTPUT_$ARG1 | grep "(m^2/g)" | head -6 | tail -1 | awk '{print $9}')
SA_Br_g=$(cat OUTPUT_$ARG1 | grep "(m^2/g)" | head -7 | tail -1 | awk '{print $9}')

echo $ARG1"   "$SA_H_v"   "$SA_C_v"   "$SA_N_v"   "$SA_O_v"   "$SA_F_v"   "$SA_S_v"   "$SA_Br_v >> results_v.dat
echo $ARG1"   "$SA_H_g"   "$SA_C_g"   "$SA_N_g"   "$SA_O_g"   "$SA_F_g"   "$SA_S_g"   "$SA_Br_g >> results_g.dat

done

