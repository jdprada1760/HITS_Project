#!/bin/bash
#$ -N JP_work3
#$ -q bridge.q
#$ -l s_rt=48:00:00
#$ -pe mvapich 32
#$ -wd /home/pradajs/Heidelberg/HITS_Project/Code
#$ -e /home/pradajs/Heidelberg/HITS_Project/Code/Jobs/
#$ -o /home/pradajs/Heidelberg/HITS_Project/Code/Jobs/
#$ -M jesus.prada@h-its.org
#$ -m bes

# Redshift dependence of semiaxes at virial radius (each redshift)
source /home/pradajs/Heidelberg/HITS_Project/Code/Jobs/bashrc
#python get_Axes.py level3_MHD 63
python get_Axes.py level3_MHD 21

