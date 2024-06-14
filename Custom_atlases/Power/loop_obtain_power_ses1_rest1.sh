#!/bin/sh

#make sure to activate your enviornment with pandas and nipype prior to running this code. See script references here for how to setup enviornment. 
#conda activate holsmon

module load anaconda
module load fsl
conda activate holsmon

for SUB in 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041 042 043 044 045 046 047 048 049 050 051 052 053 054 055 056 057 058 059 060 061 062 063 064 065 070 072 073 074 075 076 077 078 079 080 081 082 083 084 085 086 087 088 089 090 091 092; do

		python -m obtain_power_ses1_rest1.py $SUB
done
