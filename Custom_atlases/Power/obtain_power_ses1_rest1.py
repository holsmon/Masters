#!/usr/local/bin/python3.10

#NOTES: Must have an enviornment setup that has pandas and nipype install. Nipype is on python version 3.10

#Create enviornment. Rename myenv with whatever you want to call your enviornment. 
#module load anaconda 
#conda create --name myenv

#Install python librararies. Only need to do this once after creating initial enviornment above.
#conda install pandas
#conda install nipype

import os, sys, tempfile, glob, shutil, nipype, pandas
import subprocess

import nipype.interfaces.fsl as fsl
import numpy as np

#subj is defined in your loop script for the system arguments. This is your subject IDs (sub-001,002, etc.)
subj = "sub-%s" % str(sys.argv[1]).zfill(3)
#may need to change sess and task for your dataset
sess = "ses-001"
task = "task-rest1mb"

#Define your directories and define your ts.1d matrix name (stores power node matrix)
basedir = "/storage/group/nad12/default/fgh3/Users/Hollie/xcp_d/"
xcpdir = os.path.join(basedir,'xcp_d_output_no_fd_despike','output','xcp_d')
output_1d = os.path.join(basedir,'custom_atlases','work')
ts_matrix_name = "%s_%s_%s_power274_ts.1D" % (subj,sess,task)
tsname = os.path.join(output_1d, ts_matrix_name)
outputdir = os.path.join(basedir,'custom_atlases','subjects','power',sess,subj)

#This is your template directory, which contains the power atlas in MNI space you want the nodes in
templatedir = os.path.join(basedir,'custom_atlases','templates','mni152_versions')

#GET the timeseries from xcp_d preprocessed output using power ROIS (in MNI152 space) rest1. Uses preprocessed BOLD output from Xcp_d
meants = fsl.utils.ImageMeants()
meants.inputs.in_file = os.path.join(xcpdir, subj, sess, 'func', '%s_%s_%s_space-MNI152NLin6Asym_desc-denoised_bold.nii.gz' % (subj,sess,task))

#power atlas (in MNI152 space). What you want to use to extract information from preprocessed BOLD output from XCP_D 
meants.inputs.args = "--label=%s" % os.path.join(templatedir, 'power_MNI152.nii.gz')
meants.inputs.out_file = tsname
print(meants.cmdline)
res = meants.run()

#get cross corr coef of time series (ts)
timeseries = np.loadtxt(tsname,unpack=True)
ccval = np.corrcoef(timeseries)
ccval = np.nan_to_num(ccval)
np.fill_diagonal(ccval,0)
   
#save functional connectivity matrix as csv rather than 1D file   
np.savetxt(os.path.join(outputdir,'%s_%s_%s_power274_r.csv' % (subj,sess,task)),ccval,fmt='%f',delimiter=',')


