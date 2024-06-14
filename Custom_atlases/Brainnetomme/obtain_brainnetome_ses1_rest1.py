#!/usr/local/bin/python3.10

#NOTES: Must have an enviornment setup that has pandas and nipype install. Nipype is on python version 3.10

#Create enviornment. Rename myenv with whatever you want to call enviornment. 
#module load anaconda 
#conda create --name myenv

#Install python librararies. Only need to do this once after creating initial enviornment above.
#conda install pandas
#conda install nipype

import os, sys, tempfile, glob, shutil, nipype, pandas
import subprocess

import nipype.interfaces.fsl as fsl
import numpy as np

#These are defined in your look script for the system arguments. 
subj = "sub-%s" % str(sys.argv[1]).zfill(3)
sess = "ses-001"
task = "task-rest1mb"

#Define your directories
basedir = "/storage/group/nad12/default/fgh3/Users/Hollie/xcp_d/"
xcpdir = os.path.join(basedir,'xcp_d_output_no_fd_despike','output','xcp_d')
output_1d = os.path.join(basedir,'custom_atlases','work')
ts_matrix_name = "%s_%s_%s_bna274_ts.1D" % (subj,sess,task)
tsname = os.path.join(output_1d, ts_matrix_name)
outputdir = os.path.join(basedir,'custom_atlases','subjects','bna',sess,subj)

#This is your template directory
templatedir = os.path.join(basedir,'custom_atlases','templates','mni152_versions')

#GET the timeseries from xcp_d preprocessed output and bna ROIS (in MNI152 space) rest1
meants = fsl.utils.ImageMeants()
meants.inputs.in_file = os.path.join(xcpdir, subj, sess, 'func', '%s_%s_%s_space-MNI152NLin6Asym_desc-denoised_bold.nii.gz' % (subj,sess,task))

#bna atlas (in MNI152 space). uses 274 rois.
meants.inputs.args = "--label=%s" % os.path.join(templatedir, 'bna_MNI152.nii.gz')
meants.inputs.out_file = tsname
print(meants.cmdline)
res = meants.run()

#get cross corr coef of TS
timeseries = np.loadtxt(tsname,unpack=True)
ccval = np.corrcoef(timeseries)
ccval = np.nan_to_num(ccval)
np.fill_diagonal(ccval,0)
   
#save functional connectivity matrix as csv   
np.savetxt(os.path.join(outputdir,'%s_%s_%s_bna274_r.csv' % (subj,sess,task)),ccval,fmt='%f',delimiter=',')


