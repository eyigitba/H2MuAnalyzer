#! /usr/bin/env python

############################################################
####                    BatchGenerator.py               ####
####     Initialized for lepMVA working point study,    ####
####  generate python macros and HTcondor submission    ####
####  scripts for given set of lepMVA cut values.       ####
####     output can be run as ./submit_all.sh           ####
############################################################


## Basic python includes for manipulating files
import sys
import os
import array

## Specific functions used below
from shutil import rmtree
from subprocess import PIPE,Popen
from operator import itemgetter


sys.path.insert(0, '%s/python' % os.getcwd())

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'
## Directory for logs and output root files
if USER == 'abrinke1': OUT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms'
if USER == 'xzuo':     OUT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

LABEL = 'WH_ele_loose_ID_loose_iso_loose_mu_iso_v1'

PLOT_DIR = 'plots/lepMVA_scan'

LIB_NAME = 'working_point_3l'
MACRO = 'MakeMassStack'

mu1_MVAcuts = [ round(1.0*MVA/5 - 1.0 , 2) for MVA in range(10) ]
mu2_MVAcuts = [ round(1.0*MVA/5 - 1.0 , 2) for MVA in range(10) ]
lep_MVAcuts = [ -0.4, 0.4, 0.8]
POINTS_PER_JOB = 1


def WriteSingleScript( script_number, mu1_cuts, mu2_cuts, lep_cuts):
    in_file_name  = OUT_DIR + '/' + LABEL + '/' + 'all_samples.root'
    out_file_dir  = OUT_DIR + '/' + LABEL + '/' + PLOT_DIR + '/' + 'files' + '/'
    out_file_name = 'out_script_%d.root'%script_number
    script = open('batch/scripts/script_%d.py' %script_number, 'w')

    script.write('#! /usr/bin/env python\n\n')
    script.write('#############################################################\n')
    script.write('###                     python script %d                ###  \n'%script_number)    
    script.write('###   can be run directly as ./batch/scripts/script_%d  ###  \n'%script_number)
    script.write('###      or be submitted to batch via HTCondor            ###\n')
    script.write('#############################################################\n')
    script.write('\n\n')
    script.write('import sys\n')
    script.write('import os\n')
    script.write('sys.path.insert(0, "%s/python" % os.getcwd())\n')
    script.write('\n')
    script.write('from ROOT import *\n')
    script.write('from %s import %s\n'%(LIB_NAME, MACRO))    
    script.write('gROOT.SetBatch(True)\n')
    script.write('\n')
    script.write('mu1_cuts = %s\n'%mu1_cuts)
    script.write('mu2_cuts = %s\n'%mu2_cuts)
    script.write('lep_cuts = %s\n'%lep_cuts)
    script.write('def main():\n')
    script.write('    in_file_name  = "%s"\n' %in_file_name)
    script.write('    out_file_dir  = "%s"\n' %out_file_dir)
    script.write('    out_file_name = "%s"\n' %out_file_name)   
    script.write('\n')
    script.write('    for mu1_cut in mu1_cuts:\n')
    script.write('        for mu2_cut in mu2_cuts:\n')
    script.write('            for lep_cut in lep_cuts:\n')
    script.write('                MakeMassStack(in_file_name, out_file_dir, out_file_name, "WH", "lep_3l", mu1_cut, mu2_cut, lep_cut)\n')
    script.write('\n')
    script.write('\n')
    script.write('main()\n')

    print "Wrote script_%d\n" %script_number
    os.chmod(script.name, 0o777)


def WriteSingleSubmission( sub_all, run_all, sub_number):
    out_dir = OUT_DIR + '/' + LABEL + '/' + PLOT_DIR 
    submission = open('batch/condors/sub_%d.sub'%sub_number, 'w')
    
    submission.write('run_dir = /afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/ForMiniNTuple\n')
    submission.write('out_dir = %s\n' %out_dir)
    submission.write('executable      = $(run_dir)/batch/scripts/script_%d.py\n'%sub_number)
    submission.write('arguments       = $(ClusterID) $(ProcId)\n')
    submission.write('output          = $(out_dir)/out/sub_%d.$(ClusterId).$(ProcId).out\n'%sub_number)
    submission.write('error           = $(out_dir)/err/sub_%d.$(ClusterId).$(ProcId).err\n'%sub_number)
    submission.write('log             = $(out_dir)/log/sub_%d.$(ClusterId).log\n'%sub_number)
    submission.write('queue\n')

    submission.close()
    sub_all.write('condor_submit ${run_dir}/batch/condors/sub_%d.sub\n' %sub_number)
    run_all.write('${run_dir}/batch/scripts/script_%d.py\n' %sub_number)
    print "Wrote sub_%d.sub\n" %sub_number

 

def main():

    print "\n****running Batch Generator.py****\n"
    ## Create output directories for plots and submission scripts
    print 'Preparing output directories'
    if os.path.exists('batch/scripts'):
	delete_dir = raw_input('\n*** batch/scripts/ directory already exists!!! ***\nType "Y" to delete and continue, "N" to exit.\n\n')
	if delete_dir == 'Y':
            rmtree('batch/scripts')
        else:
            print 'You typed %s, not "Y" - exiting\n' % delete_dir.lower()
    if os.path.exists('batch/condors'):
        delete_dir = raw_input('\n*** batch/condors/ directory already exists!!! ***\nType "Y" to delete and continue, "N" to exit.\n\n')
        if delete_dir == 'Y':
            rmtree('batch/condors')
        else:
            print 'You typed %s, not "Y" - exiting\n' % delete_dir.lower()
    out_dir = OUT_DIR+'/'+LABEL + '/' + PLOT_DIR
    if os.path.exists(out_dir):
        delete_dir = raw_input('\n*** Directory %s already exists!!! ***\nType "Y" to delete and continue, "N" to exit.\n\n' % out_dir)
        if delete_dir == 'Y':
            rmtree(out_dir)
            print '\nDeleted %s\n' % out_dir
        else:
            print 'You typed %s, not "Y" - exiting\n' % delete_dir.lower()
    os.makedirs(out_dir)
    os.makedirs(out_dir+'/files')
    os.makedirs(out_dir+'/out')
    os.makedirs(out_dir+'/log')
    os.makedirs(out_dir+'/err')
    os.makedirs('batch/condors')
    os.makedirs('batch/scripts')

    subs_file = open('submit_all.sh', 'w') ## Master submission script
    runs_file = open('run_all.sh', 'w')    ## Master script to run all jobs locally

    ## Set directories to run the jobs from, and to output files to
    subs_file.write('\npwd_cmd="/bin/pwd"')
    runs_file.write('\npwd_cmd="/bin/pwd"')
    subs_file.write('\nrun_dir=`${pwd_cmd}`')
    runs_file.write('\nrun_dir=`${pwd_cmd}`')
    subs_file.write('\nout_dir="%s"\n' % out_dir)
    runs_file.write('\nout_dir="%s"\n' % out_dir)

    job_number = 0
    for mu1_cut in mu1_MVAcuts:
	for mu2_cut in mu2_MVAcuts:
	    for lep_cut in lep_MVAcuts:    
		mu1_cuts = [mu1_cut]
		mu2_cuts = [mu2_cut]
		lep_cuts = [lep_cut]
		WriteSingleScript( job_number, mu1_cuts, mu2_cuts, lep_cuts)
		WriteSingleSubmission( subs_file, runs_file, job_number)
		job_number += 1

    subs_file.write('\n\necho "Jobs will be output to ${out_dir}"\n')
    runs_file.write('\n\necho "Jobs will be output to ${out_dir}"\n')

    ## Write the output files
    subs_file.close()
    runs_file.close()
    print '\n\nWrote %s and %s \n' % (subs_file.name, runs_file.name)
    os.chmod(subs_file.name, 0o777) ## Render submit_all.sh executable
    os.chmod(runs_file.name, 0o777) ## Render runs_all.sh executable


main()



