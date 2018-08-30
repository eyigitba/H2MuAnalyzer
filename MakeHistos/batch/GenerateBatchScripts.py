#! /usr/bin/env python

#######################################################
###              GenerateBatchScript.py             ###
###                                                 ### 
###  Generates a script that will submit jobs that  ###
###    make histograms to the lxplus batch queue.   ###
###      Run as ./batch/GenerateBatchScript.py      ###
###  The output can then be run as ./submit_all.sh  ###
###                                                 ###
###                Andrew Brinkerhoff               ###
###                     20.08.2018                  ###
#######################################################

## Basic python includes for manipulating files
import sys
import os

## Specific functions used below
from shutil import rmtree
from subprocess import PIPE,Popen
from operator import itemgetter

## Info about data and MC NTuples from H2MuAnalyzer/MakeHistos/python/SampleDatabase.py
sys.path.insert(0, '%s/python' % os.getcwd())
from SampleDatabase import GetSamples

## Command to list files in eos on lxplus
eos_cmd = '/afs/cern.ch/project/eos/installation/ams/bin/eos.select'

MACRO = 'macros/ReadNTupleChain.C'  ## Root macro to run from each job
LOC   = 'CERN'  ## Location of input files ('CERN' or 'UF')
YEAR  = 2017    ## Dataset year (2016 or 2017)

OUT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'  ## Directory for logs and output root files
LABEL   = 'Data_Aug18_v1'  ## Unique label for this set of jobs

NJOBS   =    20  ## Maximum number of jobs to generate
JOBSIZE =   300  ## Size of input NTuples in MB, per job (default 1000)

MAX_EVT = -1    ## Maximum number of events to process per job
PRT_EVT = 1000  ## Print every Nth event in each job

DATA_ONLY = True   ## Only process data samples, not MC
MC_ONLY   = False  ## Only process MC samples, not data
SIG_ONLY  = False  ## Only process signal MC samples, no others

VERBOSE = False ## Verbose printout


## Function to write the launcher script for a single job, and add that job to the main submit_all.sh script
def WriteSingleJob(subs_file, sub_files, samp_name, in_dir_name, file_list):

    out_dir = OUT_DIR+'/'+LABEL

    job_name      = 'sub_%d_%s' % (len(sub_files), samp_name)
    launcher_name = 'batch/launchers/%s.sh' % job_name

    ## In submit_all.sh (subs_file), write a line that will submit a job (bsub) to the queue of lxplus machines
    ## that run jobs for up to 1 hour (-q 1nh), specifying the log and error output file location (-o, -e) and
    ## the script that will be run by this job (${run_dir}/batch/launchers/%s.sh)
    subs_file.write( '\nbsub -q 1nh -o ${out_dir}/log/%s.log -e ${out_dir}/err/%s.err ${run_dir}/batch/launchers/%s.sh' % (job_name, job_name, job_name) )

    sub_files.append( open(launcher_name, 'w') )

    ## Write the line defining how to run MACRO, specifying the sample and input files to run over
    run_macro  = "\nroot -b -l -q '%s/%s(" % (os.getcwd(), MACRO)
    run_macro += ('"'+samp_name+'", "'+in_dir_name+'", "'+OUT_DIR+'/'+LABEL+'/files", {"')  ## sample, in_dir, out_dir
    for in_file in file_list:  ## in_files
        run_macro += (in_file+'", "')
    run_macro  = run_macro[:-4]  ## Remove last ", "
    run_macro += '"}, "%d", ' % (len(sub_files)-1)  ## out_file_str
    run_macro += "%d, %d)'" % (MAX_EVT, PRT_EVT)  ## max_evt, prt_evt

    sub_files[-1].write('\nrun_dir="%s"' % os.getcwd())
    sub_files[-1].write('\ncd ${run_dir}')
    sub_files[-1].write('\neval `scramv1 runtime -sh`')
    sub_files[-1].write(run_macro)
    print 'Wrote file %s' % launcher_name

## End function: WriteSingleJob
                      

## Main function executed by ./batch/GenerateBatchScript.py
def main():

    print '\n\n*** Inside GenerateBatchScripts.py ***\n'

    ## Create output directories for plots and submission scripts
    print 'Preparing output directories'
    out_dir = OUT_DIR+'/'+LABEL
    if os.path.exists(out_dir):
        delete_dir = raw_input('\n*** Directory %s already exists!!! ***\nType "Y" to delete and continue, "N" to exit.\n(Will also delete current batch/launchers/ folder.)\n\n' % out_dir)
        if delete_dir == 'Y':
            rmtree(out_dir)
            rmtree('batch/launchers')
            print '\nDeleted %s and batch/launchers\n' % out_dir
        else:
            print 'You typed %s, not "Y" - exiting\n' % delete_dir.lower()
    os.makedirs(out_dir)
    os.makedirs(out_dir+'/files')
    os.makedirs(out_dir+'/log')
    os.makedirs(out_dir+'/err')
    os.makedirs('batch/launchers')

    subs_file = open('submit_all.sh', 'w') ## Master submission script
    # hadd_file = open('hadd_all.sh', 'w') ## Script to hadd output files

    ## Set directories to run the jobs from, and to output files to
    subs_file.write('\npwd_cmd="/bin/pwd"')
    subs_file.write('\nrun_dir=${pwd_cmd}')
    subs_file.write('\nout_dir="%s"\n' % out_dir)

    # hadd_file.write('\nout_dir="%s"\n' % out_dir)
    
    sub_files = [] ## Separate submission script for each job

    ## Get list of samples from SampleDatabase.py
    samples = GetSamples(LOC, YEAR)

    ## Loop over available samples
    for samp in samples:

        if (YEAR == 2017):  ## Some samples not yet available for 2017 - AWB 17.08.2018
            if ('_120' in samp.name or '_130' in samp.name): continue
            if ('ZJets' in samp.name):
                if not (samp.name == 'ZJets_AMC' or samp.name == 'ZJets_AMC_2' or samp.name == 'ZJets_m_10_50'): continue
            if ('tt_' in samp.name or 'tW_' in samp.name or samp.name == 'tZW'): continue
            if (samp.name == 'WZ_2l' or samp.name == 'ZZ_2l_2q'): continue

        if DATA_ONLY and not samp.evt_type == 'Data':
            continue
        if MC_ONLY and samp.evt_type == 'Data':
            continue
        if SIG_ONLY and not samp.evt_type == 'Sig':
            continue
        
        print '\nLooking at sample %s' % samp.name
        
        ########################################################
        ###  Get list of files, and their size and location  ###
        ########################################################
        
        in_dir_name = samp.in_dir+'/'+samp.DAS_name+'/'+samp.name
        
        versions = []  ## In case of multiple crab submissions
        print 'Running command eos ls %s' % in_dir_name
        # for ver in subprocess.check_output([eos_cmd, 'ls', in_dir_name]).splitlines():  ## Only available in Python >= 2.7
        eos_ls = Popen([eos_cmd, 'ls', in_dir_name], stdout=PIPE)
        for ver in eos_ls.communicate()[0].split():
            if 'root' in ver: continue
            if VERBOSE: print '  * Appending [%d, %d]' % ( int(ver.split('_')[0]), int(ver.split('_')[1]) )
            versions.append([int(ver.split('_')[0]), int(ver.split('_')[1])])

        if len(versions) > 0:
            versions.sort(key = itemgetter(0, 1), reverse=True)  ## Choose the latest crab submission
        else:
            print '\n\nWARNING!!!  No crab output found for sample %s, from DAS %s' % (samp.name, samp.DAS_name)
            print 'Looked in %s - maybe it is somewhere else?\n\n' % in_dir_name
            continue

        print 'Chose version %d_%d' % (versions[0][0], versions[0][1])
        in_dir_name += '/%d_%d' % (versions[0][0], versions[0][1])
        
        in_files = [] ## List of input files with their size in MB
        eos_ls = Popen([eos_cmd, 'ls', in_dir_name], stdout=PIPE)
        for subdir in eos_ls.communicate()[0].split():
            eos_ls = Popen([eos_cmd, 'ls', in_dir_name+'/'+subdir], stdout=PIPE)
            print 'Running command eos ls %s' % (in_dir_name+'/'+subdir)
            for in_file in eos_ls.communicate()[0].split():
                if 'tuple' in in_file and '.root' in in_file: ## Only look at tuple_*.root files
                    eos_du = Popen([eos_cmd, 'ls', '-l', in_dir_name+'/'+subdir+'/'+in_file], stdout=PIPE)
                    fileMB = int(eos_du.communicate()[0].split()[4]) / 1000000. ## Get file size in MB
                    in_files.append(['%s/%s' % (subdir, in_file), fileMB])

        if len(in_files) > 0:
            print 'Final file list has %d entries' % len(in_files)
        else:
            print '\n\nWARNING!!!  No output root files found for sample %s, from DAS %s' % (samp.name, samp.DAS_name)
            print 'Looked in %s - maybe it is somewhere else?\n\n' % in_dir_name
            continue


        ##############################################################
        ###  Submit jobs for groups of files with the proper size  ###
        ##############################################################

        job_size  = 0.  ## Size of jobs in each input file in MB
        job_files = []  ## Files submitted to a single job 
        for iFile in range(len(in_files)):
            if (len(sub_files) >= NJOBS - 1 and NJOBS > 0):
                break
            if (job_size > JOBSIZE and JOBSIZE > 0):
                WriteSingleJob(subs_file, sub_files, samp.name, in_dir_name, job_files)
                if VERBOSE: print 'Writing job for sample %s, %d job files from %s to %s' % (samp.name, len(job_files), job_files[0], job_files[-1])
                job_size  = 0.
                job_files = []
            job_files.append( in_files[iFile][0] )
            job_size += in_files[iFile][1]
        ## End loop: for iFile in range(len(in_files))
        WriteSingleJob(subs_file, sub_files, samp.name, in_dir_name, job_files)
        if VERBOSE: print 'Writing job for sample %s, %d job files from %s to %s' % (samp.name, len(job_files), job_files[0], job_files[-1])

        print 'We have now written a total of %d launcher files' % len(sub_files)

        if (len(sub_files) >= NJOBS and NJOBS > 0):
            break

    ## End loop: for samp in samples

    ## Write the output files
    subs_file.close()
    print 'Wrote %s' % subs_file.name
    os.chmod(subs_file.name, 0o777) ## Render submit_all.sh executable
    for sub_file in sub_files:
        sub_file.close()
        print sub_file.name
        os.chmod(sub_file.name, 0o777)

        
## End function: main()
    

if __name__ == '__main__':
    main()
