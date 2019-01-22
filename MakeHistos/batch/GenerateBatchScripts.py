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
from SampleHelper import GetNormForSample, GetSampleID

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Root macro to run from each job
# MACRO = 'macros/ReadNTupleChain.C'
# MACRO = 'macros/MC_data_comparison.C'

MACRO = 'macros/WH_lep.C'
# MACRO = 'macros/ttH_3l.C'
#MACRO = 'macros/MiniNTupliser.C'
#MACRO = 'macros/lepMVA_efficiency.C'
#MACRO = 'macros/lepMVA_variables.C'

LOC    = 'CERN'  ## Location of input files ('CERN', 'CERN_hiM', or 'UF')
#LOC   = 'CERN_lepMVA_test_v2'  ## Location of input files ('CERN', 'CERN_hiM', or 'UF', or 'CERN_lepMVA_test_v1')
YEAR   = 2017    ## Dataset year (2016 or 2017)
LUMI   = 41000   ## 36814 for 2016, 41000 for 2017
## Override default sample location from SampleDatabase.py (use IN_DIR = '' to keep default)
IN_DIR  = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1'
HADD_IN = True   ## Use pre-hadded root files (NTuple_*.root) instead of original files (tuple_*.root)

## Directory for logs and output root files
if USER == 'abrinke1': OUT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms'
if USER == 'xzuo':     OUT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

LABEL  = 'WH_lep_AWB_2019_01_19_lepMVA_test_v1'
# LABEL  = 'ttH_3l_AWB_2019_01_19_lepMVA_test_v1'
#LABEL = 'lepMVA_variables_v3'
#LABEL = 'miniNtuple_WH_2017_v5'
#LABEL   = 'WH_cat_2017_v4_v4' 

NJOBS   =   -1  ## Maximum number of jobs to generate
JOBSIZE =  100  ## Size of input NTuples in MB, per job (default 1000)

MAX_EVT = -1     ## Maximum number of events to process per job
PRT_EVT = 10000  ## Print every Nth event in each job

DATA_ONLY = False  ## Only process data samples, not MC
MC_ONLY   = False  ## Only process MC samples, not data
SIG_ONLY  = False  ## Only process signal MC samples, no others
SAMP_LIST = []  ## Leave [] empty to process multiple samples

#SAMP_LIST = ['ZJets_AMC', 'tt',              # missing single top  --XWZ 28.09.2018
#            'tZq', 'ttW','ttZ','ttH'        # tx and ttX, so far only tZq for tx, missing 'tW', 'tZW' -XWZ 27.09.2018
#            'WW', 'WZ_3l_AMC', 'ZZ_2l_2v', 'ZZ_4l',  # diboson samples, missing 'WZ_2l' and 'ZZ_2l_2q'  --XWZ 27.09.2018
#            'WWW', 'WWZ', 'WZZ', 'ZZZ',     # triboson, all the samples at hand included   - XWZ 27.09.2018
#            'H2Mu_gg', 'H2Mu_VBF', 'H2Mu_ZH', 'H2Mu_WH_pos', 'H2Mu_WH_neg', 'H2Mu_ttH']  ## for keeping track of what is used

# SAMP_LIST = ['ZJets_AMC_1j_A']  ## Name of individual samples to process ([] to process multiple samples)
# SAMP_LIST = ['H2Mu_WH_pos_125']  ## Name of individual samples to process ([] to process multiple samples)

VERBOSE = False ## Verbose printout


## Function to write the launcher script for a single job, and add that job to the main submit_all.sh script
def WriteSingleJob(subs_file, runs_file, sub_files, samp_name, in_dir_name, file_list, samp_wgt):

    out_dir = OUT_DIR+'/'+LABEL

    job_name      = 'sub_%d_%s' % (len(sub_files), samp_name)
    launcher_name = 'batch/launchers/%s.sh' % job_name

    ## In submit_all.sh (subs_file), write a line that will submit a job (bsub) to the queue of lxplus machines
    ## that run jobs for up to 1 hour (-q 1nh), specifying the log and error output file location (-o, -e) and
    ## the script that will be run by this job (${run_dir}/batch/launchers/%s.sh)
    subs_file.write( '\nbsub -q 8nh -o ${out_dir}/log/%s.log -e ${out_dir}/err/%s.err ${run_dir}/batch/launchers/%s.sh' % (job_name, job_name, job_name) )
    runs_file.write( '\n${run_dir}/batch/launchers/%s.sh' % job_name )

    sub_files.append( open(launcher_name, 'w') )

    ## Write the line defining how to run MACRO, specifying the sample and input files to run over
    run_macro  = "\nroot -b -l -q '%s/%s(" % (os.getcwd(), MACRO)
    run_macro += ('"'+samp_name+'", "'+in_dir_name+'", "'+OUT_DIR+'/'+LABEL+'/files", {"')  ## sample, in_dir, out_dir
    for in_file in file_list:  ## in_files
        run_macro += (in_file+'", "')
    run_macro  = run_macro[:-4]  ## Remove last ", "
    run_macro += '"}, "%d", ' % (len(sub_files)-1)  ## out_file_str
    run_macro += "%d, %d, %f)'" % (MAX_EVT, PRT_EVT, samp_wgt)  ## max_evt, prt_evt, sample_wgt

    sub_files[-1].write('\nrun_dir="%s"' % os.getcwd())
    sub_files[-1].write('\ncd ${run_dir}')
    sub_files[-1].write('\neval `scramv1 runtime -sh`')
    sub_files[-1].write(run_macro)
    sub_files[-1].close()
    print 'Wrote file %s' % sub_files[-1].name
    os.chmod(sub_files[-1].name, 0o777)

## End function: WriteSingleJob
                      

## Main function executed by ./batch/GenerateBatchScript.py
def main():

    print '\n\n*** Inside GenerateBatchScripts.py ***\n'

    ## Create output directories for plots and submission scripts
    print 'Preparing output directories'
    if os.path.exists('batch/launchers'):
        delete_dir = raw_input('\n*** batch/launchers/ directory already exists!!! ***\nType "Y" to delete and continue, "N" to exit.\n\n')
        if delete_dir == 'Y':
            rmtree('batch/launchers')
        else:
            print 'You typed %s, not "Y" - exiting\n' % delete_dir.lower()
    out_dir = OUT_DIR+'/'+LABEL
    if os.path.exists(out_dir):
        delete_dir = raw_input('\n*** Directory %s already exists!!! ***\nType "Y" to delete and continue, "N" to exit.\n\n' % out_dir)
        if delete_dir == 'Y':
            rmtree(out_dir)
            print '\nDeleted %s\n' % out_dir
        else:
            print 'You typed %s, not "Y" - exiting\n' % delete_dir.lower()
    os.makedirs(out_dir)
    os.makedirs(out_dir+'/files')
    os.makedirs(out_dir+'/log')
    os.makedirs(out_dir+'/err')
    os.makedirs('batch/launchers')

    subs_file = open('submit_all.sh', 'w') ## Master submission script
    runs_file = open('run_all.sh', 'w')    ## Master script to run all jobs locally
    # hadd_file = open('hadd_all.sh', 'w') ## Script to hadd output files

    ## Set directories to run the jobs from, and to output files to
    subs_file.write('\npwd_cmd="/bin/pwd"')
    runs_file.write('\npwd_cmd="/bin/pwd"')
    subs_file.write('\nrun_dir=`${pwd_cmd}`')
    runs_file.write('\nrun_dir=`${pwd_cmd}`')
    subs_file.write('\nout_dir="%s"\n' % out_dir)
    runs_file.write('\nout_dir="%s"\n' % out_dir)

    # hadd_file.write('\nout_dir="%s"\n' % out_dir)
    
    sub_files = [] ## Separate submission script for each job

    ## Get list of samples from SampleDatabase.py
    samples = GetSamples(LOC, YEAR)

    miss_samp_usr = list(SAMP_LIST)  ## List of user-specified samples that are not found in SampleDatabase.py
    miss_samp_db  = []               ## List of samples in SampleDatabase.py that are not found in crab output

    ## Loop over available samples
    for samp in samples:

        # if (YEAR == 2017):  ## Some samples not yet available for 2017 - AWB 17.08.2018
        #     if ('_120' in samp.name or '_130' in samp.name): continue
        #     if ('ZJets' in samp.name):
        #         if not (samp.name == 'ZJets_AMC' or samp.name == 'ZJets_AMC_2' or samp.name == 'ZJets_m_10_50'): continue
        #     if ('tt_' in samp.name or 'tW_' in samp.name or samp.name == 'tZW'): continue
        #     if (samp.name == 'WZ_2l' or samp.name == 'ZZ_2l_2q'): continue

        if DATA_ONLY and not samp.evt_type == 'Data':
            continue
        if MC_ONLY and samp.evt_type == 'Data':
            continue
        if SIG_ONLY and not samp.evt_type == 'Sig':
            continue
        if len(SAMP_LIST) > 0 and not samp.name in SAMP_LIST:
            continue

        print '\nLooking at sample %s' % samp.name
        miss_samp_db.append(samp.name)
        
        ########################################################
        ###  Get list of files, and their size and location  ###
        ########################################################
        
        if len(IN_DIR) != 0: in_dir_name = IN_DIR+'/'+samp.DAS_name+'/'+samp.name
        else:                in_dir_name = samp.in_dir+'/'+samp.DAS_name+'/'+samp.name
        
        versions = []  ## In case of multiple crab submissions

        print 'Running command ls %s' % in_dir_name
        # for ver in subprocess.check_output(['ls', in_dir_name]).splitlines():  ## Only available in Python >= 2.7
        ls_files = Popen(['ls', in_dir_name], stdout=PIPE)

        if not HADD_IN:  ## Using original tuple_*.root files produced by crab, not hadd-ed versions
            for ver in ls_files.communicate()[0].split():
                if 'root' in ver: continue
                if VERBOSE: print '  * Appending [%d, %d]' % ( int(ver.split('_')[0]), int(ver.split('_')[1]) )
                versions.append([int(ver.split('_')[0]), int(ver.split('_')[1]), ver])

            if len(versions) > 0:
                print versions
                print "\n"
                versions.sort(key = itemgetter(0, 1), reverse=True)  ## Choose the latest crab submission
            else:
                print '\n\nWARNING!!!  No crab output found for sample %s, from DAS %s' % (samp.name, samp.DAS_name)
                print 'Looked in %s - maybe it is somewhere else?\n\n' % in_dir_name
                continue
        
            print 'Chose version %d_%06d' % (versions[0][0], versions[0][1])
#	if samp.name is 'SingleMu_2017F':
#		in_dir_name += '/180802_164117'
#	if samp.name is 'ZJets_AMC':        # temporary for 2017 WH
#		in_dir_name += '/180802_165055'
#	else:
            in_dir_name += '/%d_%06d' % (versions[0][0], versions[0][1])
 
        in_files = [] ## List of input files with their size in MB
        ls_files = Popen(['ls', in_dir_name], stdout=PIPE)
        for subdir in ls_files.communicate()[0].split():
            if VERBOSE: print '  * subdir = %s' % subdir

            if HADD_IN and 'NTuple' in subdir and '.root' in subdir:
                    du_file = Popen(['ls', '-l', in_dir_name+'/'+subdir], stdout=PIPE)
                    fileMB  = int(du_file.communicate()[0].split()[4]) / 1000000. ## Get file size in MB
                    in_files.append(['%s' % subdir, fileMB])

            elif not HADD_IN:
                ls_files = Popen(['ls', in_dir_name+'/'+subdir], stdout=PIPE)
                print 'Running command ls %s' % (in_dir_name+'/'+subdir)
                for in_file in ls_files.communicate()[0].split():
                    if 'tuple' in in_file and '.root' in in_file: ## Only look at tuple_*.root files
                        du_file = Popen(['ls', '-l', in_dir_name+'/'+subdir+'/'+in_file], stdout=PIPE)
                        fileMB  = int(du_file.communicate()[0].split()[4]) / 1000000. ## Get file size in MB
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
        ## Get XSec / nProcessed for all files used in the sample, not only for this job
	samp_wgt = GetNormForSample(samp.name, samp.xsec, LUMI, in_dir_name, in_files)
        for iFile in range(len(in_files)):
            if (len(sub_files) >= NJOBS - 1 and NJOBS > 0):
                break
            if (job_size > JOBSIZE and JOBSIZE > 0):
                WriteSingleJob(subs_file, runs_file, sub_files, samp.name, in_dir_name, job_files, samp_wgt)
                if VERBOSE: print 'Writing job for sample %s, %d job files from %s to %s' % (samp.name, len(job_files), job_files[0], job_files[-1])
                job_size  = 0.
                job_files = []
            job_files.append( in_files[iFile][0] )
            job_size += in_files[iFile][1]
        ## End loop: for iFile in range(len(in_files))
        WriteSingleJob(subs_file, runs_file, sub_files, samp.name, in_dir_name, job_files, samp_wgt)
        if VERBOSE: print 'Writing job for sample %s, %d job files from %s to %s' % (samp.name, len(job_files), job_files[0], job_files[-1])

        print 'We have now written a total of %d launcher files' % len(sub_files)

        if (samp.name in miss_samp_usr): miss_samp_usr.remove(samp.name)  ## Remove sample name from the list of missing samples
        if (samp.name in miss_samp_db):  miss_samp_db .remove(samp.name)  ## Remove sample name from the list of missing samples

        if (len(sub_files) >= NJOBS and NJOBS > 0):
            break

    ## End loop: for samp in samples

    if len(miss_samp_usr) > 0:
        print '\n\nRequested %d samples, missing %d!!!' % ( len(SAMP_LIST), len(miss_samp_usr) )
        print miss_samp_usr
    if len(miss_samp_db) > 0:
        print '\n\nNo input ntuples found for %d samples in SampleDatabase.py!!!' % len(miss_samp_db)
        print miss_samp_db

    subs_file.write('\n\necho "Jobs will be output to ${out_dir}"\n')
    runs_file.write('\n\necho "Jobs will be output to ${out_dir}"\n')

    ## Write the output files
    subs_file.close()
    runs_file.close()
    print '\n\nWrote %s and %s :\n' % (subs_file.name, runs_file.name)
    os.chmod(subs_file.name, 0o777) ## Render submit_all.sh executable
    os.chmod(runs_file.name, 0o777) ## Render runs_all.sh executable
    for sub_file in sub_files:
        print sub_file.name

        
## End function: main()
    

if __name__ == '__main__':
    main()
