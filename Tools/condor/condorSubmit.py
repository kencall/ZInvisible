#!/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_8/external/slc6_amd64_gcc491/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

import sys
from os import system, environ
sys.path = [environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/condor/",] + sys.path

from samples import SampleCollection
import optparse 
import subprocess

mvaFileName = ""
with file(environ["CMSSW_BASE"] + "/src/ZInvisible/Tools/TopTagger.cfg") as meowttcfgFile:
    for line in meowttcfgFile:
        if "modelFile" in line:
            mvaFileName = line.split("=")[1].strip().strip("\"")
            break

#go make plots!
submitFileGMP = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakePlots.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/makeDataMCplots, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakePlots.sh, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/bTagEffHists.root, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/TTbarNoHad_NJetsISR.root, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/lepEffHists.root,  $ENV(CMSSW_BASE)/src/ZInvisible/Tools/njetWgtHists.root, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/dataMCweights.root, $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/CSVv2_ichep.csv, $ENV(CMSSW_BASE)/lib/$ENV(SCRAM_ARCH)/librecipeAUXOxbridgeMT2.so, $ENV(CMSSW_BASE)/src/opencv/lib/libopencv_core.so.3.1, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/%(trainingFile)s, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/TopTagger.cfg, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/ISRWeights.root
Output = logs/makePlots_$(Process).stdout
Error = logs/makePlots_$(Process).stderr
Log = logs/makePlots_$(Process).log
notify_user = ${LOGNAME}@FNAL.GOV
x509userproxy = $ENV(X509_USER_PROXY)

"""%{"trainingFile":mvaFileName} 

#go make top plots!
submitFileGTP = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeTopPlots.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/makeTopPlots, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakePlots.sh, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/bTagEffHists.root, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/TTbarNoHad_NJetsISR.root,  $ENV(CMSSW_BASE)/src/ZInvisible/Tools/lepEffHists.root,  $ENV(CMSSW_BASE)/src/ZInvisible/Tools/njetWgtHists.root, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/dataMCweights.root, $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/CSVv2_ichep.csv, $ENV(CMSSW_BASE)/lib/$ENV(SCRAM_ARCH)/librecipeAUXOxbridgeMT2.so, $ENV(CMSSW_BASE)/lib/$ENV(SCRAM_ARCH)/libTopTaggerTopTagger.so, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/TopTagger.cfg, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/TopTagger_noMVA.cfg, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/TopTagger_AllComb.cfg, $ENV(CMSSW_BASE)/src/opencv/lib/libopencv_core.so.3.1, $ENV(CMSSW_BASE)/src/opencv/lib/libopencv_ml.so.3.1, #$ENV(CMSSW_BASE)/src/ZInvisible/Tools/%(trainingFile)s
Output = logs/makePlots_$(Process).stdout
Error = logs/makePlots_$(Process).stderr
Log = logs/makePlots_$(Process).log
notify_user = ${LOGNAME}@FNAL.GOV
x509userproxy = $ENV(X509_USER_PROXY)


"""%{"trainingFile":mvaFileName} 

#go make lepton efficiency
submitFileGME = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeEff.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/calcEff, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeEff.sh, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/zRes.root, $ENV(CMSSW_BASE)/lib/$ENV(SCRAM_ARCH)/librecipeAUXOxbridgeMT2.so
notify_user = ${LOGNAME}@FNAL.GOV
x509userproxy = $ENV(X509_USER_PROXY)

"""

#go B Efficiency calc
submitFileGBE = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeBeff.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/beffCalc, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeBeff.sh
Output = logs/makePlots_$(Process).stdout
Error = logs/makePlots_$(Process).stderr
Log = logs/makePlots_$(Process).log
notify_user = ${LOGNAME}@FNAL.GOV
x509userproxy = $ENV(X509_USER_PROXY)

"""

#go Signal Efficiency
submitFileGSE = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeSigEff.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ZInvisible/Tools/makeSignalHistograms, $ENV(CMSSW_BASE)/src/ZInvisible/Tools/condor/goMakeSigEff.sh
Output = logs/makePlots_$(Process).stdout
Error = logs/makePlots_$(Process).stderr
Log = logs/makePlots_$(Process).log
notify_user = ${LOGNAME}@FNAL.GOV
x509userproxy = $ENV(X509_USER_PROXY)

"""

parser = optparse.OptionParser("usage: %prog [options]\n")

parser.add_option ('-n',  dest='numfile', type='int', default = 5, help="number of files per job")
parser.add_option ('-d',  dest='datasets', type='string', default = '', help="List of datasets 'ZJetsToNuNu,DYJetsToLL'")
parser.add_option ('-l',  dest='dataCollections', action='store_true', default = False, help="List all datacollections")
parser.add_option ('-L', dest='dataCollectionslong', action='store_true', default = False, help="List all datacollections and sub collections")
parser.add_option ('-r',  dest='refLumi', type='string', default = None, help="Data collection to define lumi (uses default lumi if no reference data collection is defined)")
parser.add_option ('-c',  dest='noSubmit', action='store_true', default = False, help="Do not submit jobs.  Only create condor_submit.txt.")
parser.add_option ('-e',  dest='goMakeEff', action='store_true', default = False, help="Run calcEff instead of makePlots.")
parser.add_option ('-b',  dest='goMakeBeff', action='store_true', default = False, help="Run beffCalc instead of makePlots.")
parser.add_option ('-s',  dest='goMakeSigEff', action='store_true', default = False, help="Run makeSignalHistograms instead of makePlots.")
parser.add_option ('-t',  dest='goMakeTopPlots', action='store_true', default = False, help="Run makeTopPlots instead of makePlots.")

options, args = parser.parse_args()

submitFile = ""
exeName = ""

if options.goMakeEff:
    exeName = "calcEff"
    submitFile = submitFileGME
elif options.goMakeBeff:
    exeName = "beffCalc"
    submitFile = submitFileGBE
elif options.goMakeSigEff:
    exeName = "makeSignalHistograms"
    submitFile = submitFileGSE
elif options.goMakeTopPlots:
    exeName = "makeTopPlots"
    submitFile = submitFileGTP
else:
    exeName = "makeDataMCplots"#"makePlots"
    submitFile = submitFileGMP

nFilesPerJob = options.numfile

fileParts = [submitFile]
sc = SampleCollection()
datasets = []

if options.dataCollections or options.dataCollectionslong:
    scl = sc.sampleCollectionList()
    for sampleCollection in scl:
        sl = sc.sampleList(sampleCollection)
        print sampleCollection
        if options.dataCollectionslong:
            sys.stdout.write("\t")
            for sample in sl:
                sys.stdout.write("%s  "%sample[1])
            print ""
            print ""
    exit(0)

if options.datasets:
    datasets = options.datasets.split(',')
else:
    print "No dataset specified"
    exit(0)

lumis = sc.sampleCollectionLumiList()
lumi = sc.getFixedLumi()
if options.refLumi != None:
    lumi = lumis[options.refLumi]
    print "Normalizing to %s pb-1" % (lumi)

for ds in datasets:
    ds = ds.strip()

    print ds
    for s, n in sc.sampleList(ds):
        print "\t%s"%n
        f = open(s)
        if not f == None:
            count = 0
            for l in f:
                if '.root' in l and not 'failed' in l:
                    count = count + 1
            for startFileNum in xrange(0, count, nFilesPerJob):
                fileParts.append("Arguments = %s $ENV(CMSSW_BASE) %i %i %f %s\n"%(n, nFilesPerJob, startFileNum, lumi, s))
                fileParts.append("Output = logs/%s_%s_%i.stdout\n"%(exeName, n, startFileNum))
                fileParts.append("Error = logs/%s_%s_%i.stderr\n"%(exeName, n, startFileNum))
                fileParts.append("Log = logs/%s_%s_%i.log\n"%(exeName, n, startFileNum))
                fileParts.append("Queue\n\n")

            f.close()

fout = open("condor_submit.txt", "w")
fout.write(''.join(fileParts))
fout.close()

if not options.noSubmit: 
    system('mkdir -p logs')
    system("echo 'condor_submit condor_submit.txt'")
    system('condor_submit condor_submit.txt')

