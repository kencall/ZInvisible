#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

cd $2/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$2/src/opencv/lib/

#cd /uscms_data/d3/snorberg/CMSSW_8_0_12/src/ZInvisible/Tools/


echo "xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') ."
xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .

ls

#./makePlots -st --condor -D $1 -N $3 -M $4 -L $5 -S SB_59_2016
./makePlots -s --condor -D $1 -N $3 -M $4 -L $5 -S SB_59_2016 #-E 1000

ls

rm $(echo $6 | sed 's|.*/||')

