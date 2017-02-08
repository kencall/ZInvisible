#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

cd $2/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}
#cd /uscms_data/d3/snorberg/CMSSW_8_0_12/src/ZInvisible/Tools/

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$2/src/opencv/lib/

#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$2/src/opencv/lib/

echo "xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') ."
xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .

ls

./makeDataMCplots -st --condor -D $1 -N $3 -M $4 -L $5 -S SB_v1_2017

ls

rm $(echo $6 | sed 's|.*/||')

