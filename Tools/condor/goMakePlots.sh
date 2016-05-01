#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

cd $2/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}

echo "xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') ."
xrdcp root://cmseos.fnal.gov/$(echo $6 | sed 's|/eos/uscms||') .

ls

./makePlots -st --condor -D $1 -N $3 -M $4 -L $5 -S SB_37_2015

ls
xrdcp histoutput_$1_$4.root root://cmseos.fnal.gov//store/user/sbein/HadStop/HistOutput
rm histoutput_$1_$4.root
rm $(echo $6 | sed 's|.*/||')

