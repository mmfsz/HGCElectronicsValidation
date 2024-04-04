#!/bin/bash

cmssw=${1}
input=${2}
maxEvents=${3}
outdir=${4}
fold=${5}
adcThrMIP=${6}
adcThrMIPbxm1=${7}
geom=${8}
scaleByDoseFactor=${9}
jobTag=${10}

# make outdir in /afs 
mkdir -p $outdir

pwd

name=`basename ${input}`
thrbxm1=${adcThrMIPbxm1/./p}
thrbxm1=${thrbxm1/-/m}
outfilename=${name}_geo${geom}_thr${adcThrMIP/./p}_thrbxm1${thrbxm1}.root

pwd
cd $cmssw/src
pwd
eval `scram r -sh`
cd -
pwd

cmsRun $cmssw/src/UserCode/HGCElectronicsValidation/python/hgcScintAnalyzer_cfg.py \
    input=${input}/ \
    output=output.root \
    fold=${fold} \
    adcThrMIP=${adcThrMIP} adcThrMIPbxm1=${adcThrMIPbxm1} \
    maxEvents=${maxEvents} \
    geometry=${geom} \
    scaleByDoseFactor=${scaleByDoseFactor}


# move to output dir in eos
outdirname=`basename $outdir`
eosdir=/eos/user/m/mmazza/samples/HGCAL/HGCScintAnalyzer/${outdirname}/
mkdir -p ${eosdir}
mv -v output.root ${eosdir}/${outfilename}

#python $cmssw/src/UserCode/HGCElectronicsValidation/test/scripts/prepareOccupancySummary.py ${outdir}/${finalname}.root



