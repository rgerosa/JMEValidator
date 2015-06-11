cmssw
mkdir $1
cd $1
SCRAM_ARCH=slc6_amd64_gcc481
scram project CMSSW_7_1_0_pre4
cd CMSSW_7_1_0_pre4/src
cmsenv
git cms-addpkg CondFormats/EgammaObjects
git clone https://github.com/rfriese/MVAMet
#cp /afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_7_1_0_pre4/src/* . -r
scram b -j 4
