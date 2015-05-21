JMEValidator
======

Common framework for CMS JEC / JER analyzes.

### Recipe


```sh
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_2_patch1
cd CMSSW_7_4_2_patch1/src/
cmsenv

git cms-addpkg CommonTools/PileupAlgos

# Puppi
git cms-merge-topic nhanvtran:puppi-etadep-741

# E/Gamma ID
git cms-merge-topic ikrav:egm_id_74X_v0

# Framework
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_74X
git clone git@github.com:blinkseb/TreeWrapper.git JMEAnalysis/TreeWrapper
git clone git@github.com:cms-jet/JMEValidator.git JMEAnalysis/JMEValidator -b CMSSW_7_4_X

scram b -j8

cd JMEAnalysis/JMEValidator/test
cmsRun runFramework.py
```
