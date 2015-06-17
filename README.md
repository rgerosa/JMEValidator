JMEValidator 
======

###########################################
### Recipe For Jet and Isolation Studies ##
###########################################

```sh
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_2_patch1
cd CMSSW_7_4_2_patch1/src/
cmsenv

git cms-addpkg CommonTools/PileupAlgos

# Puppi
git cms-merge-topic nhanvtran:puppi-etadep-741

# Puppi inverted
git cms-merge-topic violatingcp:pupuppi-741

# E/Gamma ID
git cms-merge-topic ikrav:egm_id_74X_v0

# Framework
git clone git@github.com:blinkseb/JetToolbox.git JMEAnalysis/JetToolbox -b jet_flavor
git clone git@github.com:blinkseb/TreeWrapper.git JMEAnalysis/TreeWrapper
git clone git@github.com:blinkseb/JMEValidator.git JMEAnalysis/JMEValidator -b refactoring_74x_next

scram b -j8

cd JMEAnalysis/JMEValidator/test
cmsRun runFramework.py
```

###################################
### Recipe For PUPPI MET Studies ##
###################################

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
git clone git@github.com:blinkseb/JetToolbox.git JMEAnalysis/JetToolbox -b jet_flavor
git clone git@github.com:blinkseb/TreeWrapper.git JMEAnalysis/TreeWrapper
git clone git@github.com:rgerosa/JMEValidator.git JMEAnalysis/JMEValidator -b PUPPETMVAMet

scram b -j8

cd JMEAnalysis/JMEValidator/test
cmsRun runFramework.py
```

Few information about the PUPPI MET version of JME Validator package:

```sh

////// To Produce Trees ////////// :

1) the cmsRun configuration code is placed in test/runFramework.py. Some parmaters can be paresed by command line when cmsRun is called:

   globalTag : specify the global tag to be used

   isMC      : specify if are data = false or MC = true

   runPuppiMuonIso : if the sequence for calculating puppi isolation for muons has to be run or not

   muonIsoCone : isolation cone for muons

   dropAnalyzerDumpEDM : produce an edm outputmodule instead of a "plain" root file via TFileService. This is used for sake of tests, to check if collections added to the event are correctly put and filled

   runMVAPUPPETAnalysis : if false, the standard JMEValidator sequence will be run (used for jet studies), if true the special sequence for Puppi met and output is produced.

   muonTypeID : In case runMVAPUPPETAnalysis is true, define which type of muon id to apply. Possibilites are : Tight (standard tight id with delta beta PU correction for isolation, TightDBeta (is tight ID + PFweighted delta beta correction), TightPuppiNoMu (is tight ID + puppi isolation calculated without muons)

   electronTypeID : Tight or Medium, up to now taken from physics 14 cutBased working points .. to be update.

   tauID : tight, medium or loose according to tau POG recipe

   standard parameter that can be passed from command line : inputFiles=<> , outputFile= , maxEvents=...

2) most of the modules are organized in python/FrameworkConfiguration.py

   These is the list of modules that are called:
   
   a) jet tool box is called for making AK4 jets from packedCandidate of MINIAOD, building the following collections: PF, PFCHS and Puppi. Proper JEC are used, for Puppi only L2 and L3 corrections are applied. Then for CHS and PF also PU JetID and QGLikelihood are evaluated. The output are patJet collections.

   b) from MuonIsolationTools.py, the muon isolation sequence is re-run for producing PFWighted DeltaB corrected value maps

   c) if runPuppiMuonIso is true : the puppi producer is re-run using a specif tune for Isolation studies, i.e. larger cone, with or without muons. The isolation maps are produced as well.

   e) if runMVAPUPPETAnalysis is True : as a function of the muonTypeID label, the muon ID is applied producing a new collection of muons which pass the ID. Muon ID is defined in python/LeptonSelectionTools_cff.py, where the patMuonIDIsoSelector producer, defined in plugins, is used. A set of defined value maps for the isolation can be used. The type of isolation considered, is ruled by typeIsoVal flag: if 0 means no PU correction, 1 means standard dBeta, 2 means -rho*ConeArea, 3 means PFWeighted dBeta (value maps have to be given), 4 means puppi (maps have to be given).

   
   f) Electron ID value maps are produced, then if runMVAPUPPETAnalysis is True, a skim of the electron collection is done according to the defined electronIDType. A new collection is produced by patElectronIDIsoSelector

   g) Also tau ID is applied, then identfied tau leptons are cleaned from muons and electrons by PATTauCleaner

   h) Jets all the jets produced by jet tool box are cleaned from identified muons, electrons and cleaned taus using PATJetCleaner.

   i) Different MET are computed : pfCHSMet + TypeI corrections, PuppiMEt + TypeI correction using ak4PFCHSL2L3Corrector. Corrected and uncorrected values stored in slimmedMet

   l) if runMVAPUPPETAnalysis is true: track met, met from PU charged particles, puppi neutral met and puppi inverted met are calculated. Inverted means that a puppi sequence inverting the probability and the neutral cut is run to obtain a neutral PU particle enriched sample.

   m) always if runMVAPUPPETAnalysis is true: Zmumu, Zee and Ztatau candidates are formed from identified leptons, asking Mll [70,110] and opposite charge. Then all the Zcandidates are merged in a single collection and exactly one candidate is required.

   n) In parallel, all the identified lepton collections are merged and exactly two leptons are required .. so that we veto the presence of other leptons in the event.

   o) finally, the mva met producer is run to produce: recoil vectors for each met type, Zboson kinematic and apply a training and prodcue the mva met

   p) Zboson kinematic, two leading jets (puppi in this case), Nvertex and recoils are stored in a plain tree by means of interface/PUPPETAnalyzer.h


3) Specific classes to run the training and read information from the output training file are committed in interface/ and src/ directories, like:
   GBRTrainer.h, GBRTrainer2D.h ... etc

4) Specific codes to execute in order to run the regression or apply it on the fly directly on root plain trees are committed in the bin/ directory.
   They can be used directly in the CMSSW area since they are compiled by scram (BuildFile is placed in the bin dir as well)

```

How to run :

```sh

cmsRun runFramework.py runMVAPUPPETAnalysis=True

It produces an output root file with some directories .. in each one a tree with the information. 

In our case PUPPET directory. 

The output file contains few information with respect to the configuration with runMVAPUPPETAnalysis=False


```

