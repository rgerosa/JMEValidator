Setup

Install MVA MET using the following recipe

    cmsrel CMSSW_7_4_16_patch2
    cd CMSSW_7_4_16_patch2/src
    cmsenv
    
    git clone git@github.com:blinkseb/TreeWrapper.git JMEAnalysis/TreeWrapper
    git clone git@github.com:rgerosa/JMEValidator.git JMEAnalysis/JMEValidator -b classicMVAMET
    
    git cms-addpkg PhysicsTools/PatUtils
    sed '/corrMET, srcMET/a \ \ \ \ outMET\.setSignificanceMatrix\(srcMET\.getSignificanceMatrix\(\)\)\;' PhysicsTools/PatUtils/plugins/CorrectedPATMETProducer.cc -i
    git cms-addpkg DataFormats/METReco
    sed '/setSignificanceMatrix/a \ \ \ \ void setSumEt\(const double \& sumEt\)\{ sumet \= sumEt\; \}\;' DataFormats/METReco/interface/MET.h -i

An exemplary analysis can be run with
    cmsRun JMEAnalysis/JMEValidator/test/runFrameworkMC.py maxEvents=1000
