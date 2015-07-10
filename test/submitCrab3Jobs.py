from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.requestName = ''
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_StartupFlat10to50bx50Raw_MCRUN2_74_V8'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_AsymptNoPURawReco_MCRUN2_74_V9A'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_AsymptFlat10to50bx25Raw_MCRUN2_74_V9'
config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9'

config.section_('JobType')
config.JobType.psetName    = 'runFrameworkMC.py'
config.JobType.pluginName  = 'Analysis'
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A', 'isMC=True', 'runMVAPUPPETAnalysis=True', 'applyJECtoPuppiJets=False', 'applyZSelections=True', 'dropAnalyzerDumpEDM=False','applyWSelections=False']
config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9', 'isMC=True', 'runMVAPUPPETAnalysis=True', 'applyJECtoPuppiJets=False', 'applyZSelections=True', 'dropAnalyzerDumpEDM=False','applyWSelections=False']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V8', 'isMC=True', 'runMVAPUPPETAnalysis=True', 'applyJECtoPuppiJets=False', 'applyZSelections=True', 'dropAnalyzerDumpEDM=False','applyWSelections=False']

config.JobType.inputFiles  = []
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-AsymptNoPURawReco_MCRUN2_74_V9A-v4/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-AsymptFlat10to50bx25Raw_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM'

config.Data.unitsPerJob = 1
config.Data.inputDBS  = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 40000
config.Data.publication = False
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A'
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_StartupFlat10to50bx50Raw_MCRUN2_74_V8'
#config.Data.outLFNDirBase  = '/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_AsymptNoPURawReco_MCRUN2_74_V9A'
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_AsymptFlat10to50bx25Raw_MCRUN2_74_V9'
config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
