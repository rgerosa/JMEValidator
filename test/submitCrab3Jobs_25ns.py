from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.requestName = ''

## MC 
config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.General.workArea = 'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Asympt25ns_MCRUN2_74_V9'
#config.General.workArea = 'TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9'
#config.General.workArea = 'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Asympt25ns_MCRUN2_74_V9'
#config.General.workArea = 'TT_TuneCUETP8M1_13TeV-powheg-pythia8_Asympt25ns_MCRUN2_74_V9-v2'

#config.General.workArea = 'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
#config.General.workArea = 'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.General.workArea = 'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.General.workArea = 'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.General.workArea = 'WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.General.workArea = 'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.General.workArea = 'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.General.workArea = 'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'

#config.General.workArea = 'WWTo2L2Nu_13TeV-powheg_Asympt25ns_MCRUN2_74_V9'
#config.General.workArea = 'WZ_TuneCUETP8M1_13TeV-pythia8_Asympt25ns_MCRUN2_74_V9'
#config.General.workArea = 'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Asympt25ns_MCRUN2_74_V9'

## DATA

config.section_('JobType')
config.JobType.psetName    = 'runFrameworkMC.py'
config.JobType.pluginName  = 'Analysis'

## MC
config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9']

#config.JobType.inputFiles  = ['PY8_RunIISpring15DR74_bx50_MC_PFCHS.db','PY8_RunIISpring15DR74_bx50_MC_Puppi.db']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')

## MC
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-AsymptNoPURawReco_MCRUN2_74_V9A-v4/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-AsymptFlat10to50bx25Raw_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset = '/WWTo2L2Nu_13TeV-powheg/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'                                                        
#config.Data.inputDataset = '/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'                                             

#config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'


config.Data.unitsPerJob = 1
config.Data.inputDBS  = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 40000
config.Data.publication = False

#MC
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9'
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Asympt25ns_MCRUN2_74_V9'
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Asympt25ns_MCRUN2_74_V9'
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/TT_TuneCUETP8M1_13TeV-powheg-pythia8_Asympt25ns_MCRUN2_74_V9-v2'
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/WWTo2L2Nu_13TeV-powheg_Asympt25ns_MCRUN2_74_V9'
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/WZ_TuneCUETP8M1_13TeV-pythia8_Asympt25ns_MCRUN2_74_V9'
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Asympt25ns_MCRUN2_74_V9'

#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'

# DATA
#config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/Run2015B_DoubleMuon_PromptReco_puppiJEC'

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_DE_DESY'
config.Data.outLFNDirBase = '/store/user/rfriese/mvamet/skimming/2015-11-09/'
config.General.workArea = '/nfs/dust/cms/user/rfriese/crab_mvamet_skim-2015-11-09'


#  LocalWords:  MINIAODSIM
