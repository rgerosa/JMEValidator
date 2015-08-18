from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.requestName = ''

## MC 
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A'
#config.General.workArea = 'WWTo2L2Nu_13TeV-powheg_Asympt50ns_MCRUN2_74_V9A'
#config.General.workArea = 'WZ_TuneCUETP8M1_13TeV-pythia8_Asympt50ns_MCRUN2_74_V9A'
#config.General.workArea = 'ZZ_TuneCUETP8M1_13TeV-pythia8_Asympt50ns_MCRUN2_74_V9A'
#config.General.workArea = 'TTTo2L2Nu_13TeV-powheg_Asympt50ns_MCRUN2_74_V9A-v2'

## DATA
#config.General.workArea = 'DoubleMuon_Run2015B-17Jul2015-v1'
config.General.workArea = 'DoubleMuon_Run2015B-PromptReco-v1'

config.section_('JobType')
config.JobType.psetName    = 'runFrameworkMC.py'
config.JobType.pluginName  = 'Analysis'

## MC
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,1.7,2.0','etaCutForMetDiagnostic=3.0']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,1.7,2.0']

## DATA
config.JobType.pyCfgParams = ['globalTag=74X_dataRun2_Prompt_v1','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,1.7,2.0','isMC=False','etaCutForMetDiagnostic=3.0']
#config.JobType.pyCfgParams = ['globalTag=74X_dataRun2_Prompt_v1','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,1.7,2.0','isMC=False']

config.JobType.inputFiles  = ['Summer15_50nsV2_DATA.db']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')

## MC
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
#config.Data.inputDataset = '/WWTo2L2Nu_13TeV-powheg/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
#config.Data.inputDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
#config.Data.inputDataset = '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
#config.Data.inputDataset = '/TTTo2L2Nu_13TeV-powheg/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'

## DATA
#config.Data.inputDataset = '/DoubleMuon/Run2015B-17Jul2015-v1/MINIAOD'
config.Data.inputDataset = '/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD'

#config.Data.runRange = '246908-251562'
config.Data.runRange = '251563-251883'

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'


#config.Data.unitsPerJob = 1
config.Data.inputDBS  = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 40000
config.Data.publication = False

#MC
config.Data.outLFNDirBase = '/store/user/rgerosa/PUPPETAnalysis/FixedEtaCut/'

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

#  LocalWords:  MINIAODSIM
