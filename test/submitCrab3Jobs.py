from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.requestName = ''
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_0p0'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_1p0'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_1p5'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_2p0'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_2p5'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_3p0'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_3p5'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_4p0'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_5p0'

#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_Muon_JEC'
#config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_jetPt20'
config.General.workArea = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_jetPt30'

config.section_('JobType')
config.JobType.psetName    = 'runFrameworkMC.py'
config.JobType.pluginName  = 'Analysis'


#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=False','useJECFromDB=True','applyJECtoPuppiJets=True']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','jetPtCut=20']
config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','jetPtCut=30']

#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,0,0','ptNeutralCutSlope=0.015,0,0']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,1,1','ptNeutralCutSlope=0.015,0,0']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,1.5,1.5','ptNeutralCutSlope=0.015,0,0']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,2,2','ptNeutralCutSlope=0.015,0,0']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,2.5,2.5','ptNeutralCutSlope=0.015,0,0']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,3,3','ptNeutralCutSlope=0.015,0,0']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,3.5,3.5','ptNeutralCutSlope=0.015,0,0']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,4,4','ptNeutralCutSlope=0.015,0,0']
#config.JobType.pyCfgParams = ['globalTag=MCRUN2_74_V9A','runPuppiNoMuon=True','useJECFromDB=True','applyJECtoPuppiJets=True','ptNeutralCut=0.1,5,5','ptNeutralCutSlope=0.015,0,0']

config.JobType.inputFiles  = ['PY8_RunIISpring15DR74_bx50_MC_PFCHS.db','PY8_RunIISpring15DR74_bx50_MC_Puppi.db']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-AsymptNoPURawReco_MCRUN2_74_V9A-v4/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-AsymptFlat10to50bx25Raw_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM'

config.Data.unitsPerJob = 1
config.Data.inputDBS  = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 40000
config.Data.publication = False
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_0p0'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_1p0'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_1p5'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_2p0'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_2p5'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_3p0'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_3p5'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_4p0'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_ptFw_centralfx_5p0'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_Muon_JEC'
#config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_jetPt20'
config.Data.outLFNDirBase = '/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_noMuon_JEC_jetPt30'

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
