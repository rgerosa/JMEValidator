import FWCore.ParameterSet.Config as cms

# Common parameters used in all modules
CommonParameters = cms.PSet(
    # record flavor information, consider both RefPt and JetPt
    doComposition   = cms.bool(True),
    doFlavor        = cms.bool(True),
    doRefPt         = cms.bool(True),
    doJetPt         = cms.bool(True),
    # MATCHING MODE: deltaR(ref,jet)
    deltaRMax       = cms.double(99.9),
    # deltaR(ref,parton) IF doFlavor is True
    deltaRPartonMax = cms.double(0.25),
    # consider all matched references
    nJetMax         = cms.uint32(0),
)
 
process = cms.Process("JRA")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.GlobalTag.globaltag = "PHYS14_25_V2::All"


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
                                    
inputFiles = cms.untracked.vstring(
        # 'root://cmsxrootd.fnal.gov//store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/1020E374-B26B-E411-8F91-E0CB4E29C513.root'
       'file:///tmp/0432E62A-7A6C-E411-87BB-002590DB92A8.root'
#        'root://cmsxrootd.fnal.gov//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root',
        # 'root://cmsxrootd.fnal.gov//store/mc/Phys14DR/Neutrino_Pt-2to20_gun/MINIAODSIM/AVE20BX25_tsg_PHYS14_25_V3-v1/00000/007B75A0-538F-E411-B87D-00259059649C.root',

       ##################
       # Mu-enrich QCD
       ##################
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/00000/7E53A298-4CA7-E411-8368-002481E154CE.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/00000/F263724C-3FA7-E411-A3BE-002590A371AE.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/325BE5B9-AAA6-E411-8371-001E673972E2.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/381E5CF2-8BA7-E411-B4ED-0025B3E05C2C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/42E3F46B-ABA6-E411-AD64-002590A887AC.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/50474FE6-8BA7-E411-AE1C-001E67397021.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/585CB6CF-ABA6-E411-96B9-001E673972E2.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/760D3B77-ABA6-E411-8F44-0025B31E3C00.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/84F27E18-ABA6-E411-A9C1-002590200938.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/D694B317-ABA6-E411-9D5B-0025B3E0648A.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/DC83F868-ACA6-E411-9C36-0025B3E0648A.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/E0DAB3C2-AAA6-E411-8BEB-002590200A7C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/10000/E4657DE7-8BA7-E411-88FD-D8D385FF35C4.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/0619CB32-D4A5-E411-B1D9-001E673969FF.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/08FDDEFE-A3A6-E411-BEB0-0025B3E063A0.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/1251DA45-A4A6-E411-999C-002590A3C978.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/1659F588-D1A5-E411-99C2-0025B3E063A8.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/1E152064-D5A5-E411-AC59-002590A80E1E.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/229AB71A-A3A6-E411-AEE2-002590A3C978.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/246D4865-D5A5-E411-8297-002481E153FA.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/2A4AE265-D5A5-E411-B3CA-0025B31E3D3C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/2C44C49D-A5A6-E411-AE53-0025B3E066A4.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/42040562-D5A5-E411-BA45-001E673972E7.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/4C98D16F-D5A5-E411-A134-0025B31E3D3C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/52880732-ECA5-E411-AEC8-001E6739723D.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/5A702077-D0A5-E411-9E9B-001E67396A63.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/76B0334C-A4A6-E411-B8C9-0025B3E05E00.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/7CB5EBFC-A3A6-E411-B109-0025B3E0638E.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/7EB89680-D0A5-E411-A666-002590A887AC.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/82F32666-D5A5-E411-9CE4-0025B3E06584.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/841896FA-D2A5-E411-8924-002590200B70.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/888EEB64-D5A5-E411-B0E9-0025B3E0657E.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/88CACCE3-D4A5-E411-B947-001E673972E7.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/96943E69-D3A5-E411-BFD7-001E67397F17.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/9EC3AA23-D5A5-E411-84F8-002590A831B6.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/A6A6263E-D1A5-E411-ABAD-001E67397454.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/AE7F4C45-A4A6-E411-81A2-002590A887EE.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/C2180236-A6A6-E411-99D9-0025902008D0.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/D00BA1B6-57A7-E411-B3C9-002590A80E1E.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/D05EFF0F-D6A5-E411-95D2-001E67397003.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/D26D4965-D5A5-E411-86DA-001E67398390.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/D2E4C41C-A3A6-E411-AC7D-0025B3E05C2C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/D8FF9D2D-D5A5-E411-BF46-001E67397F2B.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/E83EEF03-A4A6-E411-9107-D8D385FF4B32.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/20000/F0EA9D1A-A3A6-E411-A468-002590A887EE.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/0637D559-57A6-E411-8AB8-001E67396E05.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/107ED920-49A7-E411-94E0-0025B3E06388.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/1875718C-89A5-E411-84CB-002590200808.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/3AA3D149-57A6-E411-8307-001E67396DEC.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/483E68D2-57A6-E411-AAFF-002590FC5AD0.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/5E0913A9-54A6-E411-ABD7-001E67397233.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/84319A7E-59A6-E411-9318-001E67396D4C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/88277448-57A6-E411-A733-001E67397CC9.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/A055CD85-57A6-E411-BF1A-002590200B7C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/A40CA27D-8AA5-E411-9E80-0025902009E4.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/A4F467CF-98A6-E411-8F3E-002590A88802.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/A6BD4F1A-55A6-E411-92CA-002590A3C970.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/AA712A62-8AA5-E411-A3AF-001E67398CAA.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/AEBFB372-8AA5-E411-896B-001E6739801B.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/B42875F4-54A6-E411-BB29-002590A8881E.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/B4A288BE-8AA5-E411-BCB3-001E67396577.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/BC43297F-57A6-E411-A4B9-002590A37116.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/BC795964-8AA5-E411-82DF-00259020097C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/C28D0C17-55A6-E411-969A-002590A3C954.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/C66396CE-8AA5-E411-8552-001E673973D2.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/CE1738CB-8AA5-E411-9772-002590A80E1C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/D4927EA7-54A6-E411-8CDD-002590A80E1C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/D8658912-8BA5-E411-9FD0-001E67398C1E.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/DA7D8F69-59A6-E411-9E60-001E67396775.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/E2C7BE8A-57A6-E411-9106-002590A371AE.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/E2DAAC02-8AA5-E411-A0F9-002590200B68.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/E49C5416-55A6-E411-B3FD-001E67396A09.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/F8D492F7-57A6-E411-AED2-002590A88812.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/30000/FC446EBC-8AA5-E411-BD1A-002590A8882C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/02495B78-27A6-E411-B511-001E6739801B.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/087D547C-27A6-E411-AD45-001E673972CE.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/0A4D98C4-28A6-E411-92E4-0025902008F0.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/10BDAF5D-2AA6-E411-9233-002590200828.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/163D7FAE-29A6-E411-877D-002481E14E00.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/20080E21-48A7-E411-B262-002590A3CA1A.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/20C1BE5E-29A6-E411-AE37-002590200B00.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/2AC87EC7-2BA6-E411-A48D-001E6739730A.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/2C5BFB58-29A6-E411-B531-001E6739689C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/3000EDAE-29A6-E411-8832-002590A4FFB8.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/36A91FA2-31A6-E411-90FB-001E67398CA0.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/38285975-2BA6-E411-94AE-002481E14FEE.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/44B991AB-29A6-E411-AD8D-002590A4C6BE.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/50CE5232-FCA6-E411-A085-D8D385FF4A72.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/6CB4CA7F-2BA6-E411-AC49-0025B3E063A0.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/6CD7EAAD-28A6-E411-94BC-002590A88736.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/708F1677-27A6-E411-BF43-0025B3E05C60.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/72BA1078-27A6-E411-9415-001E6739811A.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/7E3FAE5F-2AA6-E411-8BF5-0025902008F0.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/803568B3-29A6-E411-91E9-002590200A18.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/86829A32-FCA6-E411-AA86-0025B3E05E0A.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/8C8BB4DF-2BA6-E411-86F9-D8D385FF778C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/AE317AC3-28A6-E411-A891-001E6739811A.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/B077A3D8-2EA6-E411-89E9-002590A8880E.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/B40B0947-50A5-E411-91BD-0025902009CC.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/B4CD8B64-29A6-E411-A334-002590200B6C.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/B63BBFBC-28A6-E411-9783-0025B3E05C60.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/B85E6ACB-2BA6-E411-95D2-002481E14E14.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/C4D94CBA-28A6-E411-9193-002590FC5ACC.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/D2113BB2-29A6-E411-B539-001E67397391.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/D8E75B76-27A6-E411-BED2-001E6739825A.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/DA4986E2-50A5-E411-B6B0-001E67398CE1.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/E01980AC-29A6-E411-B9BE-002590200948.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/EA339FAF-29A6-E411-9403-002481E14E56.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/EAA9ECB1-2CA6-E411-9E9D-002481E14FEE.root',
#  	'/store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v3/40000/EE0B41AB-28A6-E411-9000-002590A4C6BE.root'
       ##################
       #  DY
       ##################
#	 '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root',#
#	 '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06C61714-7E6C-E411-9205-002590DB92A8.root',


    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource", fileNames = inputFiles )

# Services
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName = cms.string('output.root')

process.out = cms.OutputModule("PoolOutputModule",
        outputCommands  = cms.untracked.vstring(),
        fileName       = cms.untracked.string("output_edm.root")
        )



### load default PAT sequence
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
process.load("JMEAnalysis.JMEValidator.convertPackedCandToRecoCand_cff")

process.packedPFCandidatesWoMuon  = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV>=2 && abs(pdgId)!=13 " ) )
convertedPackedPFCandidatesWoMuon = cms.EDProducer('convertCandToRecoCand',
						   src = cms.InputTag('packedPFCandidatesWoMuon')
						   )

setattr( process, 'convertedPackedPFCandidatesWoMuon', convertedPackedPFCandidatesWoMuon )
process.patseq = cms.Sequence(process.convertedPackedPFCandidates *
			      convertedPackedPFCandidatesWoMuon *
			      process.patCandidates * process.selectedPatCandidates)
process.p = cms.Path(process.patseq)

# change the input collections
process.particleFlowPtrs.src = 'convertedPackedPFCandidates'
process.pfPileUpIso.Vertices = 'offlineSlimmedPrimaryVertices'
process.pfPileUp.Vertices    = 'offlineSlimmedPrimaryVertices'

# remove unnecessary PAT modules
process.p.remove(process.makePatElectrons)
process.p.remove(process.makePatPhotons)
process.p.remove(process.makePatJets)
process.p.remove(process.makePatTaus)
process.p.remove(process.makePatMETs)
process.p.remove(process.patCandidateSummary)
process.p.remove(process.selectedPatElectrons)
process.p.remove(process.selectedPatPhotons)
process.p.remove(process.selectedPatJets)
process.p.remove(process.selectedPatTaus)
process.p.remove(process.selectedPatCandidateSummary)

### muon selection
process.selectedPatMuons.src = 'slimmedMuons'
process.selectedPatMuons.cut = 'pt>10 && abs(eta)<2.4'

# load user-defined particle collections (e.g. PUPPI)

# -- PF-Weighted
process.load('CommonTools.ParticleFlow.deltaBetaWeights_cff')

# -- PUPPI
from JMEAnalysis.JMEValidator.pfPUPPISequence_cff import *
load_pfPUPPI_sequence(process, 'pfPUPPISequence', algo = 'PUPPI',
  src_puppi = 'pfAllHadronsAndPhotonsForPUPPI',
  cone_puppi_central = 0.5
)

# change the input collections
process.pfAllHadronsAndPhotonsForPUPPI.src = 'convertedPackedPFCandidates'
process.particleFlowPUPPI.candName = 'packedPFCandidates'
process.particleFlowPUPPI.vertexName = 'offlineSlimmedPrimaryVertices'

# -- PUPPI isolation calculation without muon
load_pfPUPPI_sequence(process, 'pfNoMuonPUPPISequence', algo = 'NoMuonPUPPI',
  src_puppi = 'pfAllHadronsAndPhotonsForNoMuonPUPPI',
  cone_puppi_central = 0.5
)
process.pfAllHadronsAndPhotonsForNoMuonPUPPI.src = 'convertedPackedPFCandidatesWoMuon'
process.particleFlowNoMuonPUPPI.candName         = 'packedPFCandidatesWoMuon'
process.particleFlowNoMuonPUPPI.vertexName       = 'offlineSlimmedPrimaryVertices'

from JMEAnalysis.JMEValidator.makePUPPIJets_cff import *
load_PUPPIJet_sequence(process,"PUPPIJetSequence",[0.4,0.8])

process.p.replace(
  process.pfParticleSelectionSequence,
  process.pfParticleSelectionSequence  *
  process.pfDeltaBetaWeightingSequence *
  process.pfPUPPISequence *
  process.pfNoMuonPUPPISequence *
  process.PUPPIJetSequence *
  process.patJetPartonMatchAK4PUPPIJets *
  process.patJetGenJetMatchAK4PUPPIJets *
  process.patJetsAK4PUPPIJets *
  process.patJetPartonMatchAK8PUPPIJets *
  process.patJetGenJetMatchAK8PUPPIJets *
  process.patJetsAK8PUPPIJets
)


from JMEAnalysis.JMEValidator.MuonPFIsolationSequence_cff import *
muon_src, cone_size = 'selectedPatMuons', 0.4

process.pfCHLVForIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV>=2 && abs(charge) > 0"))
process.pfCHPUForIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV<=1 && abs(charge) > 0"))
process.pfPhotonsForIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("pdgId==22"))
process.pfNHForIso = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("pdgId!=22 && abs(charge) == 0" ))

load_muonPFiso_sequence(process, 'MuonPFIsoSequenceSTAND', algo = 'R04STAND',
  src = muon_src,
  src_charged_hadron = 'pfCHLVForIso',
  src_neutral_hadron = 'pfNHForIso',
  src_photon         = 'pfPhotonsForIso',
  src_charged_pileup = 'pfCHPUForIso',
  coneR = cone_size
)

load_muonPFiso_sequence(process, 'MuonPFIsoSequencePFWGT', algo = 'R04PFWGT',
  src = muon_src,
  src_neutral_hadron = 'pfWeightedNeutralHadrons',
  src_photon         = 'pfWeightedPhotons',
  coneR = cone_size
)

load_muonPFiso_sequence(process, 'MuonPFIsoSequencePUPPI', algo = 'R04PUPPI',
  src = muon_src,
  src_charged_hadron = 'pfPUPPIChargedHadrons',
  src_neutral_hadron = 'pfPUPPINeutralHadrons',
  src_photon         = 'pfPUPPIPhotons',
  coneR = cone_size
)

load_muonPFiso_sequence(process, 'MuonPFIsoSequenceNoMuonPUPPI', algo = 'R04NOMUONPUPPI',
  src = muon_src,
  src_charged_hadron = 'pfNoMuonPUPPIChargedHadrons',
  src_neutral_hadron = 'pfNoMuonPUPPINeutralHadrons',
  src_photon         = 'pfNoMuonPUPPIPhotons',
  coneR = cone_size
)

process.muPFIsoDepositCharged.src = 'slimmedMuons'
process.muonMatch.src = 'slimmedMuons'
process.muonMatch.matched = 'prunedGenParticles'
process.patMuons.pvSrc = 'offlineSlimmedPrimaryVertices'
process.p.remove(process.patMuons)

process.MuonPFIsoSequences = cms.Sequence(
  process.MuonPFIsoSequenceSTAND *
  process.MuonPFIsoSequencePFWGT *
  process.MuonPFIsoSequencePUPPI *
  process.MuonPFIsoSequenceNoMuonPUPPI
)

process.p.replace(
  process.selectedPatMuons,
  process.selectedPatMuons *
  process.MuonPFIsoSequences
)



# Create all needed jets collections

# jetsCollections is a dictionnary containing all the informations needed for creating a new jet collection. The format used is :
#  "name": {
#      "algo": string ; the jet algorithm to use
#      "pu_methods" : array of strings ; which PU method to use
#      "jec_payloads" : array of strings ; which JEC payload to use for making the JEC. The size must match the size of pu_methods
#      "jec_levels" : array of strings ; which JEC levels to apply
#  }

jetsCollections = {
        'AK4': {
            'algo': 'ak4',
            'pu_methods': ['Puppi', 'CHS', ''],
            'jec_payloads': ['AK4PFPUPPI', 'AK4PFchs', 'AK4PF'],
            'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
            },

        'AK8': {
            'algo': 'ak8',
            'pu_methods': ['Puppi', 'CHS', ''],
            'jec_payloads': ['AK8PFPUPPI', 'AK8PFchs', 'AK8PF'],
            'jec_levels': ['L1FastJet', 'L2Relative', 'L3Absolute']
            },
        }

from JMEAnalysis.JetToolbox.jetToolbox_cff import *

for name, params in jetsCollections.items():
    for index, pu_method in enumerate(params['pu_methods']):
        # Add the jet collection
        jetToolbox(process, params['algo'], 'dummy', 'out', PUMethod = pu_method, JETCorrPayload = params['jec_payloads'][index], JETCorrLevels = params['jec_levels'])


# Configure the analyzers

process.jmfw_analyzers = cms.Sequence()

for name, params in jetsCollections.items():
    for index, pu_method in enumerate(params['pu_methods']):

        algo = params['algo'].upper()
        jetCollection = 'selectedPatJets%sPF%s' % (algo, pu_method)

        print('Adding analyzer for jets collection \'%s\'' % jetCollection)

        analyzer = cms.EDAnalyzer('JetMETAnalyzer',
                CommonParameters,
                JetCorLabel   = cms.string(params['jec_payloads'][index]),
                JetCorLevels  = cms.vstring(params['jec_levels']),
                srcJet        = cms.InputTag(jetCollection),
                srcRho        = cms.InputTag('fixedGridRhoAllFastjet'),
                srcVtx        = cms.InputTag('offlineSlimmedPrimaryVertices'),
                srcMuons      = cms.InputTag('selectedPatMuons')
                )

        setattr(process, 'jmfw_%s' % params['jec_payloads'][index], analyzer)
        process.jmfw_analyzers += analyzer

process.puppiReader = cms.EDAnalyzer("puppiAnalyzer",
                                        treeName = cms.string("puppiTree"),
										maxEvents = cms.int32(1000),
                                        nAlgos = cms.InputTag("puppi", "PuppiNAlgos", "JRA"),
                                        rawAlphas = cms.InputTag("puppi", "PuppiRawAlphas", "JRA"),
                                        alphas = cms.InputTag("puppi", "PuppiAlphas", "JRA"),
                                        alphasMed = cms.InputTag("puppi", "PuppiAlphasMed", "JRA"),
                                        alphasRms = cms.InputTag("puppi", "PuppiAlphasRms", "JRA"),
                                        packedPFCandidates = cms.InputTag("packedPFCandidates", "", "PAT")
									)

from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMetPuppi = pfMet.clone();
process.pfMetPuppi.src = cms.InputTag('puppi')

process.selectedMuonsForZ = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
   (pfIsolationR04().sumChargedHadronPt+
    max(0.,pfIsolationR04().sumNeutralHadronEt+
    pfIsolationR04().sumPhotonEt-
    0.50*pfIsolationR04().sumPUPt))/pt < 0.20 && 
    (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))

process.leptonsAndMET = cms.EDAnalyzer("LeptonsAndMETAnalyzer",
                                       srcIsoMuons = cms.InputTag("selectedMuonsForZ"),
                                       srcMET = cms.InputTag("slimmedMETs"),
                                       srcPUPPET = cms.InputTag("pfMetPuppi"),
                                       srcVtx            = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                       srcMuons          = cms.InputTag( muon_src ),
                                       srcVMCHSTAND      = cms.InputTag('muPFIsoValueCHR04STAND'),
                                       srcVMNHSTAND      = cms.InputTag('muPFIsoValueNHR04STAND'),
                                       srcVMPhSTAND      = cms.InputTag('muPFIsoValuePhR04STAND'),
                                       srcVMPUSTAND      = cms.InputTag('muPFIsoValuePUR04STAND'),
                                       srcVMCHPFWGT      = cms.InputTag('muPFIsoValueCHR04PFWGT'),
                                       srcVMNHPFWGT      = cms.InputTag('muPFIsoValueNHR04PFWGT'),
                                       srcVMPhPFWGT      = cms.InputTag('muPFIsoValuePhR04PFWGT'),
                                       srcVMCHPUPPI      = cms.InputTag('muPFIsoValueCHR04PUPPI'),
                                       srcVMNHPUPPI      = cms.InputTag('muPFIsoValueNHR04PUPPI'),
                                       srcVMPhPUPPI      = cms.InputTag('muPFIsoValuePhR04PUPPI'),
                                       srcVMCHNOMUONPUPPI      = cms.InputTag('muPFIsoValueCHR04NOMUONPUPPI'),
                                       srcVMNHNOMUONPUPPI      = cms.InputTag('muPFIsoValueNHR04NOMUONPUPPI'),
                                       srcVMPhNOMUONPUPPI      = cms.InputTag('muPFIsoValuePhR04NOMUONPUPPI')
                                       )

process.p *= process.pfMetPuppi + process.selectedMuonsForZ + process.puppiReader + process.leptonsAndMET + process.jmfw_analyzers 


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

# schedule definition                                                                                                       
process.outpath  = cms.EndPath(process.out) 

#!
#! THAT'S ALL! CAN YOU BELIEVE IT? :-D
#!
