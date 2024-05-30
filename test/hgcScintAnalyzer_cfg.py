import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#parse command line arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('input', 
                 '/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/ttbar_D49_1120pre1_PU200_eolupdate_ter_20200713/GSD/',
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 "input directory")
options.register('adcThrMIP',
                 0.5, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.float, 
                 "threshold (in-time)")
options.register('scaleByDoseFactor',
                 1.0, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.float, 
                 "scale fluence by this factor")
options.register('adcThrMIPbxm1',
                 2.5, 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.float, 
                 "threshold (BX-1)")
options.register('fold',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "fold wafers histos x6/x3 for EE/HE")
options.register('geometry', 
                 'Extended2026D99', 
                 VarParsing.multiplicity.singleton, 
                 VarParsing.varType.string, 
                 'geometry to use')

# New inputs for HGCALMapping module 
# options.register('modules','UserCode/HGCElectronicsValidation/data/ModuleMaps/modulelocator_CEminus_V15p5.txt',mytype=VarParsing.varType.string,
#                  info="Path to module mapper. Absolute, or relative to CMSSW src directory")
# options.register('sicells','UserCode/HGCElectronicsValidation/data/CellMaps/WaferCellMapTraces.txt',mytype=VarParsing.varType.string,
#                  info="Path to Si cell mapper. Absolute, or relative to CMSSW src directory")
# options.register('sipmcells','UserCode/HGCElectronicsValidation/data/CellMaps/channels_sipmontile.hgcal.txt',mytype=VarParsing.varType.string,
#                  info="Path to SiPM-on-tile cell mapper. Absolute, or relative to CMSSW src directory")   
options.register('modules','Geometry/HGCalMapping/data/ModuleMaps/modulelocator_CEminus_V15p5.txt',mytype=VarParsing.varType.string,
                 info="Path to module mapper. Absolute, or relative to CMSSW src directory")
options.register('sicells','Geometry/HGCalMapping/data/CellMaps/WaferCellMapTraces.txt',mytype=VarParsing.varType.string,
                 info="Path to Si cell mapper. Absolute, or relative to CMSSW src directory")
options.register('sipmcells','Geometry/HGCalMapping/data/CellMaps/channels_sipmontile.hgcal.txt',mytype=VarParsing.varType.string,
                 info="Path to SiPM-on-tile cell mapper. Absolute, or relative to CMSSW src directory")                             

options.parseArguments()
print(f"options.modules = {options.modules}")   
print(f"options.sicells = {options.sicells}")   
print(f"options.sipmcells = {options.sipmcells}")         


#set geometry/global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.Geometry.Geometry%sReco_cff'%options.geometry)
process.load('Configuration.Geometry.Geometry%s_cff'%options.geometry)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 500

#ESSources/Producers for the logical mapping
#indexers
process.load('Geometry.HGCalMapping.hgCalMappingESProducer_cfi')
process.hgCalMappingESProducer.modules = cms.FileInPath(options.modules)
process.hgCalMappingESProducer.si = cms.FileInPath(options.sicells)
process.hgCalMappingESProducer.sipm = cms.FileInPath(options.sipmcells)
#cells and modules info
process.load('Configuration.StandardSequences.Accelerators_cff')
process.hgCalMappingCellESProducer = cms.ESProducer('hgcal::HGCalMappingCellESProducer@alpaka',
                                                      filelist=cms.vstring(options.sicells,options.sipmcells),
                                                      cellindexer=cms.ESInputTag('') )
process.hgCalMappingModuleESProducer = cms.ESProducer('hgcal::HGCalMappingModuleESProducer@alpaka',
                                                      filename=cms.FileInPath(options.modules),
                                                      moduleindexer=cms.ESInputTag('') )


#source to process
import os
if os.path.isdir(options.input):
    fList = ['file:'+os.path.join(options.input,f) for f in os.listdir(options.input) if '.root' in f]
else:
    fList = ['file:'+x for x in options.input.split(',')]

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(fList),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                        )

#number events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

#analyzer
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_ileakParam_toUse,HGCAL_cceParams_toUse
process.ana = cms.EDAnalyzer("HGCScintAnalyzer",
                            #  adcThrMIP       = cms.double(options.adcThrMIP),
                            #  adcThrMIPbxm1   = cms.double(options.adcThrMIPbxm1),
                            #  fold            = cms.bool(options.fold),
                            #  doseMap         = cms.string('SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-3.7.20.txt'),
                            #  scaleByDoseAlgo = cms.uint32(0),
                            #  scaleByDoseFactor = cms.double(options.scaleByDoseFactor),
                            #  ileakParam      = HGCAL_ileakParam_toUse,
                            #  cceParams       = HGCAL_cceParams_toUse
                        )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.output)
                               )

process.p = cms.Path(process.ana)
