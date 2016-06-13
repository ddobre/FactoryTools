import ROOT
# Workaround to fix threadlock issues with GUI
ROOT.PyConfig.StartGuiThread = False
import logging
logging.basicConfig(level=logging.INFO)

import commonOptions
import jetCalibConfig

parser = commonOptions.parseCommonOptions()
#you can add additional options here if you want
#parser.add_option('--verbosity', help   = "Run all algs at the selected verbosity.",choices=("info", "warning","error", "debug", "verbose"), default="error")
parser.add_option('--truth', help = "Specify running on truth jets", action="store_true", dest="isTruth")

(options, args) = parser.parse_args()
#print options

ROOT.gROOT.Macro( '$ROOTCOREDIR/scripts/load_packages.C' )

# create a new sample handler to describe the data files we use
logging.info("creating new sample handler")
sh_all = ROOT.SH.SampleHandler()

commonOptions.fillSampleHandler(sh_all, options.inputDS)

sh_all.setMetaString ("nc_tree", "CollectionTree");
#sh_all.printContent();

# this is the basic description of our job
logging.info("creating new job")
job = ROOT.EL.Job()
job.sampleHandler(sh_all)
job.useXAOD()

logging.info("creating algorithms")

outputFilename = "trees"
output = ROOT.EL.OutputStream(outputFilename);

#here we add the algorithms we want to run over
import collections
algsToRun = collections.OrderedDict()

# Basic event selection
algsToRun["basicEventSelection"]       = ROOT.BasicEventSelection()
commonOptions.configxAODAnaHelperAlg(algsToRun["basicEventSelection"] )
if options.isTruth :
  setattr(algsToRun["basicEventSelection"],"m_useMetaData",False)
  setattr(algsToRun["basicEventSelection"],"m_truthLevelOnly",True)

# Jet calibration and selection
if not options.isTruth :
  jetCalibDict = jetCalibConfig.jetCalibrationDict
  algsToRun["calibrateJets"]           = ROOT.JetCalibrator()
  commonOptions.configxAODAnaHelperAlg(algsToRun["calibrateJets"],jetCalibDict)

  jetSelectionDict = jetCalibConfig.jetSelectionDict

else :
  jetSelectionDict = jetCalibConfig.jetSelectionDict_truth

print jetSelectionDict
algsToRun["selectJets"]              = ROOT.JetSelector()
commonOptions.configxAODAnaHelperAlg(algsToRun["selectJets"],jetSelectionDict)

algsToRun["preselectDileptonicWW"]   = ROOT.PreselectDileptonicWWEvents()#todo change this if we need it
algsToRun["selectNixonResolved"]        = ROOT.SelectNixonResolvedEvents()
algsToRun["postselectDileptonicWW"]    = ROOT.PostselectDileptonicWWEvents()

algsToRun["calculateRegionVars"]                      = ROOT.CalculateRegionVars()
algsToRun["calculateRegionVars"].calculatorName       = ROOT.CalculateRegionVars.nrCalculator
# added this

print "=============================="
print options.isTruth

if options.isTruth :
  algsToRun["calculateRegionVars"].isTruth              = 1

#ROOT.CalculateRegionVars.isTruth
############

regionName = "SR"
tmpWriteOutputNtuple                       = ROOT.WriteOutputNtuple()
tmpWriteOutputNtuple.outputName            = outputFilename
tmpWriteOutputNtuple.regionName            = regionName
tmpWriteOutputNtuple.systName            = ""
algsToRun["writeOutputNtuple_"+regionName] = tmpWriteOutputNtuple

if options.doSystematics : commonOptions.doSystematics(algsToRun)

job.outputAdd(output);
commonOptions.addAlgsFromDict(job , algsToRun , options.verbosity)

if options.nevents > 0 :
    logging.info("Running " + str(options.nevents) + " events")
    job.options().setDouble (ROOT.EL.Job.optMaxEvents, float(options.nevents));

commonOptions.overwriteSubmitDir(options.submitDir , options.doOverwrite)
commonOptions.submitJob         ( job , options.driver , options.submitDir , options.gridUser , options.gridTag)
