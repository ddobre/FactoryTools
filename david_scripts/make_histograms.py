import glob
import sys,os 
import subprocess
import ROOT
from ROOT import TFile, gDirectory
from ROOT import TH1D
from optparse import OptionParser

ROOT.gROOT.LoadMacro("/home/ajhorton/atlasstyle-00-03-04/AtlasStyle.C")
ROOT.SetAtlasStyle()

parser = OptionParser()
parser.add_option('--doOverwrite', help = "Overwrite output directory if it already exists" , action = "store_true" , default=False)
options, args = parser.parse_args()

# ------------------------------------------------------- #
#                                                         #
# ---------- Declare Paramaters of Script Here ---------- #
#                                                         #
# ------------------------------------------------------- #

# Declare path of main directory here:
path = "/cluster/warehouse/ddobre/Nixon/output_directory/" 

# Declare variables you want to make histograms of (of ROOT file)
# Accessign the keys: variables.keys()
# Accessing the name: variables[variable] (second thing)

variables = {
"JetLeadPairs_MassDiff_angleMethod":{
  "title":"Mass Difference Between Jet Pairs: Angle Method",
  "xaxis":"#Delta M [GeV]"},
"JetLeadPairs_MassDiff_invMMethod":{
  "title":"Mass Difference Between Jet Pairs: InvM Method",
  "xaxis":"#Delta M [GeV]"},
"JetLeadPairs_PairAngle_angleMethod":{
  "title":"Angle Between Jet Pairs: Angle Method",
  "xaxis":"#Delta Angle [Radians]"},
"JetLeadPairs_PairAngle_invMMethod":{
  "title":"Angle Between Jet Pairs: InvM Method",
  "xaxis":"#Delta Angle [Radians]"},
"JetLeadPairs_PtDiffInPair_angleMethod":{
  "title":"Pt Difference in Jet Pairs: Angle Method",
  "xaxis":"#Delta Pt [GeV]"},
"JetLeadPairs_PtDiffInPair_invMMethod":{
  "title":"Pt Difference in Jet Pairs: InvM Method",
  "xaxis":"#Delta Pt [Gev]"},
"pt5_pt4Ratio":{
  "title":"Pt5:Pt4 Ratio",
  "xaxis":"Pt_{5}:Pt_{4} Ratio [unitless]"},
"jet5_DRNearestLead":{
  "title":"DeltaR of Jet 5 to Nearest of 4 Leading Jets",
  "xaxis":"DeltaR of Jet 5 to Nearest of Leading Jets [Radians]"},
"Jetlead4_XInvariantMass":{
  "title":"Invariant Mass of X particle",
  "xaxis":"Invariant Mass of X [GeV]"},
"Jetlead4_Ht":{
  "title":"Ht of Lead 4 Jets",
  "xaxis":"Ht [GeV]"}
}



output_directory_created = False

# ------------------------------------------------------- #
#                                                         #
# -------------- Define Functions here ------------------ #
#                                                         #
# ------------------------------------------------------- #

def makeHistogram( full_path , d , myFile , variables ):
  sample = TFile( os.path.join(full_path, d , myFile) )
  mychain = gDirectory.Get( 'trees_SR_' ) 
  entries = mychain.GetEntriesFast()
  samplename = full_path.split("/")[-1]
  histograms = {}
  canvas = ROOT.TCanvas("canvas", "canvas" ,  174, 52, 700, 500);
  
  for variable in variables.keys():
    
    mychain.Draw( variable + ">>" + variable) 
    hist = ROOT.gROOT.FindObject(variable) # manipulate the histogram 
    hist.SetTitle( variables[ variable ]["title"] )
    hist.GetYaxis().SetTitle("Events")
    hist.GetXaxis().SetTitle(variables[ variable ]["xaxis"])
    hist.SetFillColor(9)
    # hist.SetLineColor(9)
    # ROOT.gPad.Update()
    # statbox = hist.FindObject( "stats")
    # statbox.SetX1NDC(0.7)     # new x start position
    # statbox.SetX2NDC(0.9)     # new x end position
    # statbox.SetY1NDC(0.7)     # new y start position
    # statbox.SetY2NDC(0.9)     # new y end position

    canvas.SaveAs("histogram_plots/" + samplename + "/" + variable + ".eps")
  
  return 

# ------------------------------------------------------- #
#                                                         #
# -------------------- Steering Code -------------------- #
#                                                         #
# ------------------------------------------------------- #

print "---------- HISTOGRAM SCRIPT STARTED ----------"

for full_path, dirs, files in os.walk(path):

  if options.doOverwrite and output_directory_created is False:
    delete_shit = "rm -r histogram_plots"
    subprocess.call(delete_shit,shell=True)

  if output_directory_created is False:
    create_output_directory = "mkdir histogram_plots"
    subprocess.call(create_output_directory,shell=True)
    output_directory_created = True

  for d in dirs: 
    if d == "data-trees":
      sub_d = os.listdir(full_path + "/" + d)
      for myFile in sub_d:
        if ".root" in myFile:
          create_sample_directory = "mkdir histogram_plots/" + full_path.split("/")[-1]
          print "------ CREATED NEW DIRECTORY: " + full_path.split("/")[-1]
          subprocess.call(create_sample_directory,shell=True)
          makeHistogram( full_path , d, myFile, variables )
