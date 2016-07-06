import glob 
import sys,os
import subprocess
import ROOT 
from ROOT import TFile, gDirectory, TTree
from ROOT import TH1F, TH2F
from optparse import OptionParser
import pyAMI.client
import pyAMI.atlas.api as AtlasAPI

client = pyAMI.client.Client('atlas')

# Overwrite option
parser = OptionParser()
parser.add_option('--doOverwrite', help = "Overwrite output directory if it already exists" , action = "store_true" , default=False)
# Input directory option 
parser.add_option('--inputDS', help = "Specify input directory" , default= " ")
# Input file option 
parser.add_option('--inputFile', help = "Specify input File" , default= " ")
options, args = parser.parse_args()
# Useful for overwriting
output_directory_created = False

# ------------------------------------------------------- #

#                  Initialize parameters                  #

# ------------------------------------------------------- #

# Declare path of main directory:
path = "/cluster/warehouse/ddobre/grid_output/"
# print options.inputDS
path+= options.inputDS
input_root_file = options.inputFile 

variables = {
"JetLeadPairs_MassDiff_angleMethod":{
  "title":"Mass Difference Between Jet Pairs: Angle Method",
  "xaxis":"#Delta M [GeV]",
  "xmin":0,
  "xmax":13000,
  "bins":1000,
  },

"JetLeadPairs_MassDiff_invMMethod":{
  "title":"Mass Difference Between Jet Pairs: InvM Method",
  "xaxis":"#Delta M [GeV]",
  "xmin":0,
  "xmax":5000,
  "bins":1000,
  },

"JetLeadPairs_PairAngle_angleMethod":{
  "title":"Angle Between Jet Pairs: Angle Method",
  "xaxis":"#Delta Angle [Radians]",
  "xmin":0,
  "xmax":10,
  "bins":1000,
  },

"JetLeadPairs_PairAngle_invMMethod":{
  "title":"Angle Between Jet Pairs: InvM Method",
  "xaxis":"#Delta Angle [Radians]",
  "xmin":0,
  "xmax":10,
  "bins":1000,
  },

"JetLeadPairs_PtDiffInPair_angleMethod":{
  "title":"Pt Difference in Jet Pairs: Angle Method",
  "xaxis":"#Delta Pt [GeV]",
  "xmin":0,
  "xmax":7000,
  "bins":1000,
  },

"JetLeadPairs_PtDiffInPair_invMMethod":{
  "title":"Pt Difference in Jet Pairs: InvM Method",
  "xaxis":"#Delta Pt [Gev]",
  "xmin":0,
  "xmax":7000,
  "bins":1000,
  },

"pt5_pt4Ratio":{
  "title":"Pt5:Pt4 Ratio",
  "xaxis":"Pt_{5}:Pt_{4} Ratio [unitless]",
  "xmin":0,
  "xmax":1.1,
  "bins":1000,
  },

"jet5_DRNearestLead":{
  "title":"DeltaR of Jet 5 to Nearest of 4 Leading Jets",
  "xaxis":"DeltaR of Jet 5 to Nearest of Leading Jets [Radians]",
  "xmin":0,
  "xmax":4,
  "bins":500,
  },

"Jetlead4_XInvariantMass":{
  "title":"Invariant Mass of X particle",
  "xaxis":"Invariant Mass of X [GeV]",
  "xmin":0,
  "xmax":13000,
  "bins":1000,
  },

"Jetlead4_Ht":{
  "title":"Ht of Lead 4 Jets",
  "xaxis":"Ht [GeV]",
  "xmin":0,
  "xmax":13000,
  "bins":1000,
  },

"Jet_E":{
  "title":"Jet Energy",
  "xaxis":"E [GeV]",
  "xmin":0,
  "xmax":7000,
  "bins":1000,
  },

"Jet_Eta":{
  "title":"Jet Eta",
  "xaxis":"Eta [Radians]",
  "xmin":-3,
  "xmax":3,
  "bins":1000,
  },

"Jet_Phi":{
  "title":"Jet Phi",
  "xaxis":"Phi [Radians]",
  "xmin":-4,
  "xmax":4,
  "bins":1000,
  },

"Jet_Pt":{
  "title":"Jet Pt",
  "xaxis":"Pt [GeV]",
  "xmin":0,
  "xmax":7000,
  "bins":1000,
  }
}

sampleWeightDict = {
    "user.ddobre.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec":78420000000.0,
      "filtereff":0.97723,
      "events":10000},
    
    "user.ddobre.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec":78420000000.0,
      "filtereff":0.000657,
      "events":20000},

    "user.ddobre.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec":2433100000.0,
      "filtereff":0.000331,
      "events":10000},

    "user.ddobre.361023.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec":26454000.0,
      "filtereff":0.000322,
      "events":10000},

    "user.ddobre.361024.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec":254630.0,
      "filtereff":0.00053,
      "events":7993500},

    "user.ddobre.361025.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec":4553.5,
      "filtereff":0.000923,
      "events":7996500},

    "user.ddobre.361026.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ6W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec": 257.53,
      "filtereff":0.00094,
      "events":1997000},
    
    "user.ddobre.361027.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ7W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec":16.215,
      "filtereff":0.000393,
      "events":1996500},

    "user.ddobre.361028.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ8W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec":0.62502,
      "filtereff":0.010162,
      "events":2000000},

    "user.ddobre.361029.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ9W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
      "xsec":0.019639,
      "filtereff":0.012054,
      "events":2000000},

    "user.ddobre.361030.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ10W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
        "xsec":0.0011962,
        "filtereff":0.005894,
        "events":2000000},

    "user.ddobre.361031.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ11W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
        "xsec":4.2258e-05,
        "filtereff":0.002702,
        "events":1999000},

    "user.ddobre.361032.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ12W.29.06.2016_MemoryLeakTest_PythiaSamples_trees.root":{
        "xsec":1.0367e-06,
        "filtereff":0.000425,
        "events":1997500} 
}



# exit()

# ------------------------------------------------------- #

#                        Functions                        #

# ------------------------------------------------------- #

def makeWeightedROOT( root_file, full_path, variables, infoDict ):
  
  # -------------------------
  # Get the ROOT file 
  # -------------------------
  sample_path = os.path.join(full_path, root_file)
  sample = TFile( sample_path )
  sample_tree = sample.Get("trees_SR_;1")

  # -------------------------
  # Setup Event Loop
  # -------------------------
  mychain = gDirectory.Get('trees_SR_')
  entries = mychain.GetEntriesFast()
  
  # -------------------------
  # Setup output file
  # -------------------------
  weighted_file = ROOT.TFile( os.path.join( "weighted_root_files" , root_file ),"RECREATE" )
  print "Weighted root file created: " + os.path.join("weighted_root_files",root_file)
  # print "Took from full path: " + full_path 

  # -------------------------
  # Setup info dictionary
  # -------------------------
  ######### Via PyAMI #######
  # info = AtlasAPI.get_dataset_info(client,"mc15_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.merge.DAOD_EXOT2.e3569_s2576_s2132_r7725_r7676_p2669/" ) # full_path.split("/")[-1] )
  # # print info
  # print info
  # infoDict = {}
  # infoDict['crossSection'] = info['crossSection']
  #   # genFiltEff not found in the dictionary
  #   # infoDict['filterEff'] = info['genFiltEff']
  # infoDict['nEvt'] = info['totalEvents']
  # exit()
  #
  ######### Manual #########
  # print "Full path: " + full_path.split("/")[-1]
  # print "xsec: " , infoDict[full_path.split("/")[-1]]["xsec"]
  # print "filtereff: " , infoDict[full_path.split("/")[-1]]['filtereff']
  # print "events: " , infoDict[full_path.split("/")[-1]]['events']
  # section_weight = infoDict[full_path.split("/")[-1]]['xsec']*infoDict[full_path.split("/")[-1]]['filtereff'] # /infoDict[full_path.split("/")[-1]]['events']
  # print "section_weight value" , section_weight , " section_weight type: " , type(section_weight)
  # 
  # exit()
  #
  # -------------------------
  # Setup histogram dictionary
  # -------------------------
  histDict = {}
  var_names = [variable.GetName() for variable in sample_tree.GetListOfBranches()]
  for var_name in var_names:
    if var_name in variables.keys():
      print "Variable found in ROOT file: " + var_name + "... Creating histogram."
      histDict[var_name] = TH1F( var_name, var_name, variables[var_name]["bins"], variables[var_name]["xmin"] , variables[var_name]["xmax"]) 

  # -------------------------
  # Fill via eventloop 
  # -------------------------
  print "Histograms created. Now filling..."
  printcounter = 0

  for jentry in xrange( entries ):
      ientry = mychain.LoadTree( jentry )
      if ientry < 0: 
        break 

      nb = mychain.GetEntry( jentry ) 
      if nb <= 0: 
        break
      m_mcEventWeight = mychain.mcEventWeight
      
      # Fill defined variables into histograms
      for var_name in var_names:
        if var_name in variables.keys():
          container = "container_entry = mychain.{0}".format(var_name)
          exec container
          filling_hist = "histDict[var_name].Fill(mychain.{0})".format(var_name) # , m_mcEventWeight*section_weight ) ".format(var_name)

          # print "OVER HERE LOOK LOOK LOOK: mcEventWeight  " 
          # print m_mcEventWeight 
          # print " type - "
          # print type(m_mcEventWeight) 
          # print " xsection - " 
          # print infoDict['crossSection'] 
          # print " type - "  
          # print type(infoDict['crossSection'])
          # print float(infoDict['crossSection'])
          # print type(float(infoDict['crossSection']))

          if type(container_entry) is not float: 
            for item in container_entry:
              histDict[var_name].Fill(item) #, m_mcEventWeight*section_weight )
          else:
            exec filling_hist

          if printcounter % 2000000 == 0:
            print "Histogram " , root_file , "status: " , printcounter
          
          printcounter += 1

  weighted_file.Write()
  weighted_file.Close()
  
  print "Done filling histogram ---> " + root_file + " <--- filled with " , printcounter , " events." 

# ------------------------------------------------------- #

#                     Steering Code                       #

# ------------------------------------------------------- #

print "Accessing unweighted files; searching..."

# ------------------------------------------------------- #
#       Through an entire directory - not parallel        #
# ------------------------------------------------------- #
#
# for full_path, dirs, files in os.walk(path):
# 
#   if options.doOverwrite and output_directory_created is False:
#     delete_shit = "rm -r weighted_root_files"
#     subprocess.call(delete_shit,shell=True)
# 
#   if output_directory_created is False:
#     create_output_directory = "mkdir weighted_root_files"
#     subprocess.call(create_output_directory,shell=True)
#     output_directory_created = True 
#   
#   print "Accessing: " + full_path.split("/")[-1] 
# 
#   if "weighted_root_files" not in full_path.split("/")[-1]:
#     for root_file in files:
#       if ".root" in root_file:
#         makeWeightedROOT( root_file, full_path )
#
# ------------------------------------------------------- #
#         Via bash script - parallel per folder           #
# ------------------------------------------------------- #
# 
# for full_path, dirs, files in os.walk(path):
# 
#   # print "Accessing: " + full_path.split("/")[-1] 
# 
#   if "weighted_root_files" or "logs"  not in full_path.split("/")[-1]:
#     for root_file in files:
#       if ".root" in root_file:
#         print "Making a weighted copy of " +  root_file 
#         makeWeightedROOT( root_file, full_path )
# 
# ------------------------------------------------------- #
#            Via bash script - total parallel             #
# ------------------------------------------------------- # 

makeWeightedROOT(input_root_file, path, variables, sampleWeightDict)
