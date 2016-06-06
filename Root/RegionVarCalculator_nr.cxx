#include "EventLoop/StatusCode.h"
#include "EventLoop/Worker.h"
#include "xAODRootAccess/TStore.h"

#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODJet/JetAuxContainer.h"

#include "FactoryTools/RegionVarCalculator_nr.h"
#include "FactoryTools/strongErrorCheck.h"

#include <xAODAnaHelpers/HelperFunctions.h>


// this is needed to distribute the algorithm to the workers
ClassImp(RegionVarCalculator_nr)

EL::StatusCode RegionVarCalculator_nr::doInitialize(EL::Worker * worker) {
  if(m_worker != nullptr){
    std::cout << "You have called " << __PRETTY_FUNCTION__ << " more than once.  Exiting." << std::endl;
    return EL::StatusCode::FAILURE;
  }
  m_worker = worker;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RegionVarCalculator_nr::doCalculate(std::map<std::string, double>& RegionVars,
						     std::map<std::string, std::vector<double> >& VecRegionVars){

  xAOD::TStore * store = m_worker->xaodStore();//grab the store from the worker
  xAOD::TEvent* event = m_worker->xaodEvent();

  const xAOD::EventInfo* eventInfo = nullptr;
  STRONG_CHECK(event->retrieve( eventInfo, "EventInfo"));

  std::string const & regionName = eventInfo->auxdecor< std::string >("regionName");

  if      ( regionName.empty() ) {return EL::StatusCode::SUCCESS;}
  // If it hasn't been selected in any of the regions from any of the select algs, don't bother calculating anything...
  else if ( regionName == "SR" )  {return EL::StatusCode(doAllCalculations (RegionVars, VecRegionVars) == EL::StatusCode::SUCCESS &&
							 doSRCalculations  (RegionVars, VecRegionVars) == EL::StatusCode::SUCCESS);}

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RegionVarCalculator_nr::doAllCalculations(std::map<std::string, double>& RegionVars,
							 std::map<std::string, std::vector<double> > & VecRegionVars)
{/*todo*/
  xAOD::TStore * store = m_worker->xaodStore();
  xAOD::TEvent * event = m_worker->xaodEvent();

  doGeneralCalculations(RegionVars, VecRegionVars);

  // Get relevant info from the vertex container //////////////////////
  //

  const xAOD::VertexContainer* vertices = nullptr;
  STRONG_CHECK(event->retrieve( vertices, "PrimaryVertices"));
  RegionVars["NPV"] = HelperFunctions::countPrimaryVertices(vertices, 2);

  //
  /////////////////////////////////////////////////////////////////////

  auto toGeV = [](float a){return a*.001;};

  xAOD::IParticleContainer* jets_nominal(nullptr);
  STRONG_CHECK(store->retrieve(jets_nominal, "selectedJets"));

  // ---------- Initialize variables
  
  std::vector<double> Htcounter;

  std::vector<double> jetPtVec;
  std::vector<double> jetEtaVec;
  std::vector<double> jetPhiVec;
  std::vector<double> jetEVec;

  double jetHt = 0;
  double jetHt_lead = 0;

  TLorentzVector jetHt_invM;
  TLorentzVector jetHt_lead_invM;
  
  int counter = 0;

  // ---------- Loop over the container 

  for( const auto& jet : *jets_nominal) {
    jetPtVec.push_back( toGeV(jet->pt()));
    jetEtaVec.push_back( jet->p4().Eta() );
    jetPhiVec.push_back( jet->p4().Phi() );
    jetEVec.push_back( toGeV(jet->p4().E()) );

    jetHt += toGeV(jet->pt());
    jetHt_invM += jet->p4();
    // std::cout << "pt: " << toGeV(jet->pt()) << " counter: " << counter <<  std::endl;
    
    if( counter < 4 ) {
      jetHt_lead += toGeV(jet->pt());
      jetHt_lead_invM += jet->p4();
    }
  counter ++;
  }

  // ---------- Finalize 
  Htcounter.push_back(counter);
  
  VecRegionVars[ "Htcounter" ] = Htcounter;
    
  VecRegionVars[ "jetPt" ]  = jetPtVec;
  VecRegionVars[ "jetEta" ] = jetEtaVec;
  VecRegionVars[ "jetPhi" ] = jetPhiVec;
  VecRegionVars[ "jetE" ]   = jetEVec;

  RegionVars[ "JetHt" ] = jetHt;
  RegionVars[ "JetHt_lead" ] = jetHt_lead; 

  jetHt = toGeV(jetHt_invM.M());
  jetHt_lead = toGeV(jetHt_lead_invM.M());
  
  RegionVars[ "JetHt_invM" ] = jetHt; 
  RegionVars[ "JetHt_lead_invM" ] = jetHt_lead;
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode RegionVarCalculator_nr::doSRCalculations(std::map<std::string, double>& RegionVars,
							  std::map<std::string, std::vector<double> > & VecRegionVars)
{/*todo*/return EL::StatusCode::SUCCESS;}

