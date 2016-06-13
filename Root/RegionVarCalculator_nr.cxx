#include "EventLoop/StatusCode.h"
#include "EventLoop/Worker.h"
#include "xAODRootAccess/TStore.h"

#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODJet/JetAuxContainer.h"

#include "FactoryTools/RegionVarCalculator_nr.h"
#include "FactoryTools/strongErrorCheck.h"

#include <xAODAnaHelpers/HelperFunctions.h>

#include <cmath>


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

  if( !isTruth ) {
    const xAOD::VertexContainer* vertices = nullptr;
    STRONG_CHECK(event->retrieve( vertices, "PrimaryVertices"));
    RegionVars["NPV"] = HelperFunctions::countPrimaryVertices(vertices, 2);
  }
  
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

  // Jet Pairs

  std::vector<TLorentzVector> jetLead_TLor;
  
  // --- Invariant Mass  
  
  int position_invM[2] = { 0 , 5 };
  double pair_invM;
  std::vector<double> jetLeadPairs_invM;
  std::vector<TLorentzVector> jetLeadPairs;
  
  // --- Angles

  int position_angle[2] = { 0 , 5 };
  double pair_angle;
  std::vector<double> jetAnglePairs;
  std::vector<double> jetLeadPairs_angle;

  // jetHt 

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
    
    if( counter < 4 ) {
      jetLead_TLor.push_back( jet->p4() ); // I know this is unnecessary - will optimize after
      jetHt_lead += toGeV(jet->pt());
      jetHt_lead_invM += jet->p4();
    }
  counter ++;
  }

  // ---------- Calculations

  for( int i = 0; i < 3; i ++) { 
    for( int j = i + 1; j < 4; j++) {
      jetLeadPairs.push_back( jetLead_TLor[i] + jetLead_TLor[j] );
      jetAnglePairs.push_back( jetLead_TLor[i].DeltaR( jetLead_TLor[j] ) );
    } 
  } 

  //      ----------------------------------- 
  //      ----- Find jet pairs: Invariant mass method
  //      -----------------------------------

  /* 
  for( int i = 0; i < 6; i++ ) {
    std::cout << "Invariant mass of pair " << i << " = " << toGeV( jetLeadPairs[i].M() ) << std::endl;
  }  
  */

  pair_invM = toGeV( abs( (jetLeadPairs[0].M() - jetLeadPairs[5].M() ) ) );

  for( int i = 1; i < 3; i++) {
    
    if( toGeV( abs( (jetLeadPairs[i].M() - jetLeadPairs[5 - i].M() ) ) ) < pair_invM) {
      pair_invM = toGeV( abs( (jetLeadPairs[i].M() - jetLeadPairs[5 - i].M() ) ) );
      position_invM[0] = i;
      position_invM[1] = 5 - i;
    } 
  }
  
  // Now have the two jet pairs which have have the smallest invariant mass difference 
  // Can access them in the pairs array
  // Do something - lets shove the pair's invariant mass into a vector...
  
  jetLeadPairs_invM.push_back( toGeV( jetLeadPairs[ position_invM[ 0 ] ].M() ) );
  jetLeadPairs_invM.push_back( toGeV( jetLeadPairs[ position_invM[ 1 ] ].M() ) );

  //      ----------------------------------
  //      ----- Find jet pairs: Angle method 
  //      ----------------------------------

  pair_angle = ( jetAnglePairs[ 0 ] + jetAnglePairs[ 5 ] )/2;

  for( int i = 1; i < 3; i++) {
    
    if( (jetAnglePairs[ i ] + jetAnglePairs[ 5 - i ])/2 < pair_angle ) {
      pair_angle = (jetAnglePairs[ i ] + jetAnglePairs[ 5 - i ])/2;
      position_angle[0] = i;
      position_angle[1] = 5 - i;
    } 
  }
  
  // Now have the pairs from the the angle algorithm 
  // Can access them in the pairs array
  // Do something - we can toss their invariant mass into a vector and compare with the above one

  jetLeadPairs_angle.push_back( toGeV( jetLeadPairs[ position_angle[ 0 ] ].M() ) );
  jetLeadPairs_angle.push_back( toGeV( jetLeadPairs[ position_angle[ 1 ] ].M() ) );

  // ---------- Finalize 

  // Counter 
  Htcounter.push_back(counter);
  
  VecRegionVars[ "Htcounter" ] = Htcounter;
    
  // Variables
  VecRegionVars[ "jetPt" ]  = jetPtVec;
  VecRegionVars[ "jetEta" ] = jetEtaVec;
  VecRegionVars[ "jetPhi" ] = jetPhiVec;
  VecRegionVars[ "jetE" ]   = jetEVec;

  // jetHt and leading 4 jetHt
  RegionVars[ "JetHt" ] = jetHt;
  RegionVars[ "JetHt_lead" ] = jetHt_lead; 

  jetHt = toGeV(jetHt_invM.M());
  jetHt_lead = toGeV(jetHt_lead_invM.M());
  
  RegionVars[ "JetHt_invM" ] = jetHt; 
  RegionVars[ "JetHt_lead_invM" ] = jetHt_lead;

  // Leading jet pairs
  VecRegionVars[ "JetLeadPairs_invM" ] = jetLeadPairs_invM; 
  VecRegionVars[ "JetLeadPairs_angle" ] = jetLeadPairs_angle;

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode RegionVarCalculator_nr::doSRCalculations(std::map<std::string, double>& RegionVars,
							  std::map<std::string, std::vector<double> > & VecRegionVars)
{/*todo*/return EL::StatusCode::SUCCESS;}

