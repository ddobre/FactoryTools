// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "SUSYTools/SUSYObjDef_xAOD.h"

#include <FactoryTools/RegionVarCalculator_lvlv.h>
#include <FactoryTools/RegionVarCalculator_tls.h>
#include <FactoryTools/RegionVarCalculator_zl.h>
#include <FactoryTools/RegionVarCalculator_b4j.h>
#include <FactoryTools/RegionVarCalculator_nr.h>
#include <FactoryTools/CalculateRegionVars.h>

#include <FactoryTools/printDebug.h>
#include <FactoryTools/strongErrorCheck.h>



// this is needed to distribute the algorithm to the workers
ClassImp(CalculateRegionVars)



CalculateRegionVars :: CalculateRegionVars () :
calculatorName(none),//user needs to choose their calculator name
m_calculator(nullptr)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode CalculateRegionVars :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CalculateRegionVars :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CalculateRegionVars :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CalculateRegionVars :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CalculateRegionVars :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  ATH_MSG_INFO("You have configured a " << calculatorName << " calculator.  See the code for enum definitions. ");
  STRONG_CHECK( calculatorName != none);

  if(calculatorName == lvlvCalculator)
    {
      m_calculator = new RegionVarCalculator_lvlv;
      STRONG_CHECK_SC( m_calculator->initialize(wk()));
    }
  else if(calculatorName == zlCalculator){
    m_calculator  = new RegionVarCalculator_zl;
    STRONG_CHECK_SC( m_calculator->initialize(wk()));
  }
  else if(calculatorName == tlsCalculator){
    m_calculator  = new RegionVarCalculator_tls;
    STRONG_CHECK_SC( m_calculator->initialize(wk()));
  }  
  else if(calculatorName == b4jCalculator){
    m_calculator  = new RegionVarCalculator_b4j;
    STRONG_CHECK_SC( m_calculator->initialize(wk()));
  }
  else if(calculatorName == nrCalculator){
    m_calculator  = new RegionVarCalculator_nr;
    if(isTruth == 1) {
      m_calculator->isTruth = true;
    }
    std::cout << "-------------------- Printing isTruth from CalculateRegionVars.cxx: " << isTruth << std::endl;
    STRONG_CHECK_SC( m_calculator->initialize(wk()));
  }
  else {
    ATH_MSG_ERROR( "You failed to provide a proper calculator.  If you have created a new one, make sure to add it to the : " << __PRETTY_FUNCTION__ << " algorithm.");
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CalculateRegionVars :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  printDebug();

  xAOD::TStore * store = wk()->xaodStore();

  std::map<std::string,double>               * mymap    = new std::map<std::string,double>;
  std::map<std::string,std::vector<double> > * myvecmap = new std::map<std::string,std::vector<double> >;

  printDebug();
  STRONG_CHECK_SC( m_calculator->calculate(*mymap, *myvecmap   ));
  STRONG_CHECK   ( store->record( mymap    , "RegionVarsMap"   ));
  STRONG_CHECK   ( store->record( myvecmap , "VecRegionVarsMap"));
  printDebug();
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CalculateRegionVars :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CalculateRegionVars :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode CalculateRegionVars :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}
