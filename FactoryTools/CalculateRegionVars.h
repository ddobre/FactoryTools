#ifndef FactoryTools_CalculateRegionVars_H
#define FactoryTools_CalculateRegionVars_H

#include <EventLoop/Algorithm.h>

class RegionVarCalculator;

class CalculateRegionVars : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;

  enum RegionCalculatorName {
    none           = 0,
    lvlvCalculator = 1,
    zlCalculator   = 2,
    tlsCalculator  = 3,
    b4jCalculator  = 4,
    nrCalculator   = 5
  };


  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  RegionCalculatorName  calculatorName;
  int isTruth = 0;

  // this is a standard constructor
  CalculateRegionVars ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

private :
  RegionVarCalculator  * m_calculator;//!

public :
  // this is needed to distribute the algorithm to the workers
  ClassDef(CalculateRegionVars, 1);
};

#endif
