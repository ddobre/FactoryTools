#ifndef REGION_VARIABLE_CALCULATOR_nr_H
#define REGION_VARIABLE_CALCULATOR_nr_H
//author : Russell Smith
//date   : January 2015

#include "EventLoop/StatusCode.h"
#include "FactoryTools/RegionVarCalculator.h"

#include <map>
#include <iostream>

class RegionVarCalculator_nr : public RegionVarCalculator {

public :

  // bool isTruth = false;

private :
  //todo probably clean this up
  virtual EL::StatusCode doInitialize(EL::Worker * worker);
  virtual EL::StatusCode doCalculate (std::map<std::string, double>               & RegionVars,
				      std::map<std::string, std::vector<double> > & VecRegionVars);

  EL::StatusCode doAllCalculations   (std::map<std::string, double>               & RegionVars,
				      std::map<std::string, std::vector<double> > & VecRegionVars);

  EL::StatusCode doSRCalculations    (std::map<std::string, double>               & RegionVars,
				      std::map<std::string, std::vector<double> > & VecRegionVars);


public :
  // this is needed to distribute the algorithm to the workers
  ClassDef(RegionVarCalculator_nr, 1);

};

#endif //REGION_VARIABLE_CALCULATOR_nr_H

//  LocalWords:  ifndef
