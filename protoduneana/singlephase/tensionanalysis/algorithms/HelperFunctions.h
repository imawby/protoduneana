// Class header for HelperFunctions class

#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "protoduneana/singlephase/stoppingmuonfilter/algorithms/Structs.h"

namespace util{

  class HelperFunctions{

    public:

      // default constructor
      HelperFunctions();

      // default destructor
      ~HelperFunctions(){}

      // configure HelperFunctions with fhicl parameters
      double CalculateDistance(std::vector< double > vectorA, std::vector< double > vectorB);

  };

}

#endif
