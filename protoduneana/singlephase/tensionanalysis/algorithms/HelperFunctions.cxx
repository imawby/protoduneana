#ifndef HELPERFUNCTIONS_CXX
#define HELPERFUNCTIONS_CXX

#include "HelperFunctions.h"

namespace util{

  HelperFunctions::HelperFunctions()
  {
  }

  double HelperFunctions::CalculateDistance(std::vector< double > vectorA, std::vector< double > vectorB)
  {

    if (vectorA.size() != vectorB.size())
      throw std::logic_error("Trying to calculate distance bwtween two vectors of different sizes");

    double total = 0;

    for (size_t i = 0; i < vectorA.size(); i++){
      total += std::pow(vectorA.at(i)-vectorB.at(i),2);
    }

    total = std::sqrt(total);

    return total;

  }


}

#endif
