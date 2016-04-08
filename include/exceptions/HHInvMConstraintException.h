#ifndef HHInvMConstraintException_
#define HHInvMConstraintException_

//#include "exceptions/HHEnergyRangeException.h"
#include "HHEnergyRangeException.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace HHKinFit2{
  class HHInvMConstraintException: public HHEnergyRangeException
{
public:
 HHInvMConstraintException(const std::string message)
   : HHEnergyRangeException(message) {};
};
}
#endif
