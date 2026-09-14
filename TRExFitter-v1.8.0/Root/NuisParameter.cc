// Class includes
#include "TRExFitter/NuisParameter.h"

//__________________________________________________________________________________
//
NuisParameter::NuisParameter(const std::string& name) :
    fName(name),
    fTitle(name),
    fCategory(""),
    fStartValue(0),
    fFitValue(0),
    fPostFitUp(1),
    fPostFitDown(-1),
    fConstrainType(0), // 0 -> Gaussian
    fInterpCode(4) { // 4 -> Pol interp. + Exp extrap.
}
