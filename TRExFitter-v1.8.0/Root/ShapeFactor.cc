// Class include
#include "TRExFitter/ShapeFactor.h"

// Framework includes
#include "TRExFitter/Logger.h"

// -------------------------------------------------------------------------------------------------
// class ShapeFactor

//__________________________________________________________________________________
//
ShapeFactor::ShapeFactor() :
    fName(""),
    fNuisanceParameter(""),
    fTitle(""),
    fCategory(""),
    fNominal(0),
    fMin(0),
    fMax(0),
    fConst(false),
    fNbins(0) {
}

//__________________________________________________________________________________
//
ShapeFactor::ShapeFactor(const std::string& name, double nominal, double min, double max, bool isConst) :
    fName(name),
    fNuisanceParameter(name),
    fTitle(name),
    fCategory(""),
    fNominal(nominal),
    fMin(min),
    fMax(max),
    fConst(isConst),
    fNbins(0) {
}

//__________________________________________________________________________________
//
void ShapeFactor::Print() const{
    if (!fConst) LOG(DEBUG) << "        ShapeFactor: " << fName << "\t" << fNominal << ", " << fMin << ", " << fMax << "\n";
    else LOG(DEBUG) << "        ShapeFactor: " << fName << "\t" << fNominal << ", " << fMin << ", " << fMax << "(CONSTANT)\n";
}
