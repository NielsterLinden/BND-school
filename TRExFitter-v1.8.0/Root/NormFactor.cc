// Class include
#include "TRExFitter/NormFactor.h"

// Framework includes
#include "TRExFitter/Logger.h"

// -------------------------------------------------------------------------------------------------
// class NormFactor

//__________________________________________________________________________________
//
NormFactor::NormFactor() :
    fName(""),
    fNuisanceParameter(""),
    fTitle(""),
    fCategory(""),
    fSubCategory(""),
    fConst(false),
    fTau(0),
    fCorrelateUnfolding(""),
    fNominal(0),
    fMin(0),
    fMax(0)
{
}

//__________________________________________________________________________________
//
NormFactor::NormFactor(const std::string& name, double nominal, double min, double max, bool isConst, const std::string& subCategory) :
    fName(name),
    fNuisanceParameter(name),
    fTitle(name),
    fCategory(""),
    fSubCategory(subCategory),
    fConst(isConst),
    fTau(0),
    fCorrelateUnfolding(""),
    fNominal(nominal),
    fMin(min),
    fMax(max)
{
    if (nominal < min) {
        LOG(ERROR) << "NormFactor " << fName << " has Nominal value smaller than Min\n";
        exit(EXIT_FAILURE);
    }
    if (nominal > max) {
        LOG(ERROR) << "NormFactor " << fName << " has Nominal value larger than Max\n";
        exit(EXIT_FAILURE);
    }
}

//__________________________________________________________________________________
//
void NormFactor::Print() const{
    if (fConst) LOG(INFO) << fName << "\t" << fNominal << ", " << fMin << ", " << fMax << "  (CONSTANT)\n";
    else LOG(INFO) << fName << "\t" << fNominal << ", " << fMin << ", " << fMax << "\n";
}

//__________________________________________________________________________________
//
void NormFactor::SetNominalMinMax(const double nominal, const double min, const double max) {
    if (nominal < min) {
        LOG(ERROR) << "NormFactor " << fName << " has Nominal value smaller than Min\n";
        exit(EXIT_FAILURE);
    }
    if (nominal > max) {
        LOG(ERROR) << "NormFactor " << fName << " has Nominal value larger than Max\n";
        exit(EXIT_FAILURE);
    }

    fNominal = nominal;
    fMin     = min;
    fMax     = max;
}