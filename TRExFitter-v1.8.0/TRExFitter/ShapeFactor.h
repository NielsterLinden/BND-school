#ifndef SHAPEFACTOR_H
#define SHAPEFACTOR_H

/// c++ includes
#include <string>
#include <vector>

class ShapeFactor{
public:
    explicit ShapeFactor();
    
    explicit ShapeFactor(const std::string& name, double nominal=1., double min=0., double max=10., bool isConst=false);
    
    ~ShapeFactor() = default;

    ShapeFactor(const ShapeFactor& s) = delete;
    ShapeFactor(ShapeFactor&& s) = delete;
    ShapeFactor& operator=(const ShapeFactor& s) = delete;
    ShapeFactor& operator=(ShapeFactor&& s) = delete;

    void Print() const;

    std::string fName;
    std::string fNuisanceParameter;
    std::string fTitle;
    std::string fCategory;

    double fNominal;
    double fMin;
    double fMax;
    bool fConst;

    std::vector<std::string> fRegions;
    std::vector<std::string> fSamples;
    std::vector<std::string> fExclude;
    std::vector<std::pair<std::string, std::string> > fExpression;

    int fNbins;
};

#endif
