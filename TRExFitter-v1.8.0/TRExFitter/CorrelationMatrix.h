#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H

/// c++ includes
#include <string>
#include <vector>
#include <map>

class CorrelationMatrix {

public:
    explicit CorrelationMatrix();
    ~CorrelationMatrix() = default;

    CorrelationMatrix(const CorrelationMatrix& m) = delete;
    CorrelationMatrix(CorrelationMatrix&& m) = delete;
    CorrelationMatrix& operator=(const CorrelationMatrix& m) = delete;
    CorrelationMatrix& operator=(CorrelationMatrix&& m) = delete;

    //
    // Functions
    //
    void AddNuisPar(const std::string& p);
    void Resize(const int size);
    void SetCorrelation(const std::string& p0, const std::string& p1,double corr);
    void SetPlotLabel(const std::string& l){fPlotLabel = l;}
    double GetCorrelation(const std::string& p0, const std::string& p1) const;
    bool NPIsPresent(const std::string& np) const;

    /**
     * @brief Draw the correlation matrix
     *
     * @param path
     * @param useHEPDataFormat
     * @param corrMin
     * @param POIs
     */
    void Draw(const std::string& path, const bool useHEPDataFormat, const double corrMin = -1., const std::vector<std::string>& POIs={} );

    //
    // Data members
    //
    std::vector<std::string> fNuisParToHide;
    std::vector<std::string> fNuisParList;
    std::vector<std::string> fEFTParList;
    std::vector<std::vector<double> > fMatrix;
    std::string fOutFolder;
    std::string fPlotLabel;

    private:

    std::vector<std::string> fNuisParNames;
};

#endif
