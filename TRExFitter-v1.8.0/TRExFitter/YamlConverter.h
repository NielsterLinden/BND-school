#ifndef YAMLCOVERTER_H_
#define YAMLCOVERTER_H_

#include "TRExFitter/LikelihoodScanManager.h"

#include "TH2.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

#include <fstream>
#include <string>
#include <vector>

namespace YAML {
    class Emitter;
}

class TGraphAsymmErrors;

class YamlConverter {

public:
    struct RankingContainer {
        std::string name;
        double nphat;
        double nperrhi;
        double nperrlo;
        double poihi;
        double poilo;
        double poiprehi;
        double poiprelo;
    };

    struct ImpactContainer {
        std::string name;
        double error;
        double errorHi;
        double errorLo;
    };

    struct TableContainer {
        std::vector<std::string> regionNames;
        std::vector<std::string> sampleNames;
        std::vector<std::vector<double> > mcYields;
        std::vector<std::vector<double> > mcErrors;
        std::vector<std::string> dataNames;
        std::vector<std::vector<double> > dataYields;
    };

    struct PlotContainer {
        std::vector<std::string> samples;
        std::vector<std::vector<double> > signalYields;
        std::vector<std::vector<double> > backgroundYields;
        TGraphAsymmErrors* errors;
        std::vector<double> data;
        std::vector<int> blindedBins;
        std::string xAxis;
        std::string yAxis;
        std::string region;
    };

    struct SubmissionContainer {
        bool useTables;
        std::vector<std::string> unfoldingNames;
        std::string folder;
        std::vector<std::string> regionNames;
        std::vector<std::string> useLikelihoodScan;
        std::vector<std::string> use2DLikelihoodScan;
    };

    struct Systematic {
        std::string name;
        std::vector<double> impact;
        double pull;
        double constraint;
    };

    struct Measurement {
        std::vector<std::string> pois;
        std::vector<double> values;
        std::vector<Systematic> systematics;
        TMatrixDSym npCorrelation;
        TMatrixDSym statCov;
    };

    explicit YamlConverter();
    ~YamlConverter() = default;

    YamlConverter(const YamlConverter& c) = delete;
    YamlConverter& operator=(const YamlConverter& c) = delete;
    YamlConverter(const YamlConverter&& c) = delete;
    YamlConverter& operator=(const YamlConverter&& c) = delete;

    void SetLumi(const std::string& lumi){ m_lumi = lumi; }
    void SetCME(const std::string& cme){ m_cme = cme; }

    void WriteRanking(const std::vector<RankingContainer>& ranking,
                      const std::string& path,
                      const bool shiftedGOs) const;

    void WriteRankingHEPData(const std::vector<RankingContainer>& ranking,
                             const std::string& folder,
                             const std::string& suffix,
                             const bool shiftedGOs) const;

    void WriteImpact(const std::vector<ImpactContainer>& impact,
                     const std::string& path) const;

    void WriteImpactHEPData(const std::vector<ImpactContainer>& impact,
                            const std::string& folder,
                            const std::string& suffix) const;

    void WriteCorrelation(const std::vector<std::string>& np,
                          const TMatrixD& corr,
                          const std::string& path,
                          const bool isCovariance,
                          const std::string& name="CorrelationMatrix") const;

    void WriteCorrelation(const std::vector<std::string>& np,
                          const std::vector<std::vector<double> >& corr,
                          const std::string& path,
                          const bool isCovariance,
                          const std::string& name="CorrelationMatrix") const;

    void WriteCorrelationHEPData(const std::vector<std::string>& np,
                                 const TMatrixD& corr,
                                 const std::string& path,
                                 const bool isCovariance,
                                 const std::string& prefix="Correlation") const;

    void WriteMigrationResponseHEPData(const TH2* matrix,
                                       const std::string& path,
                                       const bool isResponsei,
                                       const std::string& name,
                                       const std::string& region,
                                       const bool trueHorizontal) const;

    void WriteCorrelationHEPData(const std::vector<std::string>& np,
                                 const std::vector<std::vector<double> >& corr,
                                 const std::string& path,
                                 const bool isCovariance,
                                 const std::string& prefix="Correlation") const;


    void WriteTables(const TableContainer& container,
                     const std::string& directory,
                     const bool isPostFit) const;

    void WriteTablesHEPData(const TableContainer& container,
                            const std::string& directory,
                            const bool isPostFit) const;

    void WriteUnfolding(const TGraphAsymmErrors* const graph,
                        const std::string& directory,
                        const std::string& name) const;

    void WriteUnfoldingHEPData(const TGraphAsymmErrors* const graph,
                               const std::string& xAxis,
                               const std::string& directory,
                               const std::string& name,
                               const std::vector<std::unique_ptr<TH1> >& mc,
                               const std::vector<std::string>& legendNames) const;

    void WritePlot(const PlotContainer& container,
                   const std::string& directory,
                   const std::string& subdirectory,
                   const bool isPostFit) const;

    void WritePlotUhepp(const PlotContainer& container,
                   const std::string& directory,
                   const bool isPostFit) const;

    void WritePlotHEPData(const PlotContainer& container,
                          const std::string& directory,
                          const bool isPostFit) const;

    void WriteLikelihoodScan(const std::pair<std::vector<double>, std::vector<double> >& result,
                             const std::string& path) const;

    void WriteLikelihoodScan2D(const LikelihoodScanManager::Result2D& result,
                               const std::string& path) const;

    void WriteLikelihoodScanHEPData(const std::pair<std::vector<double>, std::vector<double> >& result,
                                    const std::string& folder,
                                    const std::string& suffix) const;

    void WriteLikelihoodScan2DHEPData(const LikelihoodScanManager::Result2D& result,
                                      const std::string& folder,
                                      const std::string& suffix) const;

    void WriteHEPDataSubmission(const SubmissionContainer& container,
                                const std::vector<std::string>& pois) const;

    void AppendSubmissionWithCovariances(const std::vector<std::string>& params,
                                         const std::string& path) const;

    void WriteCombiner(const std::string& name, const Measurement& measurement, const std::string& path) const;

    /**
     * @brief Add $ delimiters to the string if it contains '\', `{`, or `}`
     *
     * @param input
     * @return std::string
     */
    static std::string AddLatexDelimiters(const std::string& input);

private:
    std::string m_lumi;
    std::string m_cme;

    void AddQualifiers(YAML::Emitter& out) const;

    void AddValueErrors(YAML::Emitter& out,
                       const double mean,
                       const double up,
                       const double down) const;

    void Write(const YAML::Emitter& out, const std::string& type, const std::string& path) const;

    bool TableContainerIsOK(const TableContainer& container) const;

    bool PlotContainerIsOK(const PlotContainer& container) const;

    void AddCorrelation(std::ofstream& file) const;

    void AddRanking(std::ofstream& file, const std::vector<std::string>& pois) const;

    void AddImpact(std::ofstream& file, const std::vector<std::string>& pois) const;

    void AddPlots(std::ofstream& file, const SubmissionContainer& container) const;

    void AddTables(std::ofstream& file) const;

    void AddLikelihoodScan(std::ofstream& file, const std::string& scan) const;

    void Add2DLikelihoodScan(std::ofstream& file, const std::string& scan) const;

    void Add(std::ofstream& o, const std::string& file, const std::string& text) const;
};

#endif
