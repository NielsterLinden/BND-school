#ifndef SIGNIFICANCETOYS_H
#define SIGNIFICANCETOYS_H

#include <string>

class RooAbsData;
namespace RooStats{
    class ModelConfig;
    class HypoTestResult;
}


/**
 * @brief Class for calculating significance with toys
 *
 */
class SignificanceToys {
public:

    /**
     * @brief Construct a new Significance Toys object
     *
     */
    explicit SignificanceToys();

    /**
     * @brief Destroy the Significance Toys object
     *
     */
    ~SignificanceToys() = default;

    /**
     * @brief Set the Ntoys object
     *
     * @param splusb number of S+B toys
     * @param bonly number of B toys
     */
    inline void SetNtoys(const int splusb, const int bonly) {
        fNtoysSplusB = splusb;
        fNtoysB = bonly;
    }

    /**
     * @brief Set the seeds for the toys
     *
     * @param seed
     */
    inline void SetToysSeed(const int seed) {fToysSeed = seed;}

    /**
     * @brief Set if the plots shouldbe produced
     *
     * @param flag
     */
    inline void SetPlot(const bool flag) {fPlot = flag;}

    /**
     * @brief Set the Output Path
     *
     * @param path
     */
    inline void SetOutputPath(const std::string& path) {fOutputPath = path;}

    /**
     * @brief Tuns the toys to estimate significance
     *
     * @param data Data
     * @param mcSplusb S+B model config
     * @param mcB B-only model config
     */
    void RunToys(RooAbsData* data,
                 RooStats::ModelConfig* mcSplusb,
                 RooStats::ModelConfig* mcB) const;

    /**
     * @brief Output as TTree
     *
     * @param result
     */
    void RootOutput(RooStats::HypoTestResult* result) const;

    /**
     * @brief Set the Significance Suffix
     *
     * @param suffix
     */
    void SetSignificanceSuffix(const std::string& suffix) {fSignificanceSuffix = suffix;}

private:
    /**
     * @brief number of S+B toys
     *
     */
    int fNtoysSplusB;

    /**
     * @brief number of B toys
     *
     */
    int fNtoysB;

    /**
     * @brief The seed for toys
     *
     */
    int fToysSeed;

    /**
     * @brief Plot?
     *
     */
    bool fPlot;

    /**
     * @brief Output path
     *
     */
    std::string fOutputPath;

    /**
     * @brief Set suffix for the root file
     *
     */
    std::string fSignificanceSuffix;
};
#endif
