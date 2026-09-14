#ifndef HISTOREADER_H_
#define HISTOREADER_H_

#include <memory>
#include <string>
#include <vector>

class Sample;
class SampleHist;
class Systematic;
class TH1;
class TRExFit;

#include <set>

/**
  * \class HistoReader
  * \brief Class that reads histograms as an input
  */

class HistoReader {

    public:
        /**
          * The constructor
          * @param fitter A pointer to TRExFitter class
          */
        explicit HistoReader(TRExFit* fitter);

        /**
          * Destructor that does nothing
          */
        ~HistoReader();

        /**
          *
          * Deleted constructors and assignment operators
          */
        HistoReader() = delete;
        HistoReader(const HistoReader& h) = delete;
        HistoReader(HistoReader&& h) = delete;
        HistoReader& operator=(const HistoReader& h) = delete;
        HistoReader& operator=(HistoReader&& h) = delete;

        /**
          * method that reads the histograms from user defined inputs
          */
        void ReadHistograms(const std::vector<std::string>& binningsToCreate);

        /**
          * method that reads the histograms produced by TREx
          */
        void ReadTRExProducedHistograms(const bool isUnfolding, const bool& checkHistogram);

    private:
        /// A pointer to the TRExFit class
        TRExFit* fFitter;

        /**
          *
          * A helper function to read singleHistogram from a file
          * @param vector of histo paths
          * @param vector of histo paths for subtraction
          * @param pointer to systematic
          * @param index of channel
          * @param index of sample
          * @param flag if the systematic is up or down
          * @param flag if we are reading MC or Data
          * @return the read histogram
          */
         std::unique_ptr<TH1> ReadSingleHistogram(const std::vector<std::string>& fullPaths,
                                                  const std::vector<std::string>& subtractPaths,
                                                  Systematic* syst,
                                                  int i_ch,
                                                  int i_smp,
                                                  bool isUp,
                                                  bool isMC);

        /**
          * A helper function to read one region for data samples
          * @param index of the region
          * @param is data
          */
        void ReadOneRegion(const int i_ch, const bool is_data, const bool recreateBinning);

        /**
          * A helper function to read shape and norm factors for samples
          * @param channel index
          * @param pointer to a sample
          */
        void ReadNormShape(SampleHist* sh,
                           const int i_ch,
                           const Sample* smp);

        /**
          * A helper function to read systeamtics
          * @param channel index
          * @param pointer to a systematic
          * @return true if systematics is overall
          */
        bool SetSystematics(const int i_ch,
                            const Sample* ismp,
                            std::shared_ptr<Systematic> syst);

        /**
         * @brief Get the Systenatic histograms
         *
         * @param i_ch
         * @param i_smp
         * @param syst
         * @param is_data
         * @param is_up
         * @return std::unique_ptr<TH1>
         */
        std::unique_ptr<TH1> GetSystHisto(const int i_ch,
                                          const std::size_t i_smp,
                                          Systematic* syst,
                                          bool is_data,
                                          bool is_up);


        /**
         * @brief Propagate the reference sample to the variation
         *
         * @param h histogram to be modified
         * @param syst systematic
         * @param i_ch channel index
         * @param i_smp sample index
         */
        void PropagateReference(TH1* h, Systematic* syst, const int i_ch, const std::size_t i_smp) const;
};

#endif
