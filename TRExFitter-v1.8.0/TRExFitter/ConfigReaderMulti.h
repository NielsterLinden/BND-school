#ifndef CONFIGREADERMULTI_H
#define CONFIGREADERMULTI_H

/// c++ includes
#include <string>

/// Forward class declaration
class ConfigParser;
class MultiFit;

/**
 * \class ConfigReaderMulti
 * \brief Class that reads config parameters
 * \author Tomas Dado
 */

class ConfigReaderMulti {
    public:

        /**
          * The default constructor.
          * @param MultiFIt A pointer to MultiFit class
          */
        explicit ConfigReaderMulti(MultiFit* multiFitter);

        /**
          * The default destructor
          */
        ~ConfigReaderMulti() = default;
        
        /**
          * Deleted constructors
          */ 
        ConfigReaderMulti() = delete;
        ConfigReaderMulti(const ConfigReaderMulti& c) = delete;
        ConfigReaderMulti(ConfigReaderMulti&& c) = delete;
        ConfigReaderMulti& operator=(const ConfigReaderMulti& c) = delete;
        ConfigReaderMulti& operator=(ConfigReaderMulti&& c) = delete;

        /**
          * Reads the config and passes parameters to TtHFit
          * @param string Config path
          * @param string Running options
          * @param string Additional options
          * @return int status code
          */
        int ReadFullConfig(const std::string& fileName, const std::string& opt, const std::string& option);

        /**
         * @brief Adjust config options after the config has been read 
         * 
         * @return int status code 
         */
        int PostConfig();

    private:

        /**
          * Helper function to read Command line settings
          * @param string config options
          * @return int status code
          */
        int ReadCommandLineOptions(const std::string &option);

        /**
          * Helper function to read JOB settings
          * @return int status code
          */
        int ReadJobOptions();

        /**
          * Helper function to read Limit settings
          * @return int status code
          */
        int ReadLimitOptions();

        /**
          * Helper function to read Significance settings
          * @return int status code
          */
        int ReadSignificanceOptions();

        /**
          * Helper function to read Fit settings
          * @param string Running options
          * @param string Option flag
          * @return int status code
          */
        int ReadFitOptions(const std::string& opt, const std::string& option);

        /**
          * Pointer to MultiFit class, set during initialization
          */
        MultiFit *fMultiFitter;

        /**
          * Instance of ConfigParser used to parse the text
          */
        ConfigParser fParser;

        /**
          * String Used for global settings
          */
        std::string fGlobalSuffix;

	/**
	  * flag for whether LHScan is performed
	  */
	bool fRunLHscan;



        /**
          * string To keept track of command like setting for LH scan
          */
	std::string fLHscanVariable;

	/**
          *  Likelihood single scan point to process
          */
        int fLHscanStep;
	int fLHscanStepY;


};

#endif
