#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

/// c++ includes
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>

// Functions
std::string Fix(const std::string& s);
std::string First(const std::string& s);
std::string Second(const std::string& s, const std::string& file_name, const unsigned int line_number);
std::string ReadValueFromConfig(const std::string& fileName,const std::string& option);

struct ConfigLine {
    std::string text;
    std::string file_name;
    unsigned int line_number;
};

// classes
class Config {
public:
    explicit Config();
    ~Config() = default;
    Config(const Config& c) = default;
    Config(Config&& c) = default;
    Config& operator=(const Config& c) = default;
    Config& operator=(Config&& c) = default;

    std::string fName;
    std::string fValue;
    std::string fFileName;
    unsigned int fLineNumber;
};

class ConfigSet {
public:
    explicit ConfigSet();
    ~ConfigSet() = default;
    ConfigSet(const ConfigSet& c) = delete;
    ConfigSet(ConfigSet&& c) = delete;
    ConfigSet& operator=(const ConfigSet& c) = delete;
    ConfigSet& operator=(ConfigSet&& c) = delete;

    void SetConfig(const std::string& name,const std::string& value, const std::string& file_name, const unsigned int line_number);
    std::string Get(const std::string& name) const;
    std::string operator[](const std::string& name) const;
    const Config& GetConfig(int i) const;
    const std::string& GetConfigName(int i) const;
    const std::string& GetConfigValue(int i) const;
    int GetN() const;
    int size() const;
    void Set(const std::string& name,const std::string& value);
    const std::string& GetName() const;
    const std::string& GetValue() const;

    int fN;
    std::string fName;
    std::string fValue;
    std::vector<Config> fConfig;
};


class ConfigParser {
public:
    explicit ConfigParser();
    ~ConfigParser() = default;
    ConfigParser(const ConfigParser& c) = delete;
    ConfigParser(ConfigParser&& c) = delete;
    ConfigParser& operator=(const ConfigParser& c) = delete;
    ConfigParser& operator=(ConfigParser&& c) = delete;

    std::unordered_map<std::string, std::vector< std::unique_ptr<ConfigSet>>> fConfSets;
    void ReadFile(const std::string& fileName);
    const std::vector< std::unique_ptr<ConfigSet> >* GetConfigSets(const std::string& name);
    const ConfigSet *GetConfigSet(const std::string& name, unsigned int i=0);
    int CheckSyntax(ConfigParser *refConfigParser);

    /**
      * Method that checks provided setting with reference config file
      * @param ConfigSet Pointer to actual configuration setting type
      * @param ConfigSet Pointer to referece config parser
      * @param string Name of the setting set
      * @param string Name of the setting
      * @return int Status code
      */
    int SettingIsValid(const ConfigSet *cs, ConfigParser *refConfigParser, const std::string &setting_set, const Config &setting) const;

    /**
      * Helper method that checks single setting with reference file
      * @param ConfigSet Pointer to actual configuration setting type
      * @param ConfigSet Pointer to referece config parser
      * @param string Name of the setting set
      * @param string Name of the setting
      * @return int Status code
      */
    int CheckSingleSetting(const ConfigSet *cs, const ConfigSet *cs_ref, const std::string &setting_set, const Config &setting) const;

    /**
      * Helper method that checks validity of provided parameters for given setting with reference file
      * @param string Setting in users config
      * @param string vector List of possible settings
      * @param string Name of the setting set
      * @param string Name of the setting
      * @return int Status code
      */
    int CheckParameters(std::string current, const std::vector<std::string> &possible_settings, const std::string& setting_set, const Config &setting) const;

    /**
      * @param string Name of the setting set
      * @param string Name of the setting
      * @param string Setting value in users config
      * @param string Possible setting from reference file
      * @param char Delimiter whn multiple parameters are provided
      * @return bool Is valid setting
      */
    bool SettingMultipleParamIsOK(const std::string &setting_set, const Config& setting, const std::string& current, const std::string& possible, const char delimiter = ',') const;

private:
    /**
     * Read the config file and replace INCLUDE lines by the content of the included config. The lines of the read config will be pushed back to lines_of_config vector
     * @param string address of the config
     * @param *vector<string> pointer to vector, where lines of the configs will be added
     * @param *vector<string> pointer to the vector with already included configs (to avoid circular include)
    */
    void ExpandConfigs(std::string config_address, std::vector<ConfigLine> *lines_of_config, std::vector<std::string> *included_configs, std::map<std::string,std::string>& fReplacement, bool doIncludes=true);

    std::map<std::string,std::string> ParseReplacements(const std::vector<ConfigLine>& lines_of_config, const std::string& fileName);
    void ReplaceFromReplacementFile(std::string* line, const std::map<std::string,std::string>& fReplacement);

    static void ResolvePath(std::string *config_address);
};

#endif
