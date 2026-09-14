// Class include
#include "TRExFitter/ConfigParser.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"

// c++ includes
#include <algorithm>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>


// using namespace std;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// -- Functions --
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------


//__________________________________________________________________________________
// Removes leading and trailing white spaces, removes everything after "%", removes '"'
std::string Fix(const std::string& s){
    if(s=="") return "";
    std::string ss = s.substr( 0, s.find_first_of('%') );
    replace( ss.begin(), ss.end(), '"', ' ');
    if (ss.find_first_not_of(' ')>=std::string::npos){
        ss = "";
    }
    else if (ss.find_first_not_of(' ')>0){
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' '));
    }
    else{
        ss=ss.substr(ss.find_first_not_of(' '),ss.find_last_not_of(' ')+1);
    }
    return ss;
}

//__________________________________________________________________________________
//
std::string First(const std::string& s){
    std::string first;
    first = s.substr( 0, s.find_first_of(':') );
    return Fix(first);
}

//__________________________________________________________________________________
//
std::string Second(const std::string& s, const std::string& file_name, const unsigned int line_number){
    std::string second;
    second = s.substr( s.find_first_of(':')+1,std::string::npos );
    second = Common::RemoveComments(second);
    if(second==""){
        LOG(ERROR) << "Error in line " << line_number << " of " << file_name << ": No value set for parameter " << First(s) << " in the config.\n";
        exit(EXIT_FAILURE);
    }
    return second;
}

//_______________________________________________________________________________________
//
void ReplaceStringInPlace(std::string& subject, const std::string& search,
                                                const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
}

//_______________________________________________________________________________________
// used to pre-read a config file to check single option values
std::string ReadValueFromConfig(const std::string& fileName,const std::string& option){
    std::string value = "";
    std::ifstream file(fileName.c_str());
    if(!file.is_open()) return value;
    std::string str;
    while(true){
        file >> str;
        if(!file.good()) break;
        if(str==(option+":")){
            file >> value;
            break;
        }
    }
    file.close();
    file.clear();
    return value;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// -- Classes --
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
// Config
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
Config::Config() :
    fName(""),
    fValue(""),
    fFileName(""),
    fLineNumber(0) {
}

//----------------------------------------------------------------------------------
// ConfigSet
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
ConfigSet::ConfigSet() :
    fN(0),
    fName(""),
    fValue("") {
}

//__________________________________________________________________________________
//
void ConfigSet::SetConfig(const std::string& name,const std::string& value, const std::string& file_name, const unsigned int line_number){
    // check if is name is there
    bool isThere = false;
    int index = -1;
    for(int i=0;i<fN;i++){
        if(fConfig[i].fName == name){
            isThere = true;
            index = i;
        }
    }
    if(!isThere){
        fConfig.emplace_back();
        fConfig.back().fName = name;
        fConfig.back().fValue = value;
        fConfig.back().fFileName = file_name;
        fConfig.back().fLineNumber = line_number;
        fN++;
    }
    else{
        fConfig[index].fName = name;
        fConfig[index].fValue = value;
        fConfig[index].fFileName = file_name;
        fConfig[index].fLineNumber = line_number;
    }
}

//__________________________________________________________________________________
//
std::string ConfigSet::Get(const std::string& name) const{
    for(int i=0;i<fN;i++){
        if(fConfig[i].fName == name){
            return fConfig[i].fValue;
        }
    }
    return "";
}

//__________________________________________________________________________________
//
std::string ConfigSet::operator[](const std::string& name) const{
    return Get(name);
}

//__________________________________________________________________________________
//
const Config& ConfigSet::GetConfig(int i) const{
    return fConfig[i];
}

//__________________________________________________________________________________
//
const std::string& ConfigSet::GetConfigName(int i) const{
    return fConfig[i].fName;
}

//__________________________________________________________________________________
//
const std::string& ConfigSet::GetConfigValue(int i) const{
    return fConfig[i].fValue;
}

//__________________________________________________________________________________
//
int ConfigSet::GetN() const{
    return fN;
}

//__________________________________________________________________________________
//
int ConfigSet::size() const{
    return fN;
}

//__________________________________________________________________________________
//
void ConfigSet::Set(const std::string& name,const std::string& value){
    fName = name;
    fValue = value;
}

//__________________________________________________________________________________
//
const std::string& ConfigSet::GetName() const{
    return fName;
}

//__________________________________________________________________________________
//
const std::string& ConfigSet::GetValue() const{
    return fValue;
}

//----------------------------------------------------------------------------------
// ConfigParser
//----------------------------------------------------------------------------------

//__________________________________________________________________________________
//
ConfigParser::ConfigParser() {}

void ConfigParser::ResolvePath(std::string *config_address) {
    std::vector<std::string> elements = Common::SplitString(*config_address, "/");
    if (elements.size() == 0)   {
        return;
    }
    const bool absolute_path = elements[0] == "";

    for (unsigned int i = 0; i < elements.size(); )  {
        if (elements[i] == "" || elements[i] == ".")    {
            elements.erase(elements.begin() + i);
            if (i > 0)  {
                i--;
            }
            continue;
        }

        if (elements[i] != "..")    {
            if (i+1 < elements.size())  {
                if (elements[i+1] == "..")  {
                    elements.erase(elements.begin() + i);
                    elements.erase(elements.begin() + i);
                    if (i > 0)  {
                        i--;
                    }
                    continue;
                }
            }
        }
        i++;
    }
    *config_address = Common::join_strings(elements, "/");
    if (absolute_path)  {
        *config_address = "/" + *config_address;
    }
}

void ConfigParser::ReplaceFromReplacementFile(std::string* line, const std::map<std::string,std::string>& fReplacement){
    std::map<std::string,std::string>::const_iterator itr=fReplacement.begin();
    for ( ;itr!=fReplacement.end();++itr) {
        std::string oldV=itr->first;
        std::string newV=itr->second;
        ReplaceStringInPlace(*line, oldV, newV);
    }
}

std::map<std::string,std::string> ConfigParser::ParseReplacements(const std::vector<ConfigLine>& lines_of_config, const std::string& fileName){
    std::string str;
    std::map<std::string,std::string> fReplacement;
    std::string replacementFileName="";
    // loop over the file once to automatically find the replacement fileName

    for (const ConfigLine &line : lines_of_config) {
        str = line.text;
        replace( str.begin(), str.end(), '\n', ' ');
        replace( str.begin(), str.end(), '\r', ' ');
        replace( str.begin(), str.end(), '\t', ' ');
        if (str.find_first_not_of(' ') != std::string::npos) {
            if (str[str.find_first_not_of(' ')] == '%') continue;
            if (str[str.find_first_not_of(' ')] == '#') continue;
        }
        if ( str.find("ReplacementFile")==std::string::npos ) continue;
        replacementFileName=Second(str, line.file_name, line.line_number);
        break;
    }
    if (replacementFileName!="" && (fileName.find("jobSchema.config") == std::string::npos)) {
        LOG(INFO) << "Opening replacement file: " << replacementFileName << " to fill the map\n";
        replacementFileName = Common::RemoveQuotes(replacementFileName);
        std::ifstream fileR(replacementFileName.c_str());
        if(!fileR.is_open()){
            LOG(ERROR) << "The replacement file: " << replacementFileName << " cannot be opened!\n";
            exit(-1);
        }
        unsigned int line_number = 0;
        while (getline(fileR, str)){
            ++line_number;
            replace( str.begin(), str.end(), '\n', ' ');
            replace( str.begin(), str.end(), '\r', ' ');
            replace( str.begin(), str.end(), '\t', ' ');
            if ( str.find("%")!=std::string::npos || str.find("#")!=std::string::npos) continue;
            if (str.find_first_not_of(' ') == std::string::npos)    continue;
            std::string key = First(str);
            std::string value = Second(str, replacementFileName, line_number);
            LOG(INFO) << "Putting in the map: [" << key << "]=" << value << "\n";
            fReplacement[key]=value;
        }
        fileR.close();
    }
    return fReplacement;
}

void ConfigParser::ExpandConfigs(std::string config_address, std::vector<ConfigLine> *lines_of_config,
                                std::vector<std::string> *included_configs, std::map<std::string,std::string>& fReplacement,
                                bool doIncludes)    {

    ResolvePath(&config_address);
    if (std::count(included_configs->begin(), included_configs->end(), config_address))    {
        const std::string configs = Common::join_strings(*included_configs, " -> ") + " -> " + config_address;
        LOG(ERROR) << "Detected circular include of configs: " << configs << "\n";
        exit(EXIT_FAILURE);
    }
    included_configs->push_back(config_address);

    std::ifstream input_file (config_address);
    if (input_file.is_open())    {
        std::string line;
        unsigned int line_number = 0;
        while ( getline (input_file,line) )        {
            ++line_number;
            Common::StripString(&line, "");

            std::string stripped_line = line;
            Common::StripString(&stripped_line);
            if (Common::StartsWith(stripped_line, "INCLUDE:")) {
                if ((line.find("XXX") != std::string::npos) && (!fReplacement.empty())) {
                    ReplaceFromReplacementFile(&line, fReplacement);
                }
                line = Common::RemoveQuotes(line);
                if (!doIncludes)   continue;
                std::vector<std::string> current_config_address_elements = Common::SplitAndStripString(config_address, "/");
                current_config_address_elements.pop_back();

                std::string new_config = line.substr(8);
                Common::StripString(&new_config);

                std::vector<std::string> new_config_elements = Common::SplitString(new_config, "/");
                if (new_config_elements.size() == 0)   {
                    LOG(ERROR) << "INCLUDE flag must be followed by a path\n";
                    exit(EXIT_FAILURE);
                }
                if ((current_config_address_elements.size() == 0)) {current_config_address_elements = {"."};}
                const bool absolute_path = new_config_elements[0] == "";


                if (!absolute_path)  {
                     new_config =  (Common::join_strings(current_config_address_elements, "/")) + "/" + new_config;
                }
                ExpandConfigs(new_config, lines_of_config, included_configs, fReplacement);
            }
            else {
                lines_of_config->push_back({line, config_address, line_number});
            }
        }
    }
    else    {
        LOG(ERROR) << "Unable to open config file \"" << config_address << "\"\n";
        exit(EXIT_FAILURE);
    }
    included_configs->pop_back();
}

//__________________________________________________________________________________
//

void ConfigParser::ReadFile(const std::string& fileName){
    if (!std::filesystem::is_regular_file(fileName)) {
        LOG(ERROR) << "Path: " << fileName << " is not a regular file path\n";
        exit(EXIT_FAILURE);
    }
    LOG(INFO) << "Reading config file: " << fileName << "\n";
    std::vector<ConfigLine> tmp_lines_of_config, lines_of_config;
    std::vector<std::string> included_configs;

    // If a replacement file is given in main config, not via an include, then parse it before parsing INCLUDES
    // this allows using replaced strings in INCLUDE: lines within the config
    std::map<std::string,std::string> fReplacement;
    // Expand the main config only, without expanding included configs
    // At this point, fReplacement is empty, no replacements are made
    ExpandConfigs(fileName, &tmp_lines_of_config, &included_configs, fReplacement, false);
    // If a replacement file was found in main config, create replacement mapping from it
    // If a replacement file was notfound in main config, it may be in an included config
    fReplacement = ParseReplacements(tmp_lines_of_config, fileName);

    // Expand all included configs -- if fReplacement was defined from main config
    // all lines willl be replaced on the fly as included configs are parsed
    ExpandConfigs(fileName, &lines_of_config, &included_configs, fReplacement);

    // if fReplacement was not defined from main config, we can look for it again
    // after we have expanded all includes
    if (fReplacement.empty())   fReplacement = ParseReplacements(lines_of_config, fileName);
    // at this point, if fReplacement is empty, there is no replacment file within main and nested configs
    // if fReplacement is not empty, it will be used to replace XXX in the fully expanded config.

    bool reading = false;
    unsigned int n = 1;
    std::string config_set_name = "";
    for (auto& config_line : lines_of_config) {
        std::string& line = config_line.text;
        replace( line.begin(), line.end(), '\n', ' ');
        replace( line.begin(), line.end(), '\r', ' ');
        replace( line.begin(), line.end(), '\t', ' ');
        if (line.find_first_not_of(' ') != std::string::npos) {
            if(line[line.find_first_not_of(' ')]=='%') continue;
            if(line[line.find_first_not_of(' ')]=='#') continue;
        }
        //
        if (line.find("XXX")!=std::string::npos) {
            LOG(INFO) << "BEFORE replacement: " << line << "\n";
            ReplaceFromReplacementFile(&line, fReplacement);
            LOG(INFO) << "AFTER replacement: " << line << "\n";
        }
        //
        if(line.find_first_not_of(' ')==std::string::npos){
            reading = false;
            continue;
        }
        const std::vector<std::string> valVec = Common::Vectorize(Second(line, config_line.file_name, config_line.line_number),';',false);
        const std::string param_name = First(line);
        if(!reading){
            n = valVec.size();
            config_set_name = param_name;
            auto& confsets = fConfSets.try_emplace(config_set_name).first->second;
            for(const auto& val: valVec){
                auto cs = std::make_unique<ConfigSet>();
                cs->Set(config_set_name, Common::RemoveSpaces(val));
                confsets.push_back(std::move(cs));
            }
            reading = true;
        }
        else{
            std::string val;
            // calculate offset to fill the last n ConfigSet objects
            const unsigned int offset = fConfSets[config_set_name].size() - n;
            for(unsigned int k=0; k<n; k++){
                if (valVec.empty()) {
                    val = "";
                } else {
                    if(k>=valVec.size()) {
                        val = valVec.back();
                    } else {
                        val = valVec[k];
                    }
                }
                fConfSets[config_set_name][k+offset]->SetConfig(param_name, Common::RemoveSpaces(val), config_line.file_name, config_line.line_number);
            }
        }
    }
}

const std::vector< std::unique_ptr<ConfigSet> >* ConfigParser::GetConfigSets(const std::string& name) {
    auto it = fConfSets.find(name);
    if(it==fConfSets.end()) return nullptr;
    return &it->second;
}

//__________________________________________________________________________________
// Returns the i-th configSet with given name
const ConfigSet *ConfigParser::GetConfigSet(const std::string& name, unsigned int i){
    const auto* confsets = GetConfigSets(name);
    if(!confsets) return nullptr;
    if(i>=confsets->size()) return nullptr;

    return (*confsets)[i].get();
}

//__________________________________________________________________________________
//
int ConfigParser::CheckSyntax(ConfigParser *refConfigParser){
    int exitStatus = 0;
    // loop over all config sets
    for(const auto& [_, confsets]: fConfSets) {
        for(const auto& cs: confsets) {
            // loop over all settings in the config set
            for(const auto& setting: cs->fConfig) {
                exitStatus += SettingIsValid(cs.get(), refConfigParser, cs->fName, setting);
            }
        }
    }

    return exitStatus;
}

//_______________________________________________________________________________________
//
int ConfigParser::SettingIsValid(const ConfigSet *cs, ConfigParser *refConfigParser, const std::string &setting_set, const Config &setting) const{
    if (refConfigParser == nullptr){
        LOG(ERROR) << "Invalid pointer to the reference ConfigParser. Please check this!\n";
        exit(EXIT_FAILURE);
    }
    const ConfigSet *cs_ref = nullptr;
    // check if the setting type exists in the reference
    auto it = refConfigParser->fConfSets.find(setting_set);;
    if(it!=refConfigParser->fConfSets.end()) {
        cs_ref = it->second.front().get();
    }
    //
    // config set is not present in the reference
    if (cs_ref == nullptr){
        LOG(ERROR) << "Error in line " << setting.fLineNumber << " of " << setting.fFileName << ": Cannot find config set '" << setting_set << "' in reference config. Please check this!\n";
        return 1;
    }
    //
    // check the validity of the single setting
    return CheckSingleSetting(cs, cs_ref, setting_set, setting);
}

//_______________________________________________________________________________________
//
int ConfigParser::CheckSingleSetting(const ConfigSet *cs, const ConfigSet *cs_ref, const std::string &setting_set, const Config &setting) const {
    std::string param = cs->Get(setting.fName);
    std::string ref_param = cs_ref->Get(setting.fName);
    //
    // this option used to exist, write a message informing the user of the new syntax
    if (setting.fName == "TtresSmoothing"){
        LOG(ERROR) << "Error in line " << setting.fLineNumber << " of " << setting.fFileName << ": The TtresSmoothing option is deprecated, use SmoothingOption: TTBARRESONANCE instead.\n";
        return 1;
    }
    else if (ref_param == ""){
        LOG(ERROR) << "Error in line " << setting.fLineNumber << " of " << setting.fFileName << ": Cannot find setting '" << setting.fName << "' for setting set " << setting_set <<  " in reference config. Please check this!\n";
        return 1;
    }
    //
    // there is nothing to check if the reference setting is simply 'std::string'
    if (ref_param == "string") return 0;
    //
    // need to check the consistency of the provided settings
    // first check the number of provided parameters
    std::vector<std::string> current_settings = Common::Vectorize(param,',');
    std::vector<std::string> possible_settings = Common::Vectorize(ref_param,'/');
    std::vector<unsigned int> possible_sizes;
    //
    for (const std::string &ioption : possible_settings){
        possible_sizes.push_back( Common::Vectorize(ioption,',').size() );
    }
    unsigned int current_setting_size = current_settings.size();
    //
    // search the vector for the possible sizes
    auto it = std::find(possible_sizes.begin(), possible_sizes.end(), current_setting_size);
    if (it == possible_sizes.end()){
        LOG(ERROR) << "Error in line " << setting.fLineNumber << " of " << setting.fFileName << ": You provided " << current_setting_size <<" parameters, for setting set '" << setting_set << "' and setting '" << setting.fName << "'\n";
        std::string tmp = "Possible setting sizes are: ";
        for (const unsigned int& i : possible_sizes){
            tmp+= std::to_string(i) + " ";
        }
        tmp+= ". Please check this!";
        LOG(ERROR) << tmp << "\n";
        return 1;
    }
    //
    // sizes are correct, not we need to check the actual input
    if(CheckParameters(param, possible_settings, setting_set, setting)) return 1;
    return 0;
}

//_______________________________________________________________________________________
//
int ConfigParser::CheckParameters(std::string current, const std::vector<std::string> &possible_settings, const std::string &setting_set, const Config &setting) const{
    if (possible_settings.size() == 1){
        if(!SettingMultipleParamIsOK(setting_set, setting, current, possible_settings.at(0))) return 1;
    }
    else {
        // multiple settings are possible
        bool isFound = false;
        std::transform(current.begin(), current.end(), current.begin(), ::toupper);
        for (const std::string& isetting : possible_settings){
            // for std::string
            if(SettingMultipleParamIsOK(setting_set, setting, current, isetting)){
                isFound = true;
                break;
            }
        }
        if (!isFound){
            LOG(ERROR) << "Error in line " << setting.fLineNumber << " of " << setting.fFileName << ": Parameter " << current <<" is not valid, for setting set '" << setting_set << "' and setting '" << setting.fName << "'\n";
            std::string tmp = "Possible values: ";
            for (const std::string &i : possible_settings){
                tmp+= i + " ";
            }
            tmp+= ". Please check this!";
            LOG(ERROR) << tmp << "\n";
            return 1;
        }
    }
    return 0;
}

//_______________________________________________________________________________________
//
bool ConfigParser::SettingMultipleParamIsOK(const std::string& setting_set, const Config& setting, const std::string& current, const std::string& possible, const char delimiter) const{
    const std::vector<std::string> current_vec = Common::Vectorize(current, delimiter);
    const std::vector<std::string> possible_vec = Common::Vectorize(possible, delimiter);
    if (current_vec.size() != possible_vec.size()) return false;
    //
    // check setting by setting
    for (std::size_t iparam = 0; iparam < possible_vec.size(); iparam++){
        if (possible_vec.at(iparam) == "string"){
            continue; // nothing to check
        } else if (possible_vec.at(iparam) == "int"){
            try {
                // cppcheck-suppress ignoredReturnValue
		std::stoi(current_vec.at(iparam));
            } catch (std::exception &e){
                LOG(ERROR) << "Error in line " << setting.fLineNumber << " of " << setting.fFileName << ": Parameter " << current_vec.at(iparam) << " is not valid, for setting set '" << setting_set << "' and setting '" << setting.fName << "', for parameter number " << (iparam+1) << ". Please check this!\n";
                return false;
            }
        } else if (possible_vec.at(iparam) == "float"){
            try {
                // cppcheck-suppress ignoredReturnValue
		std::stod(current_vec.at(iparam));
            } catch (std::exception &e){
                LOG(ERROR) << "Error in line " << setting.fLineNumber << " of " << setting.fFileName << ": Parameter " << current_vec.at(iparam) << " is not valid, for setting set '" << setting_set << "' and setting '" << setting.fName << "', for parameter number " << (iparam+1) << ". Please check this!\n";
                return false;
            }
        } else {
            if (current_vec.at(iparam) != possible_vec.at(iparam)){
                return false;
            }
        }
    }
    return true;
}
