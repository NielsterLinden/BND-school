#pragma once

#include <map>
#include <string>
#include <vector>

class ReparametrizationManager{
public:
    ReparametrizationManager() = default;
    ~ReparametrizationManager() = default;

    void ReadReparametrizationFile(const std::string& path);

    std::vector<std::pair<std::string, std::string> > GetParametrization(const std::string& region, const std::string& sample) const;
private:

    std::map<std::string, std::vector<std::pair<std::string, std::string> > > fParametrizationMap;
};
