#include "TRExFitter/YamlConverter.h"

#include "TRExFitter/Common.h"
#include "TRExFitter/Logger.h"

#include "yaml-cpp/include/yaml-cpp/yaml.h"

#include "TGraphAsymmErrors.h"
#include "TSystem.h"

#include <algorithm>
#include <ctime>
#include <fstream>

YamlConverter::YamlConverter() :
    m_lumi("139"),
    m_cme("13000") {
}

void YamlConverter::WriteRanking(const std::vector<YamlConverter::RankingContainer>& ranking,
                                 const std::string& path,
                                 const bool shiftedGOs) const {

    YAML::Emitter out;
    out << YAML::BeginSeq;
    for (const auto& irank : ranking) {
        out << YAML::BeginMap;
            out << YAML::Key << "Name";
            out << YAML::Value << irank.name;
            out << YAML::Key << "NPhat";
            out << YAML::Value << irank.nphat;
            out << YAML::Key << "NPerrHi";
            out << YAML::Value << irank.nperrhi;
            out << YAML::Key << "NPerrLo";
            out << YAML::Value << irank.nperrlo;
            out << YAML::Key << "POIup";
            out << YAML::Value << irank.poihi;
            out << YAML::Key << "POIdown";
            out << YAML::Value << irank.poilo;
            if (!shiftedGOs) {
                out << YAML::Key << "POIupPreFit";
                out << YAML::Value << irank.poiprehi;
                out << YAML::Key << "POIdownPreFit";
                out << YAML::Value << irank.poiprelo;
            }
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;

    // Write to the file
    Write(out, "ranking", path);
}

void YamlConverter::WriteImpact(const std::vector<YamlConverter::ImpactContainer>& impact,
                                const std::string& path) const {

    YAML::Emitter out;
    out << YAML::BeginSeq;
    for (const auto& iimpact : impact) {
        out << YAML::BeginMap;
            out << YAML::Key << "Category";
            out << YAML::Value << iimpact.name;
            out << YAML::Key << "Impact";
            out << YAML::Value << iimpact.error;
            out << YAML::Key << "ImpactUp";
            out << YAML::Value << iimpact.errorHi;
            out << YAML::Key << "ImpactDown";
            out << YAML::Value << iimpact.errorLo;
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;

    // Write to the file
    Write(out, "impact", path);
}

void YamlConverter::WriteRankingHEPData(const std::vector<RankingContainer>& ranking,
                                        const std::string& folder,
                                        const std::string& suffix,
                                        const bool shiftedGOs) const {

    gSystem->mkdir((folder+"/HEPData").c_str());

    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "independent_variables";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value <<  "parameter" << YAML::EndMap;
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(irank.name, '#', '\\'));
                    out << YAML::EndMap;
                }
                out << YAML::Value << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;

        // dependent variables
        out << YAML::Key << "dependent_variables";
        out << YAML::Value << YAML::BeginSeq;
            // NP value
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "NP value, error" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    AddValueErrors(out, irank.nphat, irank.nperrhi, irank.nperrlo);
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;

            // NP PosFit impact up
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "impact POI high" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    if (TRExFitter::USEHEPDATAROUNDING) {
                        out << YAML::Value << Common::KeepSignificantDigits(irank.poihi,2);
                    } else {
                        out << YAML::Value << irank.poihi;
                    }
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
            // NP Postfit imapct down
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "impact POI low" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& irank : ranking) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    if (TRExFitter::USEHEPDATAROUNDING) {
                        out << YAML::Value << Common::KeepSignificantDigits(irank.poilo,2);
                    } else {
                        out << YAML::Value << irank.poilo;
                    }
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
            if (!shiftedGOs) {
                // NP PreFit impact up
                out << YAML::BeginMap;
                    out << YAML::Key << "header";
                    out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "impact POI prefit high" << YAML::EndMap;
                    AddQualifiers(out);
                    out << YAML::Key << "values";
                    out << YAML::Value << YAML::BeginSeq;
                    for (const auto& irank : ranking) {
                        out << YAML::BeginMap;
                        out << YAML::Key << "value";
                        if (TRExFitter::USEHEPDATAROUNDING) {
                            out << YAML::Value << Common::KeepSignificantDigits(irank.poiprehi, 2);
                        } else {
                            out << YAML::Value << irank.poiprehi;
                        }
                        out << YAML::EndMap;
                    }
                    out << YAML::EndSeq;
                out << YAML::EndMap;
                // NP Prefit imapct down
                out << YAML::BeginMap;
                    out << YAML::Key << "header";
                    out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "impact POI prefit low" << YAML::EndMap;
                    AddQualifiers(out);
                    out << YAML::Key << "values";
                    out << YAML::Value << YAML::BeginSeq;
                    for (const auto& irank : ranking) {
                        out << YAML::BeginMap;
                        out << YAML::Key << "value";
                        if (TRExFitter::USEHEPDATAROUNDING) {
                            out << YAML::Value << Common::KeepSignificantDigits(irank.poiprelo,2);
                        } else {
                            out << YAML::Value << irank.poiprelo;
                        }
                        out << YAML::EndMap;
                    }
                    out << YAML::EndSeq;
                out << YAML::EndMap;
            }
        out << YAML::EndSeq;
    out << YAML::EndMap;

    // Write to the file
    Write(out, "HEPData ranking", folder + "/HEPData/Ranking"+suffix+".yaml");
}

void YamlConverter::WriteImpactHEPData(const std::vector<ImpactContainer>& impact,
                                       const std::string& folder,
                                       const std::string& suffix) const {

    gSystem->mkdir((folder+"/HEPData").c_str());

    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "independent_variables";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "category" << YAML::EndMap;
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& iimpact : impact) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(iimpact.name, '#', '\\'));
                    out << YAML::EndMap;
                }
                out << YAML::Value << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;

        // dependent variables
        out << YAML::Key << "dependent_variables";
        out << YAML::Value << YAML::BeginSeq;
            // Impact value
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Grouped uncertainty" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& iimpact : impact) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    if (TRExFitter::USEHEPDATAROUNDING) {
                        out << YAML::Value << Common::KeepSignificantDigits(iimpact.error,2);
                    } else {
                        out << YAML::Value << iimpact.error;
                    }
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;
    out << YAML::EndMap;

    // Write to the file
    Write(out, "HEPData impact", folder + "/HEPData/GroupedBreakdown"+suffix+".yaml");
}

void YamlConverter::AddQualifiers(YAML::Emitter& out) const {
    out << YAML::Key << "qualifiers";
    out << YAML::Value << YAML::BeginSeq;
    out << YAML::BeginMap;
        out << YAML::Key << "name";
        out << YAML::Value << "SQRT(s)";
        out << YAML::Key << "units";
        out << YAML::Value << "GeV";
        out << YAML::Key << "value";
        out << YAML::Value << Common::convertStoNum<int>(m_cme);
    out << YAML::EndMap;
    out << YAML::BeginMap;
        out << YAML::Key << "name";
        out << YAML::Value << "LUMINOSITY";
        out << YAML::Key << "units";
        out << YAML::Value << "fb$^{-1}$";
        out << YAML::Key << "value";
        out << YAML::Value << Common::convertStoNum<int>(m_lumi);
    out << YAML::EndMap;
    out << YAML::EndSeq;
}

void YamlConverter::AddValueErrors(YAML::Emitter& out,
                                   const double mean,
                                   const double up,
                                   const double down) const {

    double value = mean;
    if ((std::abs(down) > 1e-6) && (std::abs(up/down) > 0.9) && (std::abs(up/down) < 1.1)) {
        double error = 0.5*(std::abs(up) + std::abs(down)) ;
        int n(0);
        if (TRExFitter::USEHEPDATAROUNDING) {
            n = Common::ApplyPDGrounding(value, error);
        }
        // are symmetric
        out << YAML::Key << "value";
        if (TRExFitter::USEHEPDATAROUNDING) {
            out << YAML::Value << Form(("%."+std::to_string(n)+"f").c_str(),value);
        } else {
            out << YAML::Value << value;
        }
        out << YAML::Key << "errors";
        out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "symerror";
        if (TRExFitter::USEHEPDATAROUNDING) {
            if (n >= 0) {
                out << YAML::Value << Form(("%."+std::to_string(n)+"f").c_str(),error);
            } else {
                out << YAML::Value << Form("%.f",error);
            }
        } else {
            out << YAML::Value << error;
        }
        out << YAML::EndMap;
        out << YAML::EndSeq;
    } else {
        int n(0);
        if (TRExFitter::USEHEPDATAROUNDING) {
            double error = std::min(std::abs(up),std::abs(down));
            n = Common::ApplyPDGrounding(value, error);
        }
        out << YAML::Key << "value";
        if (TRExFitter::USEHEPDATAROUNDING) {
            out << YAML::Value << Form(("%."+std::to_string(n)+"f").c_str(),value);
        } else {
            out << YAML::Value << value;
        }
        out << YAML::Key << "errors";
        out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "asymerror";
        out << YAML::Value << YAML::BeginMap;
            out << YAML::Key << "plus";
            if (TRExFitter::USEHEPDATAROUNDING) {
                if (n >= 0) {
                    out << YAML::Key << Form(("%."+std::to_string(n)+"f").c_str(),up);
                } else {
                    out << YAML::Key << Form("%.f",up);
                }
            } else {
                out << YAML::Key << up;
            }
            out << YAML::Key << "minus";
            if (TRExFitter::USEHEPDATAROUNDING) {
                if (n >= 0) {
                    out << YAML::Key << Form(("%."+std::to_string(n)+"f").c_str(),down);
                } else {
                    out << YAML::Key << Form("%.f",down);
                }
            } else {
                out << YAML::Key << down;
            }
        out << YAML::EndMap;
        out << YAML::EndMap;
        out << YAML::EndSeq;
    }
}


void YamlConverter::Write(const YAML::Emitter& out, const std::string& type, const std::string& path) const {
    LOG(INFO) << "Writing " << type << " yaml file to: " << path << "\n";
    std::ofstream file;
    file.open(path.c_str());
    if (!file.is_open() || !file.good()) {
        LOG(ERROR) << "Cannot open yaml file at: " << path << "\n";
        return;
    }

    file << out.c_str();
    file.close();
}

void YamlConverter::WriteCorrelation(const std::vector<std::string>& np,
                                     const TMatrixD& corr,
                                     const std::string& path,
                                     const bool isCovariance,
                                     const std::string& name) const {

    const int n = corr.GetNrows();
    std::vector<std::vector<double> > corVector(n, std::vector<double>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            corVector[i][j] = corr(i,j);
        }
    }

    this->WriteCorrelation(np, corVector, path, isCovariance, name);
}

void YamlConverter::WriteCorrelation(const std::vector<std::string>& np,
                                     const std::vector<std::vector<double> >& corr,
                                     const std::string& path,
                                     const bool isCovariance,
                                     const std::string& name) const {

    const std::size_t n = np.size();
    if (corr.size() != n) {
        LOG(ERROR) << "Inconsistent inputs!\n";
        return;
    }

    YAML::Emitter out;
    out << YAML::BeginSeq;
    out << YAML::BeginMap;
    out << YAML::Key << "parameters";
    out << YAML::Value << YAML::BeginSeq;
    for (const auto& iname : np) {
        out << iname;
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;
    out << YAML::BeginMap;
    out << YAML::Key << (isCovariance ? "covariance_rows" : "correlation_rows");
    out << YAML::Value << YAML::BeginSeq;
    for (const auto& icorr : corr) {
        if (icorr.size() != n) {
            LOG(ERROR) << "Inconsistent inputs for the matrix!\n";
            return;
        }
        out << YAML::Flow << YAML::BeginSeq;
        for (const auto& i : icorr) {
            out << Form("%2.2e",i);
        }
        out << YAML::EndSeq;
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;
    out << YAML::EndSeq;

    Write(out, isCovariance ? "covariance" : "correlation", path+"/"+name+".yaml");
}

void YamlConverter::WriteCorrelationHEPData(const std::vector<std::string>& np,
                                            const TMatrixD& corr,
                                            const std::string& folder,
                                            const bool isCovariance,
                                            const std::string& name) const {

    const int n = corr.GetNrows();
    std::vector<std::vector<double> > corVector(n, std::vector<double>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            corVector[i][j] = corr(i,j);
        }
    }

    this->WriteCorrelationHEPData(np, corVector, folder, isCovariance, name);
}

void YamlConverter::WriteCorrelationHEPData(const std::vector<std::string>& np,
                                            const std::vector<std::vector<double> >& corr,
                                            const std::string& folder,
                                            const bool isCovariance,
                                            const std::string& name) const {


    gSystem->mkdir((folder+"/HEPData").c_str());

    const std::size_t n = np.size();
    if (corr.size() != n) {
        LOG(ERROR) << "Inconsistent inputs!\n";
        return;
    }

    if (n > 100) {
        LOG(INFO) << "Processing matrix of more than 100x100 elements, this may take some time...\n";
    }

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "dependent_variables";
    out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "header";
        out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << (isCovariance ? "parameter covariances" : "NP correlations") << YAML::EndMap;
        AddQualifiers(out);
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (const auto& icorr : corr) {
            if (icorr.size() != n) {
                LOG(ERROR) << "Inconsistent inputs for correlation!\n";
                return;
            }
            for (const auto& i : icorr) {
                out << YAML::BeginMap;
                out << YAML::Key << "value";
                out << YAML::Value << Form("%2.2e", i);
                out << YAML::EndMap;
            }
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    out << YAML::EndSeq;
    out << YAML::Key << "independent_variables";
    out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "header";
        out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << (isCovariance ? "parameters" : "NPs") << YAML::EndMap;
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (std::size_t inp = 0; inp < n; ++inp) {
            for (std::size_t jnp = 0; jnp < n; ++jnp) {
                out << YAML::BeginMap;
                out << YAML::Key << "value";
                out << YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(np.at(inp), '#', '\\'));
                out << YAML::EndMap;
            }
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
        out << YAML::BeginMap;
        out << YAML::Key << "header";
        out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << (isCovariance ? "parameters" : "NPs") << YAML::EndMap;
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (std::size_t inp = 0; inp < n; ++inp) {
            for (std::size_t jnp = 0; jnp < n; ++jnp) {
                out << YAML::BeginMap;
                out << YAML::Key << "value";
                out << YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(np.at(jnp), '#', '\\'));
                out << YAML::EndMap;
            }
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    out << YAML::EndSeq;
    out << YAML::EndMap;

    Write(out, isCovariance ? "HEPData covariance" : "HEPData correlation", folder+"/HEPData/"+name+".yaml");

}

void YamlConverter::WriteTables(const YamlConverter::TableContainer& container,
                                const std::string& directory,
                                const bool isPostFit) const {

    if (!YamlConverter::TableContainerIsOK(container)) {
        LOG(ERROR) << "Inconsistent inputs for tables!\n";
        return;
    }
    gSystem->mkdir((directory+"/HEPData").c_str());

    YAML::Emitter out;
    out << YAML::BeginSeq;
    for (std::size_t ireg = 0; ireg < container.regionNames.size(); ++ireg) {
        out << YAML::BeginMap;
        out << YAML::Key << "Region";
        out << YAML::Value << container.regionNames.at(ireg);
        out << YAML::Key << "Samples";
        out << YAML::Value << YAML::BeginSeq;
        for (std::size_t isample = 0; isample < container.sampleNames.size(); ++isample) {
            out << YAML::BeginMap;
                out << YAML::Key << "Sample";
                out << YAML::Value << container.sampleNames.at(isample);
                out << YAML::Key << "Yield";
                out << YAML::Value << container.mcYields.at(isample).at(ireg);
                out << YAML::Key << "Error";
                out << YAML::Value << container.mcErrors.at(isample).at(ireg);
            out << YAML::EndMap;
        }
        for (std::size_t idata = 0; idata < container.dataNames.size(); ++idata) {
            out << YAML::BeginMap;
                out << YAML::Key << "Data";
                out << YAML::Value << container.dataNames.at(idata);
                out << YAML::Key << "Yield";
                out << YAML::Value << container.dataYields.at(idata).at(ireg);
            out << YAML::EndMap;
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;
    // Write to the file
    if (isPostFit) {
        Write(out, "postfit yield tables", directory + "/Tables/Table_postfit.yaml");
    } else {
        Write(out, "prefit yield tables", directory + "/Tables/Table_prefit.yaml");
    }
}


void YamlConverter::WriteTablesHEPData(const YamlConverter::TableContainer& container,
                                       const std::string& directory,
                                       const bool isPostFit) const {

    if (!YamlConverter::TableContainerIsOK(container)) {
        LOG(ERROR) << "Inconsistent inputs for tables!\n";
        return;
    }

    gSystem->mkdir((directory+"/HEPData").c_str());

    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "independent_variables";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value <<  "process" << YAML::EndMap;
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& isample : container.sampleNames) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(isample, '#', '\\'));
                    out << YAML::EndMap;
                }
                for (const auto& idata : container.dataNames) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << idata;
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;

        // dependent variables
        out << YAML::Key << "dependent_variables";
        out << YAML::Value << YAML::BeginSeq;
        // loop over regions
        for (std::size_t ireg = 0; ireg < container.regionNames.size(); ++ireg) {
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(container.regionNames.at(ireg), '#', '\\')) << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (std::size_t ivalue = 0; ivalue < container.mcYields.size(); ++ivalue) {
                    out << YAML::BeginMap;
                    AddValueErrors(out, container.mcYields.at(ivalue).at(ireg), container.mcErrors.at(ivalue).at(ireg), container.mcErrors.at(ivalue).at(ireg));
                    out << YAML::EndMap;
                }
                for (std::size_t ivalue = 0; ivalue < container.dataYields.size(); ++ivalue) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << Form("%.f",container.dataYields.at(ivalue).at(ireg));
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        }
        out << YAML::EndSeq;
    out << YAML::EndMap;

    // Write to the file
    if (isPostFit) {
        Write(out, "HEPData postfit yield tables", directory + "/HEPData/Table_postfit.yaml");
    } else {
        Write(out, "HEPData prefit yield tables", directory + "/HEPData/Table_prefit.yaml");
    }
}

bool YamlConverter::TableContainerIsOK(const YamlConverter::TableContainer& container) const {

    const std::size_t nRegions = container.regionNames.size();
    const std::size_t nSamples = container.sampleNames.size();

    if (nRegions == 0 || nSamples == 0) return false;
    if (container.mcYields.size() != nSamples) return false;
    if (container.mcErrors.size() != nSamples) return false;
    for (const auto& ivec : container.mcYields) {
        if (ivec.size() != nRegions) return false;
    }
    for (const auto& ivec : container.mcErrors) {
        if (ivec.size() != nRegions) return false;
    }
    for (const auto& ivec : container.dataYields) {
        if (ivec.size() != nRegions) return false;
    }

    return true;
}

void YamlConverter::WriteUnfolding(const TGraphAsymmErrors* const graph,
                                   const std::string& directory,
                                   const std::string& name) const {

    const int n = graph->GetN();
    YAML::Emitter out;
    out << YAML::BeginSeq;
    for (int i = 0; i < n; ++i) {
        double x;
        double y;
        graph->GetPoint(i, x, y);

        const double x_min = graph->GetErrorXlow(i);
        const double x_max = graph->GetErrorXhigh(i);
        const double y_min = graph->GetErrorYlow(i);
        const double y_max = graph->GetErrorYhigh(i);
        out << YAML::BeginMap;
            out << YAML::Key << "range";
            out << YAML::Value << YAML::Flow << YAML::BeginSeq << x-x_min << x+x_max << YAML::EndSeq;
            out << YAML::Key << "mean";
            out << YAML::Value << y;
            out << YAML::Key << "uncertaintyUp";
            out << YAML::Value << y_max;
            out << YAML::Key << "uncertaintyDown";
            out << YAML::Value << -y_min;
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;

    Write(out, "unfolding result", directory + "/"+name+"_UnfoldingData.yaml");
}

void YamlConverter::WriteUnfoldingHEPData(const TGraphAsymmErrors* const graph,
                                          const std::string& xAxis,
                                          const std::string& directory,
                                          const std::string& name,
                                          const std::vector<std::unique_ptr<TH1> >& mc,
                                          const std::vector<std::string>& legendNames) const {

    gSystem->mkdir((directory + "/HEPData").c_str());
    const int n = graph->GetN();

    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "independent_variables";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap;
                out << YAML::Key << "name" << YAML::Value <<  YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(xAxis, '#', '\\'));
                out << YAML::EndMap;
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (int i = 0; i < n; ++i) {
                    double x;
                    double y;
                    graph->GetPoint(i, x, y);

                    const double x_min = graph->GetErrorXlow(i);
                    const double x_max = graph->GetErrorXhigh(i);
                    out << YAML::BeginMap;
                        out << YAML::Key << "high";
                        out << YAML::Value << Common::KeepSignificantDigits(x+x_max, 4);
                        out << YAML::Key << "low";
                        out << YAML::Value << Common::KeepSignificantDigits(x-x_min, 4);
                        out << YAML::Key << "value";
                        out << YAML::Value << Common::KeepSignificantDigits(x, 4);
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;

        // dependent variables
        out << YAML::Key << "dependent_variables";
        out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
            out << YAML::Key << "header";
            out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Data" << YAML::EndMap;
            AddQualifiers(out);
            out << YAML::Key << "values";
            out << YAML::Value << YAML::BeginSeq;
            for (int i = 0; i < n; ++i) {
                double x;
                double y;
                graph->GetPoint(i, x, y);

                const double y_min = graph->GetErrorYlow(i);
                const double y_max = graph->GetErrorYhigh(i);
                out << YAML::BeginMap;
                AddValueErrors(out, y, y_max, -y_min);
                out << YAML::EndMap;
            }
            out << YAML::EndSeq;
        out << YAML::EndMap;
        // MC
        for (std::size_t imc = 0; imc < mc.size(); ++imc) {
            out << YAML::BeginMap;
                out << YAML::Key;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << legendNames.at(imc).c_str() << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (int i = 1; i <= n; ++i) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << mc.at(imc)->GetBinContent(i);
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
                out << YAML::EndMap;
        }
        out << YAML::EndSeq;
    out << YAML::EndMap;

    Write(out, "HEPData unfolding result", directory + "/HEPData/"+name+"_Unfolding.yaml");
}

void YamlConverter::WritePlot(const YamlConverter::PlotContainer& container,
                              const std::string& directory,
                              const std::string& subdirectory,
                              const bool isPostFit) const {

    if (!YamlConverter::PlotContainerIsOK(container)) {
        LOG(ERROR) << "Inconsistent inputs for plots!\n";
        return;
    }

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "Samples";
    out << YAML::Value << YAML::BeginSeq;
    for (std::size_t isample = 0; isample < container.signalYields.size(); ++isample) {
        out << YAML::BeginMap;
        out << YAML::Key << "Name";
        out << YAML::Value << container.samples.at(isample);
        out << YAML::Key << "Yield";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        for (const auto& i : container.signalYields.at(isample)) {
            out << i;
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    }
    for (std::size_t isample = 0; isample < container.backgroundYields.size(); ++isample) {
        out << YAML::BeginMap;
        out << YAML::Key << "Name";
        out << YAML::Value << container.samples.at(container.signalYields.size() + isample);
        out << YAML::Key << "Yield";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        for (const auto& i : container.backgroundYields.at(isample)) {
            out << i;
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;

    out << YAML::Key << "Total";
    out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "Yield";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        for (int i = 0; i < container.errors->GetN(); ++i) {
            double x;
            double y;
            container.errors->GetPoint(i, x, y);
            out << y;
        }
        out << YAML::EndSeq;
        out << YAML::Key << "UncertaintyUp";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        for (int i = 0; i < container.errors->GetN(); ++i) {
            double x;
            double y;
            container.errors->GetPoint(i, x, y);
            const double y_max = container.errors->GetErrorYhigh(i);
            out << y_max;
        }
        out << YAML::EndSeq;
        out << YAML::Key << "UncertaintyDown";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        for (int i = 0; i < container.errors->GetN(); ++i) {
            double x;
            double y;
            container.errors->GetPoint(i, x, y);
            const double y_min = container.errors->GetErrorYlow(i);
            out << -y_min;
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;

    out << YAML::EndSeq;

    if (!container.data.empty()) {
        out << YAML::Key << "Data";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
            out << YAML::Key << "Yield";
            out << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (std::size_t i = 0; i < container.data.size(); i++) {
                if (std::find(container.blindedBins.begin(), container.blindedBins.end(), i+1) == container.blindedBins.end()) {
                    out << container.data.at(i);
                } else {
                    out << "blinded";
                }
            }
            out << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;
    }

    out << YAML::Key << "Figure";
    out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "BinEdges";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        double x;
        double y;
        for (int i = 0; i < container.errors->GetN(); ++i) {
            container.errors->GetPoint(i, x, y);
            const double x_min = container.errors->GetErrorXlow(i);
            out << x - x_min;
        }
        const double x_max = container.errors->GetErrorXlow(container.errors->GetN() - 1);
        out << x + x_max;
        out << YAML::EndSeq;
        out << YAML::Key << "XaxisLabel";
        out << YAML::Value << container.xAxis;
        out << YAML::Key << "YaxisLabel";
        out << YAML::Value << container.yAxis;
        out << YAML::EndMap;
    out << YAML::EndSeq;
    out << YAML::EndMap;


    std::string outputDirectory("");
    if (subdirectory != ""){
        outputDirectory = directory + "/Plots/" + subdirectory + "/";
    } else {
        outputDirectory = directory + "/Plots/";
    }
    if (isPostFit) {
        const std::string path = outputDirectory + container.region + "_postfit.yaml";
        Write(out, "postfit plot " + container.region, path);
    } else {
        const std::string path = outputDirectory + container.region + "_prefit.yaml";
        Write(out, "prefit plot " + container.region, path);
    }
}

void YamlConverter::WritePlotUhepp(const YamlConverter::PlotContainer& container,
                                   const std::string& directory,
                                   const bool isPostFit) const {

    if (!YamlConverter::PlotContainerIsOK(container)) {
        LOG(ERROR) << "Inconsistent inputs for plots!\n";
        return;
    }

    YAML::Emitter out;
    out << YAML::BeginMap;

    // General uhepp header
    out << YAML::Key << "version";
    out << YAML::Value << "0.4";
    out << YAML::Key << "type";
    out << YAML::Value << "histogram";

    // Metadata
    out << YAML::Key << "metadata";
    out << YAML::BeginMap;

    out << YAML::Key << "filename";
    out << YAML::Value << container.region;

    out << YAML::Key << "date";
    time_t now;
    time(&now);
    char buf[sizeof "2011-10-08T07:07:09Z"];
    strftime(buf, sizeof buf, "%FT%TZ", gmtime(&now));
    out << YAML::Value << buf;

    out << YAML::Key << "producer";
    out << YAML::Value << "TRExFitter";

    out << YAML::Key << "tags";
    out << YAML::BeginMap;
    if (isPostFit) {
      out << YAML::Key << "post-fit";
      out << YAML::Null;
    } else {
      out << YAML::Key << "pre-fit";
      out << YAML::Null;
    }
    out << YAML::EndMap;

    out << YAML::Key << "Ecm_TeV";
    out << YAML::Value << Common::convertStoNum<double>(m_cme) / 1000;

    out << YAML::Key << "lumi_ifb";
    out << YAML::Value << Common::convertStoNum<double>(m_lumi);
    out << YAML::EndMap;

    // Variable
    out << YAML::Key << "variable";
    out << YAML::BeginMap;
    out << YAML::Key << "symbol";
    out << YAML::Value << std::string("$") + container.xAxis + "$";
    out << YAML::EndMap;

    // Bins
    out << YAML::Key << "bins";
    out << YAML::BeginMap;
    out << YAML::Key << "edges";
    out << YAML::Value << YAML::Flow << YAML::BeginSeq;
    double x;
    double y;
    for (int i = 0; i < container.errors->GetN(); ++i) {
        container.errors->GetPoint(i, x, y);
        const double x_min = container.errors->GetErrorXlow(i);
        out << x - x_min;
    }
    const double x_max = container.errors->GetErrorXlow(container.errors->GetN() - 1);
    out << x + x_max;
    out << YAML::EndSeq;
    out << YAML::EndMap;

    std::map<std::string,std::vector<double>> yields;
    for (std::size_t isample = 0; isample < container.signalYields.size(); ++isample) {
        std::string sample_name = container.samples.at(isample);
        const auto& sample = container.signalYields.at(isample);
        if (yields.count(sample_name)) {
            // sample exists
            std::vector<double>& sample_yield = yields[sample_name];
            for(std::size_t i = 0; i < yields[sample_name].size(); ++i) {
              sample_yield[i] += sample.at(i);
            }
        } else {
            // New sample
            std::vector<double>& sample_yield = yields[sample_name];
            for (const auto& i : sample) {
             sample_yield.push_back(i);
            }
        }
    }
    for (std::size_t isample = 0; isample < container.backgroundYields.size(); ++isample) {
        std::string sample_name = container.samples.at(container.signalYields.size() + isample);
        const auto& sample = container.backgroundYields.at(isample);
        if (yields.count(sample_name)) {
            // sample exists
            std::vector<double>& sample_yield = yields[sample_name];
            for(std::size_t i = 0; i < yields[sample_name].size(); ++i) {
              sample_yield[i] += sample.at(i);
            }
        } else {
            // New sample
            std::vector<double>& sample_yield = yields[sample_name];
            for (const auto& i : sample) {
              sample_yield.push_back(i);
            }
        }
    }

    out << YAML::Key << "yields";
    out << YAML::Value << YAML::BeginMap;
    for (auto const& sample : yields) {
        out << YAML::Key << sample.first;
        out << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "base";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        out << 0; // Underflow
        for (size_t i = 0; i < sample.second.size(); i++) {
          out << sample.second[i];
        }
        out << 0; // Overflow
        out << YAML::Value << YAML::EndSeq;
        out << YAML::Key << "stat";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        out << 0; // Underflow
        for (size_t i = 0; i < sample.second.size(); i++) {
          out << 0;
        }
        out << 0; // Overflow
        out << YAML::Value << YAML::EndSeq;
        out << YAML::Value << YAML::EndMap;
    }

    // Total
    out << YAML::Key << "trex_total";
    out << YAML::Value << YAML::BeginMap;

    out << YAML::Key << "base";
    out << YAML::Value << YAML::Flow << YAML::BeginSeq;
    out << 0; // Underflow
    for (int i = 0; i < container.errors->GetN(); ++i) {
        container.errors->GetPoint(i, x, y);
        out << y;
    }
    out << 0; // Overflow
    out << YAML::EndSeq;


    out << YAML::Key << "var_up";
    out << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "total";
    out << YAML::Value << YAML::Flow << YAML::BeginSeq;
    out << 0; // Underflow
    for (int i = 0; i < container.errors->GetN(); ++i) {
        container.errors->GetPoint(i, x, y);
        const double y_max = container.errors->GetErrorYhigh(i);
        out << y + y_max;
    }
    out << 0; // Overflow
    out << YAML::EndSeq;
    out << YAML::EndMap;

    out << YAML::Key << "var_down";
    out << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "total";
    out << YAML::Value << YAML::Flow << YAML::BeginSeq;
    out << 0; // Underflow
    for (int i = 0; i < container.errors->GetN(); ++i) {
        container.errors->GetPoint(i, x, y);
        const double y_min = container.errors->GetErrorYlow(i);
        out << y - y_min;
    }
    out << 0; // Overflow
    out << YAML::EndSeq;
    out << YAML::EndMap; // var_down

    out << YAML::EndMap; // trex_total


    if (!container.data.empty()) {
        out << YAML::Key << "trex_data";
        out << YAML::BeginMap;

        out << YAML::Key << "base";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        out << 0;  // Underflow
        for (std::size_t i = 0; i < container.data.size(); i++) {
            if (std::find(container.blindedBins.begin(), container.blindedBins.end(), i+1) == container.blindedBins.end()) {
                out << container.data.at(i);
            } else {
                out << 0;
            }
        }
        out << 0;  // Overflow
        out << YAML::EndSeq;


        out << YAML::Key << "stat";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq;
        out << 0;  // Underflow
        for (std::size_t i = 0; i < container.data.size(); i++) {
            if (std::find(container.blindedBins.begin(), container.blindedBins.end(), i+1) == container.blindedBins.end()) {
                // TODO: Implement different error algorithm
                out << std::sqrt(container.data.at(i));
            } else {
                out << 0;
            }
        }
        out << 0;  // Overflow
        out << YAML::EndSeq;


        out << YAML::EndMap;  // Data
    }

    out << YAML::EndMap; // yields

    // Stacks
    out << YAML::Key << "stacks";
    out << YAML::Value << YAML::BeginSeq;

    // Background
    out << YAML::BeginMap;
    out << YAML::Key << "content";

    out << YAML::Value << YAML::BeginSeq;
    for (const auto& sample : yields) {
      out << YAML::BeginMap;
      out << YAML::Key << "yield";
      out << YAML::Value << YAML::Flow << YAML::BeginSeq << sample.first << YAML::EndSeq;
      out << YAML::Key << "label";
      out << YAML::Value << sample.first;
      out << YAML::EndMap;
    }
    out << YAML::EndSeq;;

    out << YAML::Key << "error";
    out << YAML::Value << YAML::DoubleQuoted << "no";
    out << YAML::Key << "type";
    out << YAML::Value << "stepfilled";
    out << YAML::EndMap;  // Background stack

    // Total
    out << YAML::BeginMap;
    out << YAML::Key << "content";

    out << YAML::Value << YAML::BeginSeq;
    out << YAML::BeginMap;
    out << YAML::Key << "yield";
    out << YAML::Value << YAML::Flow << YAML::BeginSeq << "trex_total" << YAML::EndSeq;
    out << YAML::Key << "label";
    out << YAML::Value << "Total";
    out << YAML::EndMap;
    out << YAML::EndSeq;

    out << YAML::Key << "error";
    out << YAML::Value << "env";
    out << YAML::Key << "type";
    out << YAML::Value << "step";
    out << YAML::EndMap;  // total stack

    // Data
    out << YAML::BeginMap;
    out << YAML::Key << "content";

    out << YAML::Value << YAML::BeginSeq;
    out << YAML::BeginMap;
    out << YAML::Key << "yield";
    out << YAML::Value << YAML::Flow << YAML::BeginSeq << "trex_data" << YAML::EndSeq;
    out << YAML::Key << "label";
    out << YAML::Value << "Data";
    out << YAML::EndMap;
    out << YAML::EndSeq;

    out << YAML::Key << "error";
    out << YAML::Value << "stat";
    out << YAML::Key << "type";
    out << YAML::Value << "points";
    out << YAML::EndMap;  // data stack

    out << YAML::EndSeq; // stacks

    // Ratio
    if (!container.data.empty()) {
      out << YAML::Key << "ratio";
      out << YAML::Value << YAML::BeginSeq;

      out << YAML::BeginMap;
      out << YAML::Key << "numerator";
      out << YAML::Value << YAML::Flow << YAML::BeginSeq << "trex_data" << YAML::EndSeq;
      out << YAML::Key << "denominator";
      out << YAML::Value << YAML::Flow << YAML::BeginSeq << "trex_total" << YAML::EndSeq;

      out << YAML::Key << "error";
      out << YAML::Value << "stat";
      out << YAML::Key << "den_error";
      out << YAML::Value << "env";

      out << YAML::Key << "type";
      out << YAML::Value << "points";

      out << YAML::EndMap;
      out << YAML::EndSeq;

      out << YAML::Key << "ratio_axis";
      out << YAML::Value << "Data / Total";
    }

    out << YAML::EndMap; // document

    if (isPostFit) {
        const std::string path = directory + "/Plots/" + container.region + "_postfit.uhepp.yaml";
        Write(out, "postfit uhepp plot " + container.region, path);
    } else {
        const std::string path = directory + "/Plots/" + container.region + "_prefit.uhepp.yaml";
        Write(out, "prefit uhepp plot " + container.region, path);
    }
}

void YamlConverter::WritePlotHEPData(const YamlConverter::PlotContainer& container,
                                     const std::string& directory,
                                     const bool isPostFit) const {

    if (!YamlConverter::PlotContainerIsOK(container)) {
        LOG(ERROR) << "Inconsistent inputs for plots!\n";
        return;
    }

    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "independent_variables";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap;
                out << YAML::Key << "name" << YAML::Value << YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(container.xAxis, '#', '\\'));
                out << YAML::EndMap;
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (int i = 0; i < container.errors->GetN(); ++i) {
                    double x;
                    double y;
                    container.errors->GetPoint(i, x, y);

                    const double x_min = container.errors->GetErrorXlow(i);
                    const double x_max = container.errors->GetErrorXhigh(i);
                    out << YAML::BeginMap;
                        out << YAML::Key << "high";
                        out << YAML::Value << Common::KeepSignificantDigits(x+x_max, 4);
                        out << YAML::Key << "low";
                        out << YAML::Value << Common::KeepSignificantDigits(x-x_min, 4);
                        out << YAML::Key << "value";
                        out << YAML::Value << Common::KeepSignificantDigits(x, 4);
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;

        // dependent variables
        out << YAML::Key << "dependent_variables";
        out << YAML::Value << YAML::BeginSeq;
        for (std::size_t isample = 0; isample < container.signalYields.size(); ++isample) {
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value <<
                    YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(container.samples.at(isample), '#', '\\')) << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& i : container.signalYields.at(isample)) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    if (TRExFitter::USEHEPDATAROUNDING) {
                        out << YAML::Value << Common::KeepSignificantDigits(i,2);
                    } else {
                        out << YAML::Value << i;
                    }
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        }
        for (std::size_t isample = 0; isample < container.backgroundYields.size(); ++isample) {
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value <<
                    YamlConverter::AddLatexDelimiters(Common::ReplaceAllCharacters(container.samples.at(container.signalYields.size() + isample), '#', '\\')) << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& i : container.backgroundYields.at(isample)) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    if (TRExFitter::USEHEPDATAROUNDING) {
                        out << YAML::Value << Common::KeepSignificantDigits(i,2);
                    } else {
                        out << YAML::Value << i;
                    }
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        }
        out << YAML::BeginMap;
            out << YAML::Key << "header";
            out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << "Total" << YAML::EndMap;
            AddQualifiers(out);
            out << YAML::Key << "values";
            out << YAML::Value << YAML::BeginSeq;
            for (int i = 0; i < container.errors->GetN(); ++i) {
                double x;
                double y;
                container.errors->GetPoint(i, x, y);

                const double y_min = container.errors->GetErrorYlow(i);
                const double y_max = container.errors->GetErrorYhigh(i);
                out << YAML::BeginMap;
                AddValueErrors(out, y, y_max, -y_min);
                out << YAML::EndMap;
            }
            out << YAML::EndSeq;
        out << YAML::EndMap;

        if (!container.data.empty()) {
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << "Data" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (std::size_t i = 0; i < container.data.size(); i++) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    if (std::find(container.blindedBins.begin(), container.blindedBins.end(), i+1) == container.blindedBins.end()) {
                        out << YAML::Value << Form("%.f", container.data.at(i));
                    } else {
                        out << "blinded";
                    }
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        }

        out << YAML::EndSeq;
    out << YAML::EndMap;


    if (isPostFit) {
        const std::string path = directory + "/HEPData/" + container.region + "_postfit.yaml";
        Write(out, "HEPData postfit plot " + container.region, path);
    } else {
        const std::string path = directory + "/HEPData/" + container.region + "_prefit.yaml";
        Write(out, "HEPData prefit plot " + container.region, path);
    }
}

void YamlConverter::WriteLikelihoodScan(const std::pair<std::vector<double>, std::vector<double> >& result,
                                        const std::string& path) const {

    if (result.first.size() != result.second.size()) {
        LOG(ERROR) << "Size of X and Y do not match\n";
        return;
    }

    YAML::Emitter out;
    out << YAML::BeginSeq;
    for (std::size_t i = 0; i < result.first.size(); ++i) {
        out << YAML::BeginMap;
            out << YAML::Key << "X";
            out << YAML::Value << result.first.at(i);
            out << YAML::Key << "minusdeltaNLL";
            out << YAML::Value << result.second.at(i);
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;

    // Write to the file
    Write(out, "LikelihoodScan", path);
}

void YamlConverter::WriteLikelihoodScan2D(const LikelihoodScanManager::Result2D& result,
                                          const std::string& path) const {

    if (result.x.size() != result.y.size()) {
        LOG(ERROR) << "Size of X and Y do not match\n";
        return;
    }

    if (result.x.size() != result.z.size()) {
        LOG(ERROR) << "Size of X and Z do not match\n";
        return;
    }

    YAML::Emitter out;
    out << YAML::BeginSeq;
    for (std::size_t i = 0; i < result.x.size(); ++i) {
        for (std::size_t j = 0; j < result.x.size(); ++j) {
            out << YAML::BeginMap;
                out << YAML::Key << "X";
                out << YAML::Value << result.x.at(i);
                out << YAML::Key << "Y";
                out << YAML::Value << result.y.at(j);
                out << YAML::Key << "minusdeltaNLL";
                out << YAML::Value << result.z.at(i).at(j);
            out << YAML::EndMap;
        }
    }
    out << YAML::EndSeq;

    // Write to the file
    Write(out, "LikelihoodScan2D", path);
}

void YamlConverter::WriteLikelihoodScanHEPData(const std::pair<std::vector<double>, std::vector<double> >& result,
                                               const std::string& folder,
                                               const std::string& suffix) const {

    gSystem->mkdir((folder+"/HEPData").c_str());

    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "independent_variables";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "X" << YAML::EndMap;
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& x : result.first) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    out << YAML::Value << x;
                    out << YAML::EndMap;
                }
                out << YAML::Value << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;

        // dependent variables
        out << YAML::Key << "dependent_variables";
        out << YAML::Value << YAML::BeginSeq;
            //  value
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "2deltaNLL" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& y : result.second) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "value";
                    if (TRExFitter::USEHEPDATAROUNDING) {
                        out << YAML::Value << Common::KeepSignificantDigits(y,3);
                    } else {
                        out << YAML::Value << y;
                    }
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;
    out << YAML::EndMap;

    // Write to the file
    Write(out, "HEPData LikelihoodScan", folder + "/HEPData/LikelihoodScan_"+suffix+".yaml");
}

void YamlConverter::WriteLikelihoodScan2DHEPData(const LikelihoodScanManager::Result2D& result,
                                                 const std::string& folder,
                                                 const std::string& suffix) const {

    gSystem->mkdir((folder+"/HEPData").c_str());

    YAML::Emitter out;
    out << YAML::BeginMap;
        out << YAML::Key << "independent_variables";
        out << YAML::Value << YAML::BeginSeq;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "X/Y values" << YAML::EndMap;
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (std::size_t i = 0; i < result.x.size(); ++i) {
                    for (std::size_t j = 0; j < result.x.size(); ++j) {
                        out << YAML::BeginMap;
                        out << YAML::Key << "value";
                        out << YAML::Value << result.x.at(i);
                        out << YAML::EndMap;
                    }
                }
                out << YAML::Value << YAML::EndSeq;
            out << YAML::EndMap;
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "X/Y values" << YAML::EndMap;
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (std::size_t i = 0; i < result.x.size(); ++i) {
                    for (std::size_t j = 0; j < result.x.size(); ++j) {
                        out << YAML::BeginMap;
                        out << YAML::Key << "value";
                        out << YAML::Value << result.x.at(j);
                        out << YAML::EndMap;
                    }
                }
                out << YAML::Value << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;

        // dependent variables
        out << YAML::Key << "dependent_variables";
        out << YAML::Value << YAML::BeginSeq;
            //  value
            out << YAML::BeginMap;
                out << YAML::Key << "header";
                out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "2deltaNLL" << YAML::EndMap;
                AddQualifiers(out);
                out << YAML::Key << "values";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& zVec : result.z) {
                    for (const auto& z : zVec) {
                        out << YAML::BeginMap;
                        out << YAML::Key << "value";
                        if (TRExFitter::USEHEPDATAROUNDING) {
                            out << YAML::Value << Common::KeepSignificantDigits(z,3);
                        } else {
                            out << YAML::Value << z;
                        }
                        out << YAML::EndMap;
                    }
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndSeq;
    out << YAML::EndMap;

    // Write to the file
    Write(out, "HEPData LikelihoodScan2D", folder + "/HEPData/LikelihoodScan2D_"+suffix+".yaml");
}

bool YamlConverter::PlotContainerIsOK(const YamlConverter::PlotContainer& container) const {

    const std::size_t nSamples = container.samples.size();
    const std::size_t nBins = static_cast<std::size_t>(container.errors->GetN());
    if (nBins == 0) return false;
    if ((container.signalYields.size() + container.backgroundYields.size()) != nSamples) return false;
    for (const auto& yields : container.signalYields) {
        if (yields.size() != nBins) return false;
    }

    for (const auto& yields : container.backgroundYields) {
        if (yields.size() != nBins) return false;
    }

    if (!container.data.empty() && container.data.size() != nBins) return false;

    return true;
}

void YamlConverter::WriteHEPDataSubmission(const YamlConverter::SubmissionContainer& container,
                                           const std::vector<std::string>& pois) const {

    if (pois.empty()) {
        LOG(ERROR) << "Vector of POIs is empty\n";
        return;
    }

    gSystem->mkdir((container.folder + "/HEPData").c_str());

    std::ofstream file;
    file.open((container.folder + "/HEPData/submission.yaml"));
    if (!file.is_open() || !file.good()) {
        LOG(ERROR) << "Cannot open submission.yaml file\n";
        return;
    }

    LOG(INFO) << "Creating submission.yaml file in " << container.folder << "/HEPData/\n";

    AddCorrelation(file);
    file << "\n---\n";
    AddRanking(file, pois);

    AddImpact(file, pois);

    AddPlots(file, container);

    if (container.useTables) {
        AddTables(file);
        file << "\n---\n";
    }

    for (const auto& iunfold : container.unfoldingNames) {
        Add(file, "Unfolding_"+iunfold+".yaml", "Unfolded data of " + iunfold);
        file << "\n---\n";
    }

    for (const auto& iscan : container.useLikelihoodScan) {
        AddLikelihoodScan(file, iscan);
    }

    for (const auto& iscan : container.use2DLikelihoodScan) {
        Add2DLikelihoodScan(file, iscan);
    }

    file.close();
}

void YamlConverter::AppendSubmissionWithCovariances(const std::vector<std::string>& params,
                                                    const std::string& path) const {

    LOG(INFO) << "Updating the submission file with the covariance matrices\n";

    std::ofstream file(path, std::ios_base::app);
    if (!file.good() || !file.is_open()) {
        LOG(ERROR) << "Cannot open the file at " << path << "\n";
        return;
    }

    for (const auto& iparam : params) {
        Add(file, "Covariance_postfit_"+iparam+".yaml", iparam + " postfit covariance matrix");
        file << "\n---\n";
        Add(file, "Covariance_prefit_"+iparam+".yaml", iparam + " prefit covariance matrix");
        file << "\n---\n";
    }
    Add(file, "Covariance_statOnly.yaml", "Stat only covariance matrix");
    file << "\n---\n";
    Add(file, "Combined_Covariance_postfit.yaml", "Combined postfit covariance matrix");
    file << "\n---\n";
    Add(file, "Combined_Covariance_prefit.yaml", "Combined prefit covariance matrix");
    file << "\n---\n";
    Add(file, "CovarianceMatrix.yaml", "Covariance matrix obtained directly from the fit");
    file << "\n---\n";

    file.close();
}

void YamlConverter::AddCorrelation(std::ofstream& file) const {

    Add(file, "Correlation.yaml", "NP correlation matrix");
}

void YamlConverter::AddRanking(std::ofstream& file, const std::vector<std::string>& pois) const {
    for (const auto& ipoi : pois) {
        const std::string name = "Ranking_" + ipoi + ".yaml";
        Add(file, name, "NP ranking plot for " + ipoi);
        file << "\n---\n";
    }
}

void YamlConverter::AddImpact(std::ofstream& file, const std::vector<std::string>& pois) const {
    for (const auto& ipoi : pois) {
        const std::string name = "Impact_" + ipoi + ".yaml";
        Add(file, name, "Grouped impact table for " + ipoi);
        file << "\n---\n";
    }
}


void YamlConverter::AddPlots(std::ofstream& file,
                             const YamlConverter::SubmissionContainer& container) const {

    for (const auto& ireg : container.regionNames) {
        Add(file, ireg + "_prefit.yaml", ireg + " prefit");
        file << "\n---\n";
        Add(file, ireg + "_postfit.yaml", ireg + " postfit");
        file << "\n---\n";
    }
}

void YamlConverter::AddTables(std::ofstream& file) const {
    Add(file, "Table_prefit.yaml", "Prefit yields");
    file << "\n---\n";
    Add(file, "Table_postfit.yaml", "Postfit yields");
}

void YamlConverter::AddLikelihoodScan(std::ofstream& file, const std::string& scan) const {
    Add(file, "LikelihoodScan_"+scan+".yaml", "Likelihood scan for " + scan);
    file << "\n---\n";
}

void YamlConverter::Add2DLikelihoodScan(std::ofstream& file, const std::string& scan) const {
    Add(file, "LikelihoodScan2D_"+scan+".yaml", "Two-dimensional likelihood scan for " + scan);
    file << "\n---\n";
}

void YamlConverter::Add(std::ofstream& o, const std::string& file, const std::string& text) const {

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "data_file";
    out << YAML::Value << file;
    out << YAML::Key << "description";
    out << YAML::Value << "XXX";
    out << YAML::Key << "keywords";
    out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "name";
        out << YAML::Value << "reactions";
        out << YAML::Key << "values";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq << "XXX" << YAML::EndSeq;
        out << YAML::EndMap;
        out << YAML::BeginMap;
        out << YAML::Key << "name";
        out << YAML::Value << "phrases";
        out << YAML::Key << "values";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq << "XXX" << YAML::EndSeq;
        out << YAML::EndMap;
        out << YAML::BeginMap;
        out << YAML::Key << "name";
        out << YAML::Value << "cmenergies";
        out << YAML::Key << "values";
        out << YAML::Value << YAML::Flow << YAML::BeginSeq << "13000" << YAML::EndSeq;
        out << YAML::EndMap;
    out << YAML::EndSeq;
    out << YAML::Key << "location";
    out << YAML::Value << "XXX";
    out << YAML::Key << "name";
    out << YAML::Value << text;
    out << YAML::EndMap;

    o << out.c_str();
}

void YamlConverter::WriteCombiner(const std::string& name, const YamlConverter::Measurement& measurement, const std::string& path) const {
    if (measurement.pois.size() != measurement.values.size()) {
        LOG(ERROR) << "Sizes of POIs and their values do not match\n";
        return;
    }

    if (measurement.npCorrelation.GetNrows() != static_cast<int>(measurement.systematics.size())) {
        LOG(ERROR) << "Size of the systematics and the NPcorrelations do not match\n";
        return;
    }

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << name;
        out << YAML::Value << YAML::BeginMap;
            out << YAML::Key << "poiNames";
            out << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (const auto& ipoi : measurement.pois) {
                out << ipoi;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "values";
            out << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (const auto ivalue : measurement.values) {
                if (TRExFitter::USEHEPDATAROUNDING) {
                    out << Common::KeepSignificantDigits(ivalue,3);
                } else {
                    out << ivalue;
                }
            }
            out << YAML::EndSeq;
            if (measurement.statCov.GetNrows() > 0) {
                out << YAML::Key << "statCovMatrix";
                out << YAML::Value << YAML::BeginSeq;
                for (int i = 0; i < measurement.statCov.GetNrows(); ++i) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "row";
                    out << YAML::Value << YAML::Flow << YAML::BeginSeq;
                    for (int j = 0; j < i+1; ++j) {
                        if (TRExFitter::USEHEPDATAROUNDING) {
                            out << Common::KeepSignificantDigits(measurement.statCov(i,j),3);
                        } else {
                            out << measurement.statCov(i,j);
                        }
                    }
                    out << YAML::EndSeq;
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            }
            out << YAML::Key << "systematics";
                out << YAML::Value << YAML::BeginSeq;
                for (const auto& isyst : measurement.systematics) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "parameter";
                    out << YAML::Value << isyst.name;
                    out << YAML::Key << "impact";
                    out << YAML::Value << YAML::Flow << YAML::BeginSeq;
                    for (const auto iimpact : isyst.impact) {
                        if (TRExFitter::USEHEPDATAROUNDING) {
                            out << Common::KeepSignificantDigits(iimpact,3);
                        } else {
                            out << iimpact;
                        }
                    }
                    out << YAML::EndSeq;
                    out << YAML::Key << "constraint";
                    if (TRExFitter::USEHEPDATAROUNDING) {
                        out << YAML::Value << Common::KeepSignificantDigits(isyst.constraint,3);
                    } else {
                        out << YAML::Value << isyst.constraint;
                    }
                    out << YAML::Key << "pull";
                    if (TRExFitter::USEHEPDATAROUNDING) {
                        out << YAML::Value << Common::KeepSignificantDigits(isyst.pull,3);
                    } else {
                        out << YAML::Value << isyst.pull;
                    }
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::Key << "nuisanceParameters";
                out << YAML::Value << YAML::BeginMap;
                out << YAML::Key << "names";
                out << YAML::Value << YAML::Flow << YAML::BeginSeq;
                for (const auto& inp : measurement.systematics) {
                    out << inp.name;
                }
                out << YAML::EndSeq;
                out << YAML::Key << "correlationMatrix";
                out << YAML::Value << YAML::BeginSeq;
                for (int i = 0; i < measurement.npCorrelation.GetNrows(); ++i) {
                    out << YAML::BeginMap;
                    out << YAML::Key << "row";
                    out << YAML::Value << YAML::Flow << YAML::BeginSeq;
                    for (int j = 0; j < i+1; ++j) {
                        if (TRExFitter::USEHEPDATAROUNDING) {
                            out << Common::KeepSignificantDigits(measurement.npCorrelation(i,j),3);
                        } else {
                            out << measurement.npCorrelation(i,j);
                        }
                    }
                    out << YAML::EndSeq;
                    out << YAML::EndMap;
                }
                out << YAML::EndSeq;
            out << YAML::EndMap;
        out << YAML::EndMap;
    out << YAML::EndMap;
    Write(out, "Combiner config", path);
}


std::string YamlConverter::AddLatexDelimiters(const std::string& input) {
    bool hasLatex(false);

    if (input.find('\\') != std::string::npos) hasLatex = true;
    if (!hasLatex && input.find('{') != std::string::npos) hasLatex = true;
    if (!hasLatex && input.find('}') != std::string::npos) hasLatex = true;

    if (!hasLatex) return input;

    std::string result = "$";
    result += input;
    result += "$";

    return result;
}

void YamlConverter::WriteMigrationResponseHEPData(const TH2* matrix,
                                                  const std::string& path,
                                                  const bool isResponse,
                                                  const std::string& name,
                                                  const std::string& region,
                                                  const bool trueHorizontal) const {

    gSystem->mkdir((path+"/HEPData").c_str());

    std::string xLabel;
    std::string yLabel;
    if(trueHorizontal) {
        xLabel = "True-level";
        yLabel = "Detector-level";
    }
    else {
        xLabel = "Detector-level";
        yLabel = "True-level";
    }

    const int binsX = matrix->GetNbinsX();
    const int binsY = matrix->GetNbinsY();

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "dependent_variables";
    out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "header";
        out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << (isResponse ? "response matrix element" : "migration matrix element") << YAML::EndMap;
        AddQualifiers(out);
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (int ibin = 1; ibin <= binsX; ++ibin) {
            for (int jbin = 1; jbin <= binsY; ++jbin) {
                out << YAML::BeginMap;
                out << YAML::Key << "value";
                out << YAML::Value << Form("%2.4e", matrix->GetBinContent(ibin, jbin));
                out << YAML::EndMap;
            }
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    out << YAML::EndSeq;
    out << YAML::Key << "independent_variables";
    out << YAML::Value << YAML::BeginSeq;
        out << YAML::BeginMap;
        out << YAML::Key << "header";
        out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << xLabel << YAML::EndMap;
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (int ibin = 1; ibin <= binsX; ++ibin) {
            for (int jbin = 1; jbin <= binsY; ++jbin) {
                out << YAML::BeginMap;
                out << YAML::Key << "value";
                out << static_cast<std::string>(region + " bin " + ibin);
                out << YAML::EndMap;
            }
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
        out << YAML::BeginMap;
        out << YAML::Key << "header";
        out << YAML::Value << YAML::BeginMap << YAML::Key << "name" << YAML::Value << yLabel << YAML::EndMap;
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (int ibin = 1; ibin <= binsX; ++ibin) {
            for (int jbin = 1; jbin <= binsY; ++jbin) {
                out << YAML::BeginMap;
                out << YAML::Key << "value";
                out << static_cast<std::string>(region + " bin " + jbin);
                out << YAML::EndMap;
            }
        }
        out << YAML::EndSeq;
        out << YAML::EndMap;
    out << YAML::EndSeq;
    out << YAML::EndMap;

    const std::string nameFinal = (isResponse ? "response_matrix_" : "migration_matrix_") + name;
    Write(out, isResponse ? "HEPData response matrix" : "HEPData migration matrix", path+"/HEPData/"+nameFinal+".yaml");
}
