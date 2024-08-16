#include "efac.h"
#include <iostream>

void efac_t::set_parameter(const string_t& name, const parameter_t& param,
                           const rapidjson::Value& param_json)
{
    if (name == "global") {
        global = param;
    } else if (name == "per_flag") {

        flag = param_json["flag"].GetString();
        flag_indices = Eigen::VectorXi::Zero(globals::pulsar->nobs);
        flag_values.clear();

        for (int o = 0; o < globals::pulsar->nobs; o++) {
            bool found = false;
            for (int f = 0; f < globals::pulsar->obsn[o].nFlags; f++) {
                string_t obs_flag(globals::pulsar->obsn[o].flagID[f]);
                if (obs_flag == flag) {

                    string_t flag_value(globals::pulsar->obsn[o].flagVal[f]);
                    auto it = std::find(flag_values.begin(), flag_values.end(), flag_value);
                    if (it != flag_values.end()) {
                        auto index = std::distance(flag_values.begin(), it);
                        flag_indices(o) = index;
                    } else {

                        std::cout << "Found new EFAC " << flag << " "
                                  << globals::pulsar->obsn[o].flagVal[f] << std::endl;

                        flag_values.push_back(globals::pulsar->obsn[o].flagVal[f]);
                        flag_indices(o) = flag_values.size() - 1;
                    }
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw std::runtime_error("No flag found for EFAC on observation " +
                                         std::to_string(o));
            }
        }

        per_flag = param;

    } else {
        throw std::runtime_error("Invalid parameter name for EFAC: " + name);
    }
}

bool efac_t::is_fully_specified() const
{
    return true;
}

bool efac_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {"global", "per_flag"};
    return valid_params.find(param_name) != valid_params.end();
}

void efac_t::print() const
{
    std::cout << "EFAC Element:" << std::endl;
    if (global.has_value()) {
        std::cout << "global: ";
        global->print();
    }
    if (per_flag.has_value()) {
        std::cout << "per_flag: ";
        per_flag->print();
    }
}

int efac_t::get_fitted_dims()
{
    int fitted_dims = 0;
    if (global.has_value()) {
        fitted_dims += 1;
    }
    if (per_flag.has_value()) {
        fitted_dims += flag_values.size();
    }
    return fitted_dims;
}

string_t efac_t::get_name() const
{
    return "EFAC";
}

void efac_t::apply(double* Cube, Eigen::VectorXd& noise, double& prior_term, int& p_index)
{
    if (global.has_value()) {
        double multiplier = global.value().get_exp_value(Cube[p_index]);

        if (global.value().prior_type == prior_type_t::uniform) {
            prior_term += log(multiplier);
        }

        p_index++;

        noise = (noise * multiplier).array().square();
    }

    if (per_flag.has_value()) {
        Eigen::VectorXd multipliers = Eigen::VectorXd::Ones(flag_values.size());
        for (size_t i = 0; i < flag_values.size(); i++) {
            multipliers(i) = per_flag.value().get_exp_value(Cube[p_index++]);

            if (per_flag.value().prior_type == prior_type_t::uniform) {
                prior_term += log(multipliers(i));
            }
        }

        noise = noise.cwiseProduct(multipliers(flag_indices)).array().square();
    }
}