#include "efac.h"
#include <iostream>

void efac_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json)
{
    if (name == "global") {
        parameter_map_[name] = param;
    } else if (name == "per_flag") {

        flag = param_json.get_value<string_t>("flag");
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


                        flag_values.push_back(globals::pulsar->obsn[o].flagVal[f]);
                        flag_indices(o) = flag_values.size() - 1;
                    }
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw std::runtime_error("No flag found for EFAC on observation " + std::to_string(o));
            }
        }

        // Now create a parameter for each unique flag value
        for (size_t i = 0; i < flag_values.size(); i++) {
            string_t param_name = "per_flag::" + flag + "::" + flag_values[i];
            parameter_t new_param = param;
            new_param.set_index(param.get_index() + i);
            parameter_map_[param_name] = new_param;
        }
    } else {
        throw std::runtime_error("Invalid parameter name for EFAC: " + name);
    }
}

bool efac_t::is_fully_specified() const
{
    if (get_optional_parameter("global").has_value() == flag_values.empty()) {
        return true;
    }
    return false;
}

bool efac_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {"global", "per_flag"};
    return valid_params.find(param_name) != valid_params.end();
}

void efac_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const
{
    if (auto param = get_optional_parameter("global")) {
        fprintf(par_file, "TNGLobalEF %g\n", parameters[param.value()->get_index()]);
    }
    if (auto param = get_optional_parameter("per_flag")) {
        for (size_t f = 0; f < flag_values.size(); f++) {
            fprintf(par_file, "TNEF %s %s %g\n", flag.c_str(), flag_values[f].c_str(), parameters[param.value()->get_index()]);
        }
    }
}

string_t efac_t::get_name() const
{
    return "EFAC";
}

void efac_t::apply(const std::vector<double>& parameter_values, Eigen::VectorXd& noise, double& prior_term) const
{
    if (auto param = get_optional_parameter("global")) {
        double efac = param.value()->get_exp_value(parameter_values);

        if (param.value()->prior_type == prior_type_t::uniform) {
            prior_term += log(efac);
        }

        // Apply EFAC^2 to the error variance (noise is sigma, not sigma^2)
        // Legacy: 1/(EFAC^2 * sigma^2 + EQUAD + ...) 
        // So we multiply noise by EFAC: noise = sigma * EFAC
        noise = noise * efac;
    }

    if (!flag_values.empty()) {
        Eigen::VectorXd efacs = Eigen::VectorXd::Ones(flag_values.size());
        for (size_t i = 0; i < flag_values.size(); i++) {
            string_t param_name = "per_flag::" + flag + "::" + flag_values[i];
            auto param = get_parameter(param_name);
            efacs(i) = param->get_exp_value(parameter_values);

            if (param->prior_type == prior_type_t::uniform) {
                prior_term += log(efacs(i));
            }
        }

        // Apply EFAC to each observation based on its flag
        // Legacy: 1/(EFAC^2 * sigma^2 + EQUAD + ...)
        // So we multiply noise by EFAC: noise = sigma * EFAC
        for (int o = 0; o < noise.size(); o++) {
            noise(o) = noise(o) * efacs(flag_indices(o));
        }
    }
}