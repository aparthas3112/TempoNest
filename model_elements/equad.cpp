#include "equad.h"
#include <iostream>

void equad_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json)
{
    std::cout << "Setting EQUAD parameter: " << name << std::endl;

    // json_loader::print_node(param_json);
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

                        std::cout << "Found new EQUAD " << flag << " " << globals::pulsar->obsn[o].flagVal[f] << std::endl;

                        flag_values.push_back(globals::pulsar->obsn[o].flagVal[f]);
                        flag_indices(o) = flag_values.size() - 1;
                    }
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw std::runtime_error("No flag found for EQUAD on observation " + std::to_string(o));
            }
        }

        // Now create a parameter for each unique flag value
        for (size_t i = 0; i < flag_values.size(); i++) {
            string_t param_name = "per_flag::" + flag + "::" + std::to_string(i);
            parameter_t new_param = param;
            new_param.set_index(param.get_index() + i);
            parameter_map_[param_name] = new_param;
        }

    } else {
        throw std::runtime_error("Invalid parameter name for EFAC: " + name);
    }
}

bool equad_t::is_fully_specified() const
{
    if (get_optional_parameter("global").has_value() == flag_values.empty()) {
        return true;
    }
    return false;
}

bool equad_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {"global", "per_flag"};
    return valid_params.find(param_name) != valid_params.end();
}

void equad_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const {}

string_t equad_t::get_name() const
{
    return "EQUAD";
}

void equad_t::apply(const std::vector<double>& parameter_values, Eigen::VectorXd& noise, double& prior_term) const
{
    if (auto param = get_optional_parameter("global")) {
        double value = param.value()->get_exp_value(parameter_values);
        double equad = value * value;

        if (param.value()->prior_type == prior_type_t::uniform) {
            prior_term += log(value);
        }

        noise = noise.array() + equad;
    }

    if (!flag_values.empty()) {
        Eigen::VectorXd equad_values = Eigen::VectorXd::Zero(flag_values.size());
        for (size_t i = 0; i < flag_values.size(); i++) {
            string_t param_name = "per_flag::" + flag + "::" + std::to_string(i);
            auto param = get_parameter(param_name);
            double value = param->get_exp_value(parameter_values);
            equad_values(i) = value * value;

            if (param->prior_type == prior_type_t::uniform) {
                prior_term += log(value);
            }
        }

        noise += equad_values(flag_indices);
    }
}