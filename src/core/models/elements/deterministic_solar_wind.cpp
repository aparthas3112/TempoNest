#include "deterministic_solar_wind.h"
#include "../../utils/settings.h"

void deterministic_solar_wind_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json)
{
    if (!is_valid_parameter(name)) {
        die("Invalid parameter name: " + name + " for " + get_name());
    }

    parameter_map_[name] = param;
}

bool deterministic_solar_wind_t::is_fully_specified() const
{
    // The element is fully specified if it has no parameters (not used) 
    // or if it has the electron_density parameter
    if (parameter_map_.empty()) {
        return true;  // Element not being used
    }
    
    return parameter_map_.find("electron_density") != parameter_map_.end();
}

bool deterministic_solar_wind_t::is_valid_parameter(const string_t& param_name) const
{
    return param_name == "electron_density";
}

string_t deterministic_solar_wind_t::get_name() const
{
    return "Deterministic Solar Wind";
}

void deterministic_solar_wind_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const
{
    if (parameter_map_.empty()) {
        return;
    }

    auto electron_density_param = get_optional_parameter("electron_density");
    if (electron_density_param) {
        int index = (*electron_density_param)->get_index();
        fprintf(par_file, "# Deterministic Solar Wind - Electron Density: %.6f +/- %.6f\n", 
                parameters[index], uncertainties[index]);
    }
}

void deterministic_solar_wind_t::apply(const std::vector<double>& parameter_values, Eigen::VectorXd& residuals) const
{
    if (parameter_map_.empty()) {
        return;  // Element not being used
    }

    auto electron_density_param = get_optional_parameter("electron_density");
    if (!electron_density_param) {
        return;  // Parameter not found
    }

    int param_index = (*electron_density_param)->get_index();
    double ne_fitted = parameter_values[param_index];
    
    // Get reference solar wind electron density from pulsar data
    double ne_reference = globals::pulsar->ne_sw;
    
    // Apply deterministic correction to residuals
    // Formula: residuals[i] -= (ne_fitted - ne_reference) × tdis2[i]
    for (int i = 0; i < globals::pulsar->nobs; i++) {
        if (globals::pulsar->obsn[i].deleted == 0) {  // Only apply to non-deleted observations
            double tdis2 = globals::pulsar->obsn[i].tdis2;
            residuals[i] -= (ne_fitted - ne_reference) * tdis2;
        }
    }
}