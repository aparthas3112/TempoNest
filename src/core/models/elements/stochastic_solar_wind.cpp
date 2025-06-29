#include "stochastic_solar_wind.h"
#include "../../utils/settings.h"
#include <cmath>

void stochastic_solar_wind_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json)
{
    if (!is_valid_parameter(name)) {
        die("Invalid parameter name: " + name + " for " + get_name());
    }

    parameter_map_[name] = param;
}

bool stochastic_solar_wind_t::is_fully_specified() const
{
    // The element is fully specified if it has no parameters (not used) 
    // or if it has the log_amplitude parameter
    if (parameter_map_.empty()) {
        return true;  // Element not being used
    }
    
    return parameter_map_.find("log_amplitude") != parameter_map_.end();
}

bool stochastic_solar_wind_t::is_valid_parameter(const string_t& param_name) const
{
    return param_name == "log_amplitude";
}

string_t stochastic_solar_wind_t::get_name() const
{
    return "Stochastic Solar Wind";
}

void stochastic_solar_wind_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const
{
    if (parameter_map_.empty()) {
        return;
    }

    auto log_amplitude_param = get_optional_parameter("log_amplitude");
    if (log_amplitude_param) {
        int index = (*log_amplitude_param)->get_index();
        fprintf(par_file, "# Stochastic Solar Wind - Log Amplitude: %.6f +/- %.6f\n", 
                parameters[index], uncertainties[index]);
    }
}

void stochastic_solar_wind_t::apply(const std::vector<double>& parameter_values, Eigen::VectorXd& noise) const
{
    if (parameter_map_.empty()) {
        return;  // Element not being used
    }

    auto log_amplitude_param = get_optional_parameter("log_amplitude");
    if (!log_amplitude_param) {
        return;  // Parameter not found
    }

    int param_index = (*log_amplitude_param)->get_index();
    double log_amplitude = parameter_values[param_index];
    
    // Convert from log10 space: WhiteSolarWind = 10^(log_amplitude)
    double white_solar_wind = std::pow(10.0, log_amplitude);
    
    // Get reference solar wind electron density from pulsar data
    double ne_reference = globals::pulsar->ne_sw;
    
    // Apply stochastic solar wind noise to noise vector
    // Formula: SWTerm = WhiteSolarWind × tdis2[i] / ne_sw
    // noise[i] += SWTerm² (added to inverse noise variance)
    for (int i = 0; i < globals::pulsar->nobs; i++) {
        if (globals::pulsar->obsn[i].deleted == 0) {  // Only apply to non-deleted observations
            double tdis2 = globals::pulsar->obsn[i].tdis2;
            double sw_term = white_solar_wind * tdis2 / ne_reference;
            noise[i] += sw_term * sw_term;  // Add variance component
        }
    }
}