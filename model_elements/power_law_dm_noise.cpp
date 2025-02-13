#include "power_law_dm_noise.h"
#include <iostream>
#include <stdexcept>

pl_dm_noise_t::pl_dm_noise_t() : num_freqs(33)
{
    frequencies = Eigen::VectorXd::Zero(num_freqs);
    for (int i = 0; i < frequencies.size(); i++) {
        frequencies[i] = i + 1;
    }
}

void pl_dm_noise_t::set_parameter(const string_t& name, const parameter_t& param, const json_node_t& param_json)
{
    if (name == "amplitude") {
        amplitude = param;
        parameters_.push_back(&amplitude);
    } else if (name == "spectral_index") {
        spectral_index = param;
        parameters_.push_back(&spectral_index);
    } else if (name == "num_coeffs") {
        num_coeffs = param;
        parameters_.push_back(&num_coeffs.value());
    } else {
        throw std::runtime_error("Invalid parameter name for Power Law DM Noise: " + name);
    }
}

bool pl_dm_noise_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {"amplitude", "spectral_index", "num_coeffs"};
    return valid_params.find(param_name) != valid_params.end();
}

bool pl_dm_noise_t::is_fully_specified() const
{
    if (amplitude.include == spectral_index.include) {
        return true;
    }
    if (num_coeffs.has_value() && num_coeffs.value().include && amplitude.include) {
        return true;
    }
    return false;
}

void pl_dm_noise_t::print() const
{
    std::cout << "Power Law DM Noise Element:" << std::endl;
    std::cout << "Amplitude: ";
    amplitude.print();
    std::cout << "Spectral Index: ";
    spectral_index.print();
    if (num_coeffs.has_value()) {
        std::cout << "Cutoff: ";
        num_coeffs->print();
    }
}

void pl_dm_noise_t::write_to_par_file(FILE* par_file, const std::vector<double>& parameters, const std::vector<double>& uncertainties) const {}

int pl_dm_noise_t::get_fitted_dims()
{
    return 2;
}

string_t pl_dm_noise_t::get_name() const
{
    return "Power Law DM Noise";
}

void pl_dm_noise_t::apply(const std::vector<double>& parameter_values, Eigen::VectorXd& powercoeff, int& start_pos, double maxtspan, double& uniform_prior, double& freq_det) const
{

    double dm_amp = amplitude.get_exp_value(parameter_values);
    double dm_index = spectral_index.get_value(parameter_values);

    Eigen::VectorXd dm_coeffs = (frequencies * 365.25 / maxtspan).array().pow(-dm_index);

    if (amplitude.prior_type == prior_type_t::uniform) {
        uniform_prior += log(dm_amp);
    }

    double f1yr = 1.0 / 3.16e7;

    double pl_amp = (dm_amp * dm_amp) * pow(f1yr, (-3)) / (maxtspan * 24 * 60 * 60);

    dm_coeffs *= pl_amp;

    powercoeff.segment(start_pos, num_freqs) += dm_coeffs;
    powercoeff.segment(start_pos + num_freqs, num_freqs) += dm_coeffs;

    freq_det += 2 * powercoeff.segment(start_pos, num_freqs).array().log().sum();
    start_pos += 2 * num_freqs;
}