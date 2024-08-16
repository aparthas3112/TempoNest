#include "power_law_red_noise.h"
#include <iostream>
#include <stdexcept>

pl_red_noise_t::pl_red_noise_t() : num_freqs(33)
{
    frequencies = Eigen::VectorXd::Zero(num_freqs);
    for (int i = 0; i < frequencies.size(); i++) {
        frequencies[i] = i + 1;
    }
}

void pl_red_noise_t::set_parameter(const string_t& name, const parameter_t& param,
                                   const rapidjson::Value& param_json)
{
    if (name == "amplitude") {
        amplitude = param;
    } else if (name == "spectral_index") {
        spectral_index = param;
    } else if (name == "num_coeffs") {
        num_coeffs = param;
    } else {
        throw std::runtime_error("Invalid parameter name for Power Law Red Noise: " + name);
    }
}

bool pl_red_noise_t::is_valid_parameter(const string_t& param_name) const
{
    static const std::unordered_set<string_t> valid_params = {"amplitude", "spectral_index",
                                                              "num_coeffs"};
    return valid_params.find(param_name) != valid_params.end();
}

bool pl_red_noise_t::is_fully_specified() const
{
    if (amplitude.include == spectral_index.include) {
        return true;
    }
    if (num_coeffs.has_value() && num_coeffs.value().include && amplitude.include) {
        return true;
    }
    return false;
}

void pl_red_noise_t::print() const
{
    std::cout << "Power Law Red Noise Element:" << std::endl;
    std::cout << "Amplitude: ";
    amplitude.print();
    std::cout << "Spectral Index: ";
    spectral_index.print();
    if (num_coeffs.has_value()) {
        std::cout << "Cutoff: ";
        num_coeffs->print();
    }
}

int pl_red_noise_t::get_fitted_dims()
{
    return 2;
}

string_t pl_red_noise_t::get_name() const
{
    return "Power Law Red Noise";
}

void pl_red_noise_t::apply(double* Cube, Eigen::VectorXd& powercoeff, int& p_count, int& start_pos,
                           double maxtspan, double& uniform_prior, double& freq_det)
{
    double red_amp = amplitude.get_exp_value(Cube[p_count++]);
    double red_index = spectral_index.get_value(Cube[p_count++]);

    Eigen::VectorXd red_coeffs = (frequencies * 365.25 / maxtspan).array().pow(-red_index);

    double Tspan = maxtspan;
    double f1yr = 1.0 / 3.16e7;

    if (amplitude.prior_type == prior_type_t::uniform) {
        uniform_prior += log(red_amp);
    }

    double pl_amp =
        (red_amp * red_amp / 12.0 / (M_PI * M_PI)) * pow(f1yr, (-3)) / (Tspan * 24 * 60 * 60);

    red_coeffs *= pl_amp;

    powercoeff.segment(start_pos, num_freqs) += red_coeffs;
    powercoeff.segment(start_pos + num_freqs, num_freqs) += red_coeffs;

    freq_det += 2 * powercoeff.segment(start_pos, num_freqs).array().log().sum();
    start_pos += 2 * num_freqs;
}