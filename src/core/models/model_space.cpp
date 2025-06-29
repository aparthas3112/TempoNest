#include "model_space.h"
#include "../likelihood/gpu_functions.h"
#include "../utils/output_formatter.h"
#include "elements/timing_model.h"
#include "../utils/logger.h"
#include "t2fit.h"

model_space_t::model_space_t()
{
    load_model();

    // check that the model is fully specified
    if (!is_fully_specified()) {
        die("Model is not fully specified");
    }

    update_array_size_info();

    formBatsAll(globals::pulsar, globals::num_pulsars);
    formResiduals(globals::pulsar, globals::num_pulsars, 1);

    store_total_matrix();

    gpu_data::initialize(total_matrix_);
}

void model_space_t::load_model()
{
    json_node_array_t json_elements = globals::config.get_array<json_node_t>("elements");

    bool have_timing_model = false;
    int param_index = 0;
    int element_index = 1;

    // Print header for model elements section
    output_formatter::print_header("MODEL ELEMENTS");

    for (size_t i = 0; i < json_elements.size(); i++) {
        const auto& json_element = json_elements[i];
        string_t element_name = json_element.get_value<string_t>("name");

        // Check if any parameters are included or if marginalise flag is set
        bool any_param_included = false;
        const auto& json_params = json_element.get_array<json_node_t>("parameters");
        for (size_t j = 0; j < json_params.size(); j++) {
            if (json_params[j].get_value<bool>("include")) {
                any_param_included = true;
                break;
            }
        }

        // For timing model, check timing mode configuration
        bool should_marginalise = false;
        bool should_include_all_params = false;
        bool should_fit_all_params = false;
        double sigma_multiplier = 10.0;
        
        if (element_name == "Timing Model") {
            // Try to parse as string first (new modes: "all", "manual")
            auto marginalise_str_opt = json_element.get_optional_value<string_t>("marginalise");
            if (marginalise_str_opt.has_value()) {
                string_t marginalise_mode = marginalise_str_opt.value();
                if (marginalise_mode == "all") {
                    should_marginalise = true;
                    // Note will be shown in print_formatted
                } else if (marginalise_mode == "manual") {
                    // Note will be shown in print_formatted
                } else {
                    throw std::runtime_error("Invalid marginalise mode: " + marginalise_mode + ". Use 'all' or 'manual'");
                }
            } else {
                // Check for fit_all mode
                auto marginalise_bool_opt = json_element.get_optional_value<bool>("marginalise");
                auto fit_all_opt = json_element.get_optional_value<bool>("fit_all");
                
                if (marginalise_bool_opt.has_value() && !marginalise_bool_opt.value() && 
                    fit_all_opt.has_value() && fit_all_opt.value()) {
                    should_fit_all_params = true;
                    auto sigma_opt = json_element.get_optional_value<double>("sigma_multiplier");
                    if (sigma_opt.has_value()) {
                        sigma_multiplier = sigma_opt.value();
                    }
                    std::cout << "Note: Timing Model will fit all parameters with ±" << sigma_multiplier << "σ ranges" << std::endl;
                } else {
                    throw std::runtime_error("Invalid timing model configuration. Use marginalise='all'/'manual' or marginalise=false with fit_all=true");
                }
            }
        }

        if (!any_param_included && !should_marginalise && !should_include_all_params && !should_fit_all_params) {
            std::cout << "Note: " << element_name << " element not included as all parameters are excluded." << std::endl;
            continue;
        }

        auto element = create_element(element_name, json_element);
        auto* element_ptr = element.get();  // Get raw pointer before moving
        
        // Calculate frequencies for noise models before printing (requires temporary tspan calculation)
        if (element_name == "Power Law Red Noise" || element_name == "Power Law DM Noise") {
            // We need to calculate maxtspan early for proper display
            // This is a temporary calculation - will be done again later with final values
            double start = 1e20, end = -1e20;
            for (int i = 0; i < globals::pulsar->nobs; i++) {
                if ((double)globals::pulsar->obsn[i].bat < start) start = (double)globals::pulsar->obsn[i].bat;
                if ((double)globals::pulsar->obsn[i].bat > end) end = (double)globals::pulsar->obsn[i].bat;
            }
            double temp_maxtspan = 1 * (end - start);
            
            if (element_name == "Power Law Red Noise") {
                auto* red_noise = static_cast<pl_red_noise_t*>(element_ptr);
                red_noise->calculate_frequencies(temp_maxtspan);
            } else if (element_name == "Power Law DM Noise") {
                auto* dm_noise = static_cast<pl_dm_noise_t*>(element_ptr);
                dm_noise->calculate_frequencies(temp_maxtspan);
            }
        }

        // Auto-include all timing parameters for fit_all mode
        if (element_name == "Timing Model" && should_fit_all_params) {
            auto* timing_model = static_cast<timing_model_t*>(element_ptr);
            
            // Safety check: ensure timing model is properly initialized
            if (timing_model->t2_fitted_labels.empty()) {
                std::cout << "Warning: Timing model has no fitted parameters to auto-include" << std::endl;
            } else {
                for (size_t i = 0; i < timing_model->t2_fitted_labels.size(); i++) {
                    const string_t& param_name = timing_model->t2_fitted_labels[i];
                    
                    // Skip JUMP parameters as they're not supported for fitting
                    if (param_name.find("JUMP") != std::string::npos) {
                        continue;
                    }
                    
                    // Create parameter with configurable σ range
                    parameter_t param;
                    param.include = true;
                    param.prior_type = prior_type_t::uniform;
                    param.min_value = -sigma_multiplier;
                    param.max_value = sigma_multiplier;
                    param.name = param_name;
                    
                    // Store the current number of fitted dimensions before adding new parameter
                    int previous_dims = element->get_fitted_dims();
                    
                    // Set the initial index for this parameter
                    param.set_index(param_index);
                    param.set_parent(element_ptr);
                    
                    // Add the parameter to the timing model 
                    // For fit_all mode, we need at least one JSON parameter template
                    if (json_params.empty()) {
                        throw std::runtime_error("fit_all mode requires at least one parameter template in JSON configuration");
                    }
                    json_node_t dummy_json = json_params[0];
                    element->set_parameter(param_name, param, dummy_json);
                    
                    // Calculate how many new parameters were added
                    int new_dims = element->get_fitted_dims() - previous_dims;
                    param_index += new_dims;
                    
                    std::cout << "Auto-included timing parameter: " << param_name << " with ±" << sigma_multiplier << "σ range" << std::endl;
                }
            }
        }

        for (size_t j = 0; j < json_params.size(); j++) {
            const auto& json_param = json_params[j];
            string_t param_name = json_param.get_value<string_t>("name");
            if (element->is_valid_parameter(param_name)) {
                parameter_t param = parse_parameter(json_param);
                if (!param.include) {
                    continue;
                }
                // Store the current number of fitted dimensions before adding new parameter(s)
                int previous_dims = element->get_fitted_dims();

                // Set the initial index for this parameter
                param.set_index(param_index);
                param.set_parent(element_ptr);  // Set the parent before adding parameter

                // set_parameter might create multiple parameters internally
                element->set_parameter(param_name, param, json_param);
                // Calculate how many new parameters were added
                int new_dims = element->get_fitted_dims() - previous_dims;

                // Increment param_index by the number of new parameters
                param_index += new_dims;

            } else {
                std::cout << "Warning: Ignoring invalid parameter '" << param_name << "' for " << element->get_name() << std::endl;
            }
        }

        if (!element->is_fully_specified()) {
            die("Element " + element_name + " is not fully specified.");
        }

        elements_[element_name] = std::move(element);
        elements_[element_name]->print_formatted(element_index);
        element_index++;
        
        if (element_name == "Timing Model") {
            have_timing_model = true;
        }
    }

    // If we haven't loaded a timing model, create an empty one
    if (!have_timing_model) {
        auto element = create_element("Timing Model");
        elements_["Timing Model"] = std::move(element);
        elements_["Timing Model"]->print_formatted(element_index);
    }
    
    output_formatter::print_separator();
}

parameter_t model_space_t::parse_parameter(const json_node_t& json_param)
{
    parameter_t param;
    param.load_from_json(json_param);
    return param;
}

element_t model_space_t::create_element(const string_t& type, const std::optional<json_node_t>& element_json)
{
    if (type == "Power Law Red Noise") {
        return std::make_unique<pl_red_noise_t>(element_json);
    }
    if (type == "Power Law DM Noise") {
        return std::make_unique<pl_dm_noise_t>(element_json);
    }
    if (type == "EFAC") {
        return std::make_unique<efac_t>();
    }
    if (type == "EQUAD") {
        return std::make_unique<equad_t>();
    }
    if (type == "ECORR") {
        return std::make_unique<ecorr_t>();
    }
    if (type == "Timing Model") {
        return std::make_unique<timing_model_t>();
    }
    if (type == "Deterministic Solar Wind") {
        return std::make_unique<deterministic_solar_wind_t>();
    }
    if (type == "Stochastic Solar Wind") {
        return std::make_unique<stochastic_solar_wind_t>();
    }

    // Add other element types here...

    die("Unknown element type: " + type);
}

int model_space_t::get_fitted_dims() const
{
    return get_sampling_parameters().size();
}

bool model_space_t::is_fully_specified() const
{
    for (const auto& [name, element] : elements_) {
        if (!element->is_fully_specified()) {
            std::cout << "Element '" << name << "' is not fully specified" << std::endl;
            return false;
        }
    }
    return true;
}

std::vector<const parameter_t*> model_space_t::get_sampling_parameters() const
{
    std::vector<const parameter_t*> parameters;
    for (const auto& [name, element] : elements_) {
        auto element_params = element->get_parameters();
        parameters.insert(parameters.end(), element_params.begin(), element_params.end());
    }

    std::sort(parameters.begin(), parameters.end(), [](const parameter_t* a, const parameter_t* b) { return a->get_index() < b->get_index(); });

    return parameters;
}

void model_space_t::update_array_size_info()
{
    auto& timing_model = get_element<timing_model_t>("Timing Model");
    design_size_ = timing_model.design_size;

    // Calculate time span
    double start, end;
    int go = 0;
    for (int i = 0; i < globals::pulsar->nobs; i++) {
        if (globals::pulsar->obsn[i].deleted == 0) {
            if (go == 0) {
                go = 1;
                start = (double)globals::pulsar->obsn[i].bat;
                end = start;
            } else {
                if (start > (double)globals::pulsar->obsn[i].bat)
                    start = (double)globals::pulsar->obsn[i].bat;
                if (end < (double)globals::pulsar->obsn[i].bat)
                    end = (double)globals::pulsar->obsn[i].bat;
            }
        }
    }

    max_tspan_ = 1 * (end - start);

    // Calculate frequencies for noise models based on time span
    output_formatter::print_header("FREQUENCY COEFFICIENTS");
    std::cout << "    Time span: " << max_tspan_ << " days" << std::endl;
    
    if (auto pl_red_opt = get_optional_element<pl_red_noise_t>("Power Law Red Noise")) {
        auto& pl = pl_red_opt->get();
        pl.calculate_frequencies(max_tspan_);
        std::cout << "    Red Noise: " << pl.num_freqs << " frequencies (" << pl.days_per_coeff 
                  << " days/coeff, total coeffs: " << (2 * pl.num_freqs) << ")" << std::endl;
    }

    if (auto pl_dm_opt = get_optional_element<pl_dm_noise_t>("Power Law DM Noise")) {
        auto& pl = pl_dm_opt->get();
        pl.calculate_frequencies(max_tspan_);
        std::cout << "    DM Noise: " << pl.num_freqs << " frequencies (" << pl.days_per_coeff 
                  << " days/coeff, total coeffs: " << (2 * pl.num_freqs) << ")" << std::endl;
    }

    if (auto ecorr_opt = get_optional_element<ecorr_t>("ECORR")) {
        auto& ecorr = ecorr_opt->get();
        std::cout << "    ECORR: " << ecorr.get_num_coefficients() << " epochs" << std::endl;
    }

    // Calculate total coefficients
    int totCoeff = 0;

    if (auto pl_red_opt = get_optional_element<pl_red_noise_t>("Power Law Red Noise")) {
        auto& pl = pl_red_opt->get();
        totCoeff += 2 * pl.num_freqs;
    }

    if (auto pl_dm_opt = get_optional_element<pl_dm_noise_t>("Power Law DM Noise")) {
        auto& pl = pl_dm_opt->get();
        totCoeff += 2 * pl.num_freqs;
    }

    if (auto ecorr_opt = get_optional_element<ecorr_t>("ECORR")) {
        auto& ecorr = ecorr_opt->get();
        totCoeff += ecorr.get_num_coefficients();
    }

    noise_size_ = totCoeff;
    total_size_ = design_size_ + totCoeff;
}

void model_space_t::getEigenDVectorLike(Eigen::MatrixXd& design_matrix)
{
    auto& timing_model = get_element<timing_model_t>("Timing Model");

    int imargin = 0;

    // Use our safe parameter scanner data instead of unsafe fitinfo
    for (size_t iparam = 0; iparam < timing_model.t2_fitted_labels.size(); ++iparam) {
        // check if this is something we want to marginalise over
        if (!timing_model.marginalised[iparam]) {
            continue;
        }

        // Get parameter type and index from our validated data
        param_label p = timing_model.t2_fit_indices[iparam].first;
        const int k = timing_model.t2_fit_indices[iparam].second;

        // Get the parameter derivative function from Tempo2's fitting system
        // We need to find the corresponding fitinfo entry for the derivative function
        FitInfo* fitinfo = &(globals::pulsar->fitinfo);
        double (*param_deriv_func)(pulsar*, int, double, int, param_label, int) = nullptr;
        
        // Find matching parameter in fitinfo for derivative function
        for (int j = 0; j < fitinfo->nParams; ++j) {
            if (fitinfo->paramIndex[j] == p && fitinfo->paramCounters[j] == k) {
                param_deriv_func = fitinfo->paramDerivs[j];
                break;
            }
        }
        
        // If we found the derivative function, compute design matrix elements
        if (param_deriv_func != nullptr) {
            for (int iobs = 0; iobs < globals::pulsar->nobs; ++iobs) {
                const double x = globals::pulsar->obsn[iobs].bat - globals::pulsar->param[param_pepoch].val[0];
                design_matrix(iobs, imargin) = param_deriv_func(globals::pulsar, 0, x, iobs, p, k);
            }
            ++imargin;
        } else {
            // This shouldn't happen if our parameter detection is correct
            logger::log_error("No derivative function found for parameter: " + timing_model.t2_fitted_labels[iparam]);
            // Skip this parameter to avoid matrix bounds error
        }
    }
}

void model_space_t::store_total_matrix()
{
    total_matrix_ = Eigen::MatrixXd::Zero(globals::pulsar->nobs, total_size_);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Form the Design Matrix////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Forming Design Matrix " << design_size_ << std::endl;
    int TimetoMargin = design_size_;
    if (TimetoMargin > 0) {
        design_matrix_ = Eigen::MatrixXd::Zero(globals::pulsar->nobs, TimetoMargin);
        getEigenDVectorLike(design_matrix_);

        // Perform SVD
        Eigen::BDCSVD<Eigen::MatrixXd> svd(design_matrix_, Eigen::ComputeThinU | Eigen::ComputeThinV);

        Eigen::MatrixXd U = svd.matrixU();
        for (int i = 0; i < globals::pulsar->nobs; i++) {
            for (int j = 0; j < TimetoMargin; j++) {
                total_matrix_(i, j) = U(i, j);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////Set up Coefficients///////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    double* freqs = new double[noise_size_];
    double* DMVec = new double[globals::pulsar->nobs];

    double DMKappa = 2.410 * std::pow(10.0, -16);
    int startpos = 0;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////Red Noise///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (auto pl_red_opt = get_optional_element<pl_red_noise_t>("Power Law Red Noise")) {
        auto& pl = pl_red_opt->get();

        for (int i = 0; i < pl.num_freqs; i++) {
            freqs[startpos + i] = pl.frequencies[i] / max_tspan_;
            freqs[startpos + i + pl.num_freqs] = freqs[startpos + i];
        }

        for (int i = 0; i < pl.num_freqs; i++) {
            for (int k = 0; k < globals::pulsar->nobs; k++) {
                double time = (double)globals::pulsar->obsn[k].bat;
                total_matrix_(k, i + TimetoMargin + startpos) = cos(2 * M_PI * freqs[i] * time);
                total_matrix_(k, i + pl.num_freqs + TimetoMargin + startpos) = sin(2 * M_PI * freqs[i] * time);
            }
        }

        startpos += 2 * pl.num_freqs;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////DM Variations////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (auto pl_dm_opt = get_optional_element<pl_dm_noise_t>("Power Law DM Noise")) {
        auto& pl = pl_dm_opt->get();

        for (int o = 0; o < globals::pulsar->nobs; o++) {
            DMVec[o] = 1.0 / (DMKappa * std::pow((double)globals::pulsar->obsn[o].freqSSB, 2));
        }

        for (int i = 0; i < pl.num_freqs; i++) {
            freqs[startpos + i] = pl.frequencies[i] / max_tspan_;
            freqs[startpos + i + pl.num_freqs] = freqs[startpos + i];

            for (int k = 0; k < globals::pulsar->nobs; k++) {
                double time = (double)globals::pulsar->obsn[k].bat;
                total_matrix_(k, i + TimetoMargin + startpos) = cos(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];
                total_matrix_(k, i + pl.num_freqs + TimetoMargin + startpos) = sin(2 * M_PI * freqs[startpos + i] * time) * DMVec[k];
            }
        }

        startpos += 2 * pl.num_freqs;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////ECORR Epochs///////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (auto ecorr_opt = get_optional_element<ecorr_t>("ECORR")) {
        auto& ecorr = ecorr_opt->get();
        
        // Add ECORR quantization matrices to the design matrix
        ecorr.get_design_matrix(total_matrix_, TimetoMargin + startpos);
        
        startpos += ecorr.get_num_coefficients();
    }

    delete[] DMVec;
    delete[] freqs;
}
