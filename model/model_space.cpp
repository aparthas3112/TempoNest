#include "model_space.h"

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
}

void model_space_t::load_model()
{
    json_node_array_t json_elements = globals::config.get_array<json_node_t>("elements");

    bool have_timing_model = false;

    for (size_t i = 0; i < json_elements.size(); i++) {
        const auto& json_element = json_elements[i];
        string_t element_name = json_element.get_value<string_t>("name");

        // Check if any parameters are included
        bool any_param_included = false;
        const auto& json_params = json_element.get_array<json_node_t>("parameters");
        for (size_t j = 0; j < json_params.size(); j++) {
            if (json_params[j].get_value<bool>("include")) {
                any_param_included = true;
                break;
            }
        }

        if (!any_param_included) {
            std::cout << "Note: " << element_name << " element not included as all parameters are excluded." << std::endl;
            continue;
        }

        auto element = create_element(element_name);

        for (size_t j = 0; j < json_params.size(); j++) {
            const auto& json_param = json_params[j];
            string_t param_name = json_param.get_value<string_t>("name");
            if (element->is_valid_parameter(param_name)) {
                parameter_t param = parse_parameter(json_param);
                if (!param.include) {
                    continue;
                }
                element->set_parameter(param_name, param, json_param);
            } else {
                std::cout << "Warning: Ignoring invalid parameter '" << param_name << "' for " << element->get_name() << std::endl;
            }
        }

        if (!element->is_fully_specified()) {
            die("Element " + element_name + " is not fully specified.");
        }

        elements_[element_name] = std::move(element);
        elements_[element_name]->print();
        if (element_name == "Timing Model") {
            have_timing_model = true;
        }
    }

    // If we haven't loaded a timing model, create an empty one
    if (!have_timing_model) {
        auto element = create_element("Timing Model");
        elements_["Timing Model"] = std::move(element);
    }
}

parameter_t model_space_t::parse_parameter(const json_node_t& json_param)
{
    parameter_t param;
    param.load_from_json(json_param);
    return param;
}

element_t model_space_t::create_element(const string_t& type)
{
    if (type == "Power Law Red Noise") {
        return std::make_unique<pl_red_noise_t>();
    }
    if (type == "EFAC") {
        return std::make_unique<efac_t>();
    }
    if (type == "EQUAD") {
        return std::make_unique<equad_t>();
    }
    if (type == "Timing Model") {
        return std::make_unique<timing_model_t>();
    }

    // Add other element types here...

    die("Unknown element type: " + type);
}

int model_space_t::get_fitted_dims() const
{
    int total_dims = 0;
    for (const auto& [key, element] : elements_) {
        total_dims += element->get_fitted_dims();
    }
    return total_dims;
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

    noise_size_ = totCoeff;
    total_size_ = design_size_ + totCoeff;
}

void model_space_t::getEigenDVectorLike(Eigen::MatrixXd& design_matrix)
{

    auto& timing_model = get_element<timing_model_t>("Timing Model");

    int imargin = 0;

    FitInfo* fitinfo = &(globals::pulsar->fitinfo);

    for (int iparam = 0; iparam < fitinfo->nParams; ++iparam) {
        // check if this is something we want to marginalise over
        if (!timing_model.marginalised[iparam]) {
            continue;
        }

        param_label p = fitinfo->paramIndex[iparam];

        const int k = fitinfo->paramCounters[iparam];
        for (int iobs = 0; iobs < globals::pulsar->nobs; ++iobs) {
            const double x = globals::pulsar->obsn[iobs].bat - globals::pulsar->param[param_pepoch].val[0];

            design_matrix(iobs, imargin) = fitinfo->paramDerivs[iparam](globals::pulsar, 0, x, iobs, p, k);
        }
        ++imargin;
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

    delete[] DMVec;
    delete[] freqs;
}
