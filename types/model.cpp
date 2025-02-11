#include "model.h"
#include "../json_loader.h"

namespace model {

int total_size;
int design_size;
int noise_size;
double max_tspan;

// matrices, can be used to speed up likelihood if they are constant
Eigen::MatrixXd design_matrix;
Eigen::MatrixXd total_matrix;

element_t timing_model;

optional_element_t pl_red_noise;
optional_element_t pl_dm_noise;
optional_element_t efac;
optional_element_t equad;

parameter_t parse_parameter(const json_node_t& json_param)
{
    parameter_t param;
    param.load_from_json(json_param);
    return param;
}

element_t create_model_element(const string_t& element_name)
{
    if (element_name == "Power Law Red Noise") {
        std::cout << "loading power law red noise " << std::endl;
        return std::make_unique<pl_red_noise_t>();
    } else if (element_name == "EFAC") {
        std::cout << "loading efac " << std::endl;

        return std::make_unique<efac_t>();
    } else if (element_name == "EQUAD") {
        std::cout << "loading equad " << std::endl;

        return std::make_unique<equad_t>();
    } else if (element_name == "Timing Model") {
        std::cout << "loading timing model " << std::endl;

        return std::make_unique<timing_model_t>();
    }
    // Add other model elements as needed
    std::cout << "Unknown model element: " << element_name << std::endl;
    throw std::runtime_error("Unknown model element: " + element_name);
}

void load_model(const string_t& filename)
{
    std::optional<json_node_t> model = globals::config.get_optional_value<json_node_t>("model");

    if (model.has_value()) {
        bool have_timing_model = false;
        const auto& json_model = model.value();

        const auto& json_elements = json_model.get_array<json_node_t>("elements");
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

            auto element = create_model_element(element_name);

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
                throw std::runtime_error("Element " + element_name + " is not fully specified.");
            }

            if (element_name == "Power Law Red Noise") {
                model::pl_red_noise = std::move(element);
                model::pl_red_noise.value()->print();
            } else if (element_name == "Power Law DM Noise") {
                model::pl_dm_noise = std::move(element);
                model::pl_dm_noise.value()->print();
            } else if (element_name == "EFAC") {
                model::efac = std::move(element);
                model::efac.value()->print();
            } else if (element_name == "EQUAD") {
                model::equad = std::move(element);
                model::equad.value()->print();
            } else if (element_name == "Timing Model") {
                have_timing_model = true;
                model::timing_model = std::move(element);
                model::timing_model->print();
            }
        }

        // if we havn't loaded a timing model, create an empty one
        if (!have_timing_model) {
            auto element = create_model_element("Timing Model");
            model::timing_model = std::move(element);
        }
    }
}

int get_model_dims()
{
    int dims = 0;

    if (pl_red_noise.has_value()) {
        dims += pl_red_noise.value()->get_fitted_dims();
    }

    if (pl_dm_noise.has_value()) {
        dims += pl_dm_noise.value()->get_fitted_dims();
    }

    if (efac.has_value()) {
        dims += efac.value()->get_fitted_dims();
    }

    if (equad.has_value()) {
        dims += equad.value()->get_fitted_dims();
    }

    dims += timing_model->get_fitted_dims();

    return dims;
}
}  // namespace model