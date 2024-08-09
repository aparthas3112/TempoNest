#include "json_loader.h"
#include <fstream>
#include <stdexcept>
#include "rapidjson/istreamwrapper.h"
#include "types/model.h"
#include "types/model_element.h"
#include "types/parameter.h"

parameter_t json_loader::parse_parameter(const rapidjson::Value& json_param)
{
    parameter_t param;
    param.name = json_param["name"].GetString();
    param.description = json_param["description"].GetString();

    param.prior_type = string_t(json_param["prior_type"].GetString()) == "uniform"
                           ? prior_type_t::uniform
                           : prior_type_t::log_uniform;

    param.include = json_param["include"].GetBool();
    param.fit = json_param["fit"].GetBool();
    param.min_value = json_param["min_value"].GetDouble();
    param.max_value = json_param["max_value"].GetDouble();
    return param;
}

element_t json_loader::create_model_element(const string_t& element_name)
{
    if (element_name == "Power Law Red Noise") {
        std::cout << "loading power law red noise " << std::endl;
        return std::make_unique<pl_red_noise_element>();
    } else if (element_name == "EFAC") {
        std::cout << "loading efac " << std::endl;

        return std::make_unique<efac_element>();
    } else if (element_name == "EQUAD") {
        std::cout << "loading equad " << std::endl;

        return std::make_unique<equad_element>();
    }
    // Add other model elements as needed
    throw std::runtime_error("Unknown model element: " + element_name);
}

rapidjson::Document json_loader::parse_json_file(const string_t& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    rapidjson::IStreamWrapper isw(file);
    rapidjson::Document doc;

    // Enable comments and trailing commas
    doc.ParseStream<rapidjson::kParseCommentsFlag | rapidjson::kParseTrailingCommasFlag>(isw);

    if (doc.HasParseError()) {
        throw std::runtime_error("JSON parse error");
    }

    return doc;
}

void json_loader::load_from_json(const string_t& filename)
{
    rapidjson::Document doc = parse_json_file(filename);

    // Load global settings
    if (doc.HasMember("global_settings")) {
        // Parse global settings as needed
    }

    // Load model elements
    if (doc.HasMember("elements")) {
        const auto& json_elements = doc["elements"];
        for (rapidjson::SizeType i = 0; i < json_elements.Size(); i++) {
            const auto& json_element = json_elements[i];
            string_t element_name = json_element["name"].GetString();

            // Check if any parameters are included
            bool any_param_included = false;
            const auto& json_params = json_element["parameters"];
            for (rapidjson::SizeType j = 0; j < json_params.Size(); j++) {
                if (json_params[j]["include"].GetBool()) {
                    any_param_included = true;
                    break;
                }
            }

            if (!any_param_included) {
                std::cout << "Note: " << element_name
                          << " element not included as all parameters are excluded." << std::endl;

                continue;
            }

            auto element = create_model_element(element_name);

            for (rapidjson::SizeType j = 0; j < json_params.Size(); j++) {
                const auto& json_param = json_params[j];
                string_t param_name = json_param["name"].GetString();
                if (element->is_valid_parameter(param_name)) {
                    element->set_parameter(param_name, parse_parameter(json_param));
                } else {
                    std::cout << "Warning: Ignoring invalid parameter '" << param_name << "' for "
                              << element->get_name() << std::endl;
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
            }
        }
    }
}