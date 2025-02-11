#include "json_node.h"
#include <fstream>
#include <iostream>
#include "json_loader.h"

bool json_node_t::has_member(const string_t& path) const
{
    return navigate_to_path(path) != nullptr;
}

void json_node_t::print(std::ostream& out) const
{
    if (!value_)
        return;
    print_value(*value_, out);
}

void json_node_t::print_value(const rapidjson::Value& value, std::ostream& out)
{
    rapidjson::OStreamWrapper osw(out);
    rapidjson::PrettyWriter<rapidjson::OStreamWrapper> writer(osw);
    value.Accept(writer);
    out << std::endl;
}

void json_node_t::mark_as_used(const string_t& path) const
{
    if (loader_) {
        loader_->mark_as_used(path);
    }
}

const rapidjson::Value* json_node_t::navigate_to_path(const string_t& inner_path) const
{
    if (!value_ || !value_->IsObject()) {
        return nullptr;
    }

    if (inner_path.empty()) {
        return value_;
    }

    const rapidjson::Value* current = value_;
    std::size_t start = 0;
    string_t current_path = path_;

    // Navigate through nested objects until our target
    while (true) {
        // Find next key in path
        std::size_t pos = inner_path.find('.', start);

        // When pos is npos, substr will automatically take until the end
        string_t key = (pos == string_t::npos) ? inner_path.substr(start) : inner_path.substr(start, pos - start);

        current_path = make_json_path(current_path, key);

        if (!current->IsObject()) {
            return nullptr;
        }

        if (!current->HasMember(key.c_str())) {
            return nullptr;
        }

        current = &(*current)[key.c_str()];

        // If this was the final key, return the value
        if (pos == string_t::npos) {
            return current;
        }

        start = pos + 1;
    }
}

std::optional<json_node_t> json_node_t::get_node_at_path(const string_t& inner_path) const
{
    string_t full_path = make_json_path(path_, inner_path);
    mark_as_used(full_path);

    const rapidjson::Value* found = navigate_to_path(inner_path);
    if (!found) {
        return std::nullopt;
    }

    return json_node_t(found, *loader_, full_path);
}
