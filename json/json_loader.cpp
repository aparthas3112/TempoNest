#include "json_loader.h"
#include <unistd.h>
#include <fstream>
#include <iostream>
#include "../utils.h"
#include "logger.h"

/**
 * Creates a temporary file containing the given JSON content.
 * Uses mkstemp which requires a template filename ending in "XXXXXX".
 * These X's are replaced with random characters to ensure a unique filename.
 * Returns the generated filename which should be unlinked after use.
 */
std::string create_temp_json_file(const std::string& content)
{
    // Template must end in XXXXXX - these will be replaced with random characters
    char filename_template[] = "json_test_XXXXXX";
    int fd = mkstemp(filename_template);
    if (fd == -1) {
        throw std::runtime_error("Failed to create temporary file");
    }

    // Write the content to the file
    if (write(fd, content.c_str(), content.size()) != static_cast<ssize_t>(content.size())) {
        close(fd);
        unlink(filename_template);  // Delete the file
        throw std::runtime_error("Failed to write to temporary file");
    }

    close(fd);
    return std::string(filename_template);
}

bool load_temp_json_file(const std::string& content, json_loader_t& loader)
{
    std::string filename = create_temp_json_file(content);
    bool success = loader.load_json(filename);
    unlink(filename.c_str());
    return success;
}

bool json_loader_t::load_json(
    const string_t& filename,
    std::optional<std::reference_wrapper<std::unordered_set<string_t>>> included_files)
{
    // If no set provided, create an empty one
    std::unordered_set<string_t> empty_set;
    std::unordered_set<string_t>& local_included_files =
        included_files.has_value() ? included_files->get() : empty_set;

    if (!local_included_files.insert(filename).second) {
        die("File included multiple times: " + filename);
    }

    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        logger::log_error("Could not open file for reading: " + filename);
        return false;
    }

    rapidjson::IStreamWrapper isw(ifs);
    doc.ParseStream(isw);

    if (doc.HasParseError()) {
        logger::log_error("Error parsing JSON file: " + filename);
        return false;
    }

    // Clear used paths
    used_paths.clear();

    // Initialize included_files set and start processing
    string_t base_path = filename.substr(0, filename.find_last_of("/\\") + 1);
    process_includes(base_path, local_included_files);

    return true;
}

std::optional<json_node_t> json_loader_t::get_node_at_path(const string_t& path) const
{
    return get_root_node().get_node_at_path(path);
}

void json_loader_t::print(std::ostream& out) const
{
    json_node_t(&doc, *this, "").print(out);
}

bool json_loader_t::print_node(const string_t& path, std::ostream& out) const
{
    std::optional<json_node_t> node = get_node_at_path(path);
    if (!node.has_value()) {
        return false;
    }
    node->print(out);
    return true;
}

void json_loader_t::mark_as_used(const string_t& path) const
{
    used_paths.insert(path);
}

void json_loader_t::collect_unused_paths(const string_t& path_to_value,
                                         const rapidjson::Value& value,
                                         std::vector<string_t>& unused_paths) const
{
    if (value.IsObject()) {
        for (auto it = value.MemberBegin(); it != value.MemberEnd(); ++it) {
            string_t key = make_json_path(path_to_value, it->name.GetString());
            collect_unused_paths(key, it->value, unused_paths);
        }
    } else if (value.IsArray()) {
        // If this is an array of objects, check each element with its index
        if (value.Size() > 0 && value[0].IsObject()) {
            for (rapidjson::SizeType i = 0; i < value.Size(); i++) {
                string_t array_path = path_to_value + "[" + std::to_string(i) + "]";
                collect_unused_paths(array_path, value[i], unused_paths);
            }
        }
        // For arrays of primitives, just check the array path itself
        // this is because we will just get the array with get_array<T> and so from that point won't
        // need to interact with the json loader any further
        else {
            if (used_paths.find(path_to_value) == used_paths.end()) {
                unused_paths.push_back(path_to_value);
            }
        }
    } else {
        // This is a leaf value, check if its path was used
        if (used_paths.find(path_to_value) == used_paths.end()) {
            unused_paths.push_back(path_to_value);
        }
    }
}

bool json_loader_t::validate() const
{
    std::vector<string_t> unused_paths;
    collect_unused_paths("", doc, unused_paths);

    if (!unused_paths.empty()) {
        logger::log_error("Error: Found unused configuration entries:");
        for (const string_t& key : unused_paths) {
            logger::log_error("  " + key);
        }
        return false;
    }
    return true;
}

void json_loader_t::process_includes(const string_t& base_path,
                                     std::unordered_set<string_t>& included_files)
{
    // Only process objects - this check is important for recursive includes
    // Look for the INCLUDE directive at the root level

    if (!doc.IsObject() || !doc.HasMember("INCLUDE")) {
        return;
    }

    const auto& include = doc["INCLUDE"];

    // INCLUDE can be either a string (single file) or array of strings (multiple files)
    if (!include.IsString() && !include.IsArray()) {
        die("INCLUDE must be either a string or array of strings");
    }

    // Convert both single string and array cases to a vector for uniform processing
    std::vector<string_t> include_paths;
    if (include.IsString()) {
        include_paths.push_back(include.GetString());
    } else {
        // For array case, validate each element is a string
        for (const auto& inc : include.GetArray()) {
            if (!inc.IsString()) {
                die("All INCLUDE array elements must be strings");
            }
            include_paths.push_back(inc.GetString());
        }
    }

    // Remove the INCLUDE directive since it's a processing instruction
    // and shouldn't be part of the final merged document
    doc.RemoveMember("INCLUDE");

    // Process each include path in order
    for (const string_t& include_path : include_paths) {
        // Make path relative to the including file's location if it's not absolute
        string_t full_path = include_path;
        if (!full_path.empty() && full_path[0] != '/' && full_path[0] != '\\') {
            full_path = base_path + full_path;
        }

        json_loader_t temp_loader;
        if (!temp_loader.load_json(full_path, std::ref(included_files))) {
            die("Failed to load included file: " + full_path);
        }

        // Merge the fully processed included document into our document
        merge_document(temp_loader.doc);
    }
}

void json_loader_t::merge_document(const rapidjson::Document& other)
{
    // doc and other must be objects at this point so we call merge_objects
    merge_objects(doc, other, "");
}

void json_loader_t::merge_objects(rapidjson::Value& target, const rapidjson::Value& source,
                                  const string_t& path)
{
    if (!target.IsObject() || !source.IsObject()) {
        die("Can only merge object values at: " + path);
    }

    auto& alloc = doc.GetAllocator();

    // Iterate through all members of source
    for (auto it = source.MemberBegin(); it != source.MemberEnd(); ++it) {
        const char* key = it->name.GetString();
        const string_t current_path = make_json_path(path, key);

        // Check if the key exists in the target object
        if (!target.HasMember(key)) {
            // We need to create a new key string that's owned by our target document
            rapidjson::Value new_key(
                key, alloc);  // Creates copy of the key string using target's allocator

            // Create a new value object that will hold the copied value
            rapidjson::Value new_value;

            // Deep copy the value from source to our new value, using target's allocator
            new_value.CopyFrom(it->value, alloc);

            // Add the new key-value pair to the target object
            target.AddMember(new_key, new_value, alloc);

            continue;  // Move to next key-value pair
        }

        // If the key does exist in the target, we need to handle the merge recursively
        auto& target_value = target[key];

        // Handle different merge scenarios based on types
        if (target_value.IsArray() && it->value.IsArray()) {
            merge_arrays(target_value, it->value, current_path);
        } else if (target_value.IsObject() && it->value.IsObject()) {
            merge_objects(target_value, it->value, current_path);
        } else if (target_value.GetType() != it->value.GetType()) {
            die("Type mismatch at " + current_path);
        } else {
            die("Duplicate node is neither array nor object at " + current_path);
        }
    }
}

void check_array_types(const rapidjson::Value& arr, const string_t& arr_name,
                       rapidjson::Type& element_type, const string_t& path)
{
    for (const auto& elem : arr.GetArray()) {
        if (element_type == rapidjson::kNullType) {
            element_type = elem.GetType();
        } else if (elem.GetType() != element_type) {
            die("Mixed types in array at " + path + " in " + arr_name);
        }
    }
}

void json_loader_t::merge_arrays(rapidjson::Value& target, const rapidjson::Value& source,
                                 const string_t& path)
{

    // First verify these are both arrays
    if (!target.IsArray() || !source.IsArray()) {
        die("Can only merge array values at: " + path);
    }

    // Then verify all elements in both arrays are of the same type
    rapidjson::Type element_type = rapidjson::kNullType;

    check_array_types(target, "target", element_type, path);
    check_array_types(source, "source", element_type, path);

    auto& alloc = doc.GetAllocator();

    // Now merge the arrays
    for (const auto& elem : source.GetArray()) {
        rapidjson::Value new_elem;         // Create new container
        new_elem.CopyFrom(elem, alloc);    // Deep copy using target's allocator
        target.PushBack(new_elem, alloc);  // Add to target array
    }
}
