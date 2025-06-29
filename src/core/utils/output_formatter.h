#pragma once

#include <iostream>
#include <string>
#include <iomanip>

namespace output_formatter {

    // Box drawing characters for clean sections
    const std::string BOX_TOP_LEFT = "┌";
    const std::string BOX_TOP_RIGHT = "┐";
    const std::string BOX_BOTTOM_LEFT = "└";
    const std::string BOX_BOTTOM_RIGHT = "┘";
    const std::string BOX_HORIZONTAL = "─";
    const std::string BOX_VERTICAL = "│";
    const std::string TREE_BRANCH = "├─";
    const std::string TREE_LAST = "└─";
    const std::string TREE_VERTICAL = "│";

    // Colors (optional - can be disabled on systems that don't support them)
    const std::string RESET = "\033[0m";
    const std::string BOLD = "\033[1m";
    const std::string GREEN = "\033[32m";
    const std::string YELLOW = "\033[33m";
    const std::string BLUE = "\033[34m";
    const std::string RED = "\033[31m";

    inline void print_header(const std::string& title, int width = 70) {
        std::cout << BOX_TOP_LEFT;
        for (int i = 0; i < width - 2; i++) std::cout << BOX_HORIZONTAL;
        std::cout << BOX_TOP_RIGHT << std::endl;
        
        int padding = (width - title.length() - 2) / 2;
        std::cout << BOX_VERTICAL << std::string(padding, ' ') 
                  << BOLD << title << RESET 
                  << std::string(width - title.length() - padding - 2, ' ') 
                  << BOX_VERTICAL << std::endl;
        
        std::cout << BOX_BOTTOM_LEFT;
        for (int i = 0; i < width - 2; i++) std::cout << BOX_HORIZONTAL;
        std::cout << BOX_BOTTOM_RIGHT << std::endl << std::endl;
    }

    inline void print_section(const std::string& title) {
        std::cout << BOLD << title << ":" << RESET << std::endl;
    }

    inline void print_setting(const std::string& name, const std::string& value, bool enabled = true) {
        std::string status_color = enabled ? GREEN : RED;
        std::cout << "  • " << name << ": " << status_color << value << RESET << std::endl;
    }

    inline void print_element_header(int index, const std::string& name, const std::string& extra = "") {
        std::cout << std::endl << "[" << index << "] " << BOLD << name << RESET;
        if (!extra.empty()) {
            std::cout << " " << YELLOW << "(" << extra << ")" << RESET;
        }
        std::cout << std::endl;
    }

    inline void print_parameter(const std::string& name, const std::string& prior_type, 
                               double min_val, double max_val, bool included, bool last = false) {
        std::string branch = last ? TREE_LAST : TREE_BRANCH;
        std::string status = included ? (GREEN + " ✓ included" + RESET) : (RED + " ✗ excluded" + RESET);
        
        std::cout << "    " << branch << " " << BOLD << name << RESET << ": " 
                  << prior_type << " [" << min_val << ", " << max_val << "]" << status << std::endl;
    }

    inline void print_marginalised_note() {
        std::cout << "    " << TREE_LAST << " " << YELLOW << "All fitted parameters from .par file will be marginalised" << RESET << std::endl;
    }

    inline void print_separator() {
        std::cout << std::endl;
    }

} // namespace output_formatter