#include "tests.h"
#include <iostream>

const char* SUCCESS_ART = R"(
   _____  _    _  _____ _____ ______  _____ _____
  / ____|| |  | |/ ____/ ____|  ____|/ ____/ ____|
 | (___  | |  | | |   | |    | |__  | (___| (___
  \___ \ | |  | | |   | |    |  __|  \___ \\___ \
  ____) || |__| | |___| |____| |____ ____) |___) |
 |_____/  \____/ \_____\_____|______|_____/_____/

)";

const char* FAILED_ART = R"(
  ______      _____ _      ______ _____  
 |  ____/\   |_   _| |    |  ____|  __ \ 
 | |__ /  \    | | | |    | |__  | |  | |
 |  __/ /\ \   | | | |    |  __| | |  | |
 | | / ____ \ _| |_| |____| |____| |__| |
 |_|/_/    \_\_____|______|______|_____/ 
                                         
)";

// Function to print in color
void print_color(const char* message, const char* color)
{
    printf("%s%s%s", color, message, "\033[0m");
}

void run_tests()
{
    bool success = true;

    success &= run_likelihood_tests();

    if (success) {
        print_color(SUCCESS_ART, "\033[0;32m");  // Green color
    } else {
        print_color(FAILED_ART, "\033[0;31m");  // Red color
    }
}