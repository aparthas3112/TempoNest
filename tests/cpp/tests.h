#pragma once
#include <memory>
#include "../../src/core/models/model.h"

void run_tests(const std::shared_ptr<model_t> model);
bool run_likelihood_tests(const std::shared_ptr<model_t> model);