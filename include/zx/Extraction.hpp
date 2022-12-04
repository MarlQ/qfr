/*
* This file is part of MQT QFR library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/quantum/ for more information.
*/

#pragma once

#include "Definitions.hpp"
#include "Rational.hpp"
#include "Utils.hpp"
#include "ZXDiagram.hpp"
#include "dd/Control.hpp"
#include "zx/FunctionalityConstruction.hpp"
#include "QuantumComputation.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>



namespace zx {

    void extract(qc::QuantumComputation& circuit, zx::ZXDiagram& diag);
    void testExtraction();
    bool updateFrontier(zx::ZXDiagram& diag, std::vector<zx::Vertex>& frontier, std::vector<zx::Vertex>& outputs, qc::QuantumComputation& circuit);
    void processFrontier(zx::ZXDiagram& diag, std::vector<zx::Vertex>& frontier, std::vector<zx::Vertex>& outputs, qc::QuantumComputation& circuit);
    std::vector<std::pair<int, int>> permutation_as_swaps(std::map<int, int> perm);
}

