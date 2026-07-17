// Original black-box conformance driver for ViennaLS v5.8.3.
//
// This file contains no Vienna implementation code. It uses the public ViennaLS
// API to evolve the same analytic 2-D keyhole cross-section used by petch's
// extruded 3-D manufactured audit. Build provenance is supplied by the Python
// adapter; event output is plain CSV.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsVelocityField.hpp>

namespace ls = viennals;

namespace {

constexpr int D = 2;
using Numeric = double;
using Domain = ls::Domain<Numeric, D>;

class ConstantNormalVelocity final : public ls::VelocityField<Numeric> {
  Numeric speed_;

public:
  explicit ConstantNormalVelocity(Numeric speed) : speed_(speed) {}

  Numeric getScalarVelocity(const ls::Vec3D<Numeric> &, int,
                            const ls::Vec3D<Numeric> &,
                            unsigned long) final {
    return speed_;
  }

  ls::Vec3D<Numeric> getVectorVelocity(const ls::Vec3D<Numeric> &, int,
                                       const ls::Vec3D<Numeric> &,
                                       unsigned long) final {
    return {0.0, 0.0, 0.0};
  }
};

ls::SmartPointer<Domain>
makeBox(const Numeric *bounds, Domain::BoundaryType *boundary,
        Numeric dx, Numeric x0, Numeric z0, Numeric x1, Numeric z1) {
  auto result = ls::SmartPointer<Domain>::New(bounds, boundary, dx);
  Numeric lower[D] = {x0, z0};
  Numeric upper[D] = {x1, z1};
  ls::MakeGeometry<Numeric, D>(
      result, ls::SmartPointer<ls::Box<Numeric, D>>::New(lower, upper))
      .apply();
  return result;
}

void unite(ls::SmartPointer<Domain> target, ls::SmartPointer<Domain> addition) {
  ls::BooleanOperation<Numeric, D>(
      target, addition, ls::BooleanOperationEnum::UNION)
      .apply();
}

ls::SmartPointer<Domain> buildKeyhole(Numeric dx, bool sealed = false) {
  const Numeric bounds[2 * D] = {0.0, 0.60, 0.0, 1.00};
  Domain::BoundaryType boundary[D] = {
      // The central keyhole is separated from both x boundaries by solid.
      // Reflective x is therefore physically equivalent for T1 and avoids
      // conflating ViennaLS periodic graph bookkeeping with normal motion.
      Domain::BoundaryType::REFLECTIVE_BOUNDARY,
      Domain::BoundaryType::INFINITE_BOUNDARY,
  };

  // Construct one solid half-space and subtract a rounded chamber joined to a
  // narrow open neck. The same analytic circle/box union is used by the petch
  // black-box comparator and avoids grid-sensitive rectangular shoulder voids.
  auto solid = ls::SmartPointer<Domain>::New(bounds, boundary, dx);
  Numeric planeOrigin[D] = {0.0, 0.75};
  Numeric planeNormal[D] = {0.0, 1.0};
  ls::MakeGeometry<Numeric, D>(
      solid,
      ls::SmartPointer<ls::Plane<Numeric, D>>::New(
          planeOrigin, planeNormal))
      .apply();
  auto gas = ls::SmartPointer<Domain>::New(bounds, boundary, dx);
  Numeric chamberCenter[D] = {0.30, 0.35};
  ls::MakeGeometry<Numeric, D>(
      gas,
      ls::SmartPointer<ls::Sphere<Numeric, D>>::New(chamberCenter, 0.18))
      .apply();
  // The independent reverse-motion fixture has a prescribed 0.10 um cap.
  // Reversing the first discretely closed state would compare different cap
  // thicknesses because event localization is solver- and grid-dependent.
  const Numeric neckTop = sealed ? 0.65 : 1.00 + dx;
  unite(gas, makeBox(bounds, boundary, dx, 0.25, 0.35, 0.35, neckTop));
  ls::BooleanOperation<Numeric, D>(
      solid, gas, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();
  return solid;
}

struct TopologyReceipt {
  std::size_t components;
  std::size_t voidPoints;
  std::size_t activePoints;
};

TopologyReceipt topology(ls::SmartPointer<Domain> levelSet) {
  ls::MarkVoidPoints<Numeric, D> marker(levelSet);
  marker.setVoidTopSurface(ls::VoidTopSurfaceEnum::LEX_HIGHEST);
  marker.setSaveComponentIds(true);
  marker.apply();
  const auto *voidMarkers = levelSet->getPointData().getScalarData(
      ls::MarkVoidPoints<Numeric, D>::voidPointLabel, true);
  std::size_t voidPoints = 0;
  if (voidMarkers != nullptr) {
    voidPoints = static_cast<std::size_t>(std::count_if(
        voidMarkers->begin(), voidMarkers->end(),
        [](Numeric value) { return value > 0.5; }));
  }
  const TopologyReceipt receipt = {
      marker.getNumberOfComponents(), voidPoints, levelSet->getNumberOfPoints()};
  // MarkVoidPoints stores diagnostic arrays on the sparse points. ViennaCore
  // v2.2.1's expansion path cannot safely translate those arrays after the
  // topology has changed, so the black-box adapter consumes and clears them
  // before the next advection. They are diagnostics, not process state.
  levelSet->getPointData().setNumberOfScalarData(0);
  levelSet->getPointData().setNumberOfVectorData(0);
  return receipt;
}

Numeric advance(ls::SmartPointer<Domain> levelSet, Numeric speed, Numeric dt,
                bool ignoreVoids) {
  auto velocity = ls::SmartPointer<ConstantNormalVelocity>::New(speed);
  ls::Advect<Numeric, D> kernel(levelSet, velocity);
  kernel.setAdvectionTime(dt);
  kernel.setSpatialScheme(ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
  kernel.setTemporalScheme(ls::TemporalSchemeEnum::FORWARD_EULER);
  kernel.setIgnoreVoids(ignoreVoids);
  kernel.apply();
  return kernel.getAdvectedTime();
}

void write(std::ofstream &stream, const std::string &phase, Numeric phaseTime,
           Numeric physicalTime, const TopologyReceipt &receipt) {
  stream << phase << ',' << std::setprecision(17) << phaseTime << ','
         << physicalTime << ',' << receipt.components << ','
         << receipt.voidPoints << ',' << receipt.activePoints << '\n';
  stream.flush();
}

Numeric positiveArgument(const char *text, const char *name) {
  const Numeric value = std::stod(text);
  if (!std::isfinite(value) || value <= 0.0)
    throw std::invalid_argument(std::string(name) + " must be positive");
  return value;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: " << argv[0]
              << " DX_UM MAX_COAT_S MAX_ETCH_S OUTPUT.csv\n";
    return 64;
  }
  try {
    const Numeric dx = positiveArgument(argv[1], "dx");
    const Numeric maximumCoat = positiveArgument(argv[2], "maximum coat time");
    const Numeric maximumEtch = positiveArgument(argv[3], "maximum etch time");
    const Numeric dt = 5.0 * dx;
    const Numeric coatSpeed = 0.025; // um / s
    const Numeric etchSpeed = -0.050; // um / s
    auto levelSet = buildKeyhole(dx);
    if (std::getenv("PETCH_VIENNALS_DEBUG_VTK") != nullptr) {
      auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
      ls::ToSurfaceMesh<Numeric, D>(levelSet, mesh).apply();
      ls::VTKWriter<Numeric>(mesh, std::string(argv[4]) + ".initial.vtp").apply();
    }
    std::ofstream output(argv[4]);
    if (!output)
      throw std::runtime_error("could not open output CSV");
    output << "phase,phase_time_s,physical_time_s,components,void_points,"
              "active_points\n";

    Numeric physicalTime = 0.0;
    Numeric phaseTime = 0.0;
    auto receipt = topology(levelSet);
    write(output, "initial", phaseTime, physicalTime, receipt);
    bool enclosed = receipt.voidPoints > 0;
    for (std::size_t step = 0;
         !enclosed && phaseTime + 0.5 * dt <= maximumCoat; ++step) {
      const Numeric accepted = advance(levelSet, coatSpeed, dt, false);
      phaseTime += accepted;
      physicalTime += accepted;
      receipt = topology(levelSet);
      write(output, "coat", phaseTime, physicalTime, receipt);
      enclosed = receipt.voidPoints > 0;
    }
    if (!enclosed) {
      std::cerr << "keyhole did not enclose before coat ceiling\n";
      return 2;
    }

    // Restart the reverse branch from the identical analytic sealed geometry,
    // rather than from this solver's grid-snapped first-closure state.
    levelSet = buildKeyhole(dx, true);
    physicalTime = 0.0;
    receipt = topology(levelSet);
    if (receipt.voidPoints == 0)
      throw std::runtime_error("analytic sealed keyhole contains no resolved void");
    write(output, "sealed_reference", 0.0, physicalTime, receipt);

    phaseTime = 0.0;
    bool reopened = false;
    for (std::size_t step = 0;
         !reopened && phaseTime + 0.5 * dt <= maximumEtch; ++step) {
      // T1 is a geometry-only reverse normal motion. Voids remain active, which
      // is ViennaLS's default and matches the petch T1 comparator mode.
      const Numeric accepted = advance(levelSet, etchSpeed, dt, false);
      phaseTime += accepted;
      physicalTime += accepted;
      receipt = topology(levelSet);
      write(output, "etch", phaseTime, physicalTime, receipt);
      reopened = receipt.voidPoints == 0;
    }
    if (!reopened) {
      std::cerr << "keyhole did not reopen before etch ceiling\n";
      return 3;
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "viennals keyhole probe: " << error.what() << '\n';
    return 70;
  }
}
