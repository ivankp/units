#include "Units.hpp"

using namespace units::quantities;
using namespace units::literals;

Energy<double> KineticEnergy(Mass<double> m, Speed<double> v) {
    return 0.5 * m * (v * v);
}

int main() {
    return KineticEnergy(3*"kg"_u, 7*"m/s"_u) / "J"_u > 10.;
}
