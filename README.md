# Overview
The Units library implements a type-based physical units system.

The goal is to provide easy-to-use means of units specification and conversion
with compile-time checks for correctness.

Requires C++20 or later standard.

## Benefits
* Type system checks units correctness.
  * Use of incompatible quantity is a compile-time error.
* Seamless conversion between unit systems.
  * Switching between units is no different syntactically than staying in one unit system.
* Terse syntax for unit literals.
* Self-documenting code.
  * Units have to be specified where quantities are used.
  * No implicit conversion between quantities with and without dimensions.
* Dimensionless quantities are implicitly convertible.
* All metric prefixes are supported.

## Integration
Copy one header file, [Units.hpp](include/Units.hpp), into your project.

# Examples
## Use in a function
Calculate kinetic energy for a body with a given mass and speed.
```c++
using namespace units::quantities;
using namespace units::literals;

Energy<double> KineticEnergy(Mass<double> m, Speed<double> v) {
    return 0.5 * m * Pow<2>(v);
}

// Later in the program
{
    // Only dimensionless quantities are implicitly convertible.
    // In this way, the type system ensures that units have been converted.

    // metric units
    double EnergyInJoules = KineticEnergy(3*"kg"_u, 7*"m/s"_u) / "J"_u;

    // arbitrary units
    double EnergyInkWh = KineticEnergy(3*"t"_u /* metric tonnes */, 7*"mi/h"_u) / "kWh"_u;
}
```

## Conversion factors
Determine a conversion factor.
```c++
using namespace units::literals;

constexpr double MPa_to_Pa = "MPa / Pa"_u; // same as 1e6
```
These literals can be used directly in place of constant conversion factors that are not type safe
with respect to units.

# Documentation
**Coming soon**

# Measurements
## Compilation timing
Measurements of compilation time
(see [`tests/compilation-timing`](tests/compilation-timing))
show that using this library added only about half as much
compilation time as simply including the `<iostream>` header.

The following table summarizes compilation timing measurement from 200
compilation runs with `g++ -std=c++20 -Wall -O3`.

| Test               | Mean Time, ms | Stdev, ms |
|--------------------|--------------:|----------:|
| **Without units**  |       43.150  |  1.701    |
|   With    units    |      155.905  |  1.222    |
|   Difference       |      112.755  |  2.094    |
|   Include iostream |      241.655  |  1.664    |
|   Difference       |      198.505  |  2.379    |

## Optimization
In a test with a complex dimensionless literal cast to a POD value
(see [`tests/assembly`](tests/assembly)),
both `g++` and `clang++` optimize the entire abstraction away at `-O1`,
and generate assembly with a single instruction to return an `int` value.
