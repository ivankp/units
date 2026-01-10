#pragma once

#include <cstdint>
#include <cmath>

namespace units {

template <int8_t kg, int8_t m, int8_t s, typename T>
class Quantity {
    T value;

    template <int8_t, int8_t, int8_t, typename>
    friend class Quantity;

public:
    constexpr Quantity() noexcept: value{} { }
    constexpr Quantity(T value) noexcept: value(value) { }

    static constexpr bool dimensionless = kg == 0 && m == 0 && s == 0;

    constexpr operator T() const noexcept requires dimensionless {
        return value;
    }

    constexpr Quantity operator+(Quantity o) const noexcept {
        return value + o.value;
    }

    constexpr Quantity operator-(Quantity o) const noexcept {
        return value - o.value;
    }

    template <int8_t Rkg, int8_t Rm, int8_t Rs, typename RT>
    constexpr Quantity<kg + Rkg, m + Rm, s + Rs, decltype(T{} * RT{})>
    operator*(Quantity<Rkg, Rm, Rs, RT> o) const noexcept {
        return value * o.value;
    }

    template <int8_t Rkg, int8_t Rm, int8_t Rs, typename RT>
    constexpr Quantity<kg - Rkg, m - Rm, s - Rs, decltype(T{} / RT{})>
    operator/(Quantity<Rkg, Rm, Rs, RT> o) const noexcept {
        return value / o.value;
    }

    template <int8_t n>
    friend constexpr Quantity<kg*n, m*n, s*n, T> pow(Quantity q) noexcept {
        return std::pow(q.value, n);
    }
};

template <typename T>
Quantity(T x) -> Quantity<0, 0, 0, T>;

} // namespace units
