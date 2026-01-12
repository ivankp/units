#pragma once

#include <cstdint>
#include <cmath>
#include <type_traits>

namespace units {

template <typename T, int8_t kg, int8_t m, int8_t s>
class Quantity;

template <typename T>
struct IsQuantity : std::false_type { };

template <typename T, int8_t kg, int8_t m, int8_t s>
struct IsQuantity<Quantity<T, kg, m, s>> : std::true_type { };

template <typename T>
concept AQuantity = IsQuantity<T>::value;

template <typename T>
concept NotAQuantity = !IsQuantity<T>::value;

template <typename T, int8_t kg, int8_t m, int8_t s>
class Quantity {
    T value;

    template <typename, int8_t, int8_t, int8_t>
    friend class Quantity;

public:
    constexpr Quantity() noexcept: value{} { }
    constexpr Quantity(T value) noexcept: value(value) { }

    static constexpr bool dimensionless = kg == 0 && m == 0 && s == 0;

    constexpr operator T() const noexcept requires dimensionless {
        return value;
    }

    // Comparison --------------------------------------------------------------

    template <typename R>
    constexpr auto operator<=>(Quantity<R, kg, m, s> r) const noexcept {
        return value <=> r.value;
    }

    template <typename R>
    constexpr bool operator==(Quantity<R, kg, m, s> r) const noexcept {
        return value == r.value;
    }

    // Addition ----------------------------------------------------------------

    template <typename R>
    friend constexpr auto operator+(Quantity l, Quantity<R, kg, m, s> r) noexcept
    -> Quantity<decltype(l.value + r.value), kg, m, s>
    {
        return l.value + r.value;
    }

    template <typename R>
    constexpr Quantity& operator+=(Quantity<R, kg, m, s> r) noexcept {
        value += r.value;
        return *this;
    }

    // Subtraction -------------------------------------------------------------

    template <typename R>
    friend constexpr auto operator-(Quantity l, Quantity<R, kg, m, s> r) noexcept
    -> Quantity<decltype(l.value - r.value), kg, m, s>
    {
        return l.value - r.value;
    }

    template <typename R>
    constexpr Quantity& operator-=(Quantity<R, kg, m, s> r) noexcept {
        value -= r.value;
        return *this;
    }

    // Multiplication ----------------------------------------------------------

    template <typename R, int8_t Rkg, int8_t Rm, int8_t Rs>
    friend constexpr auto operator*(Quantity l, Quantity<R, Rkg, Rm, Rs> r) noexcept
    -> Quantity<decltype(l.value * r.value), kg + Rkg, m + Rm, s + Rs>
    {
        return l.value * r.value;
    }

    friend constexpr auto operator*(Quantity l, NotAQuantity auto r) noexcept
    -> Quantity<decltype(l.value * r), kg, m, s>
    {
        return l.value * r;
    }

    friend constexpr auto operator*(NotAQuantity auto l, Quantity r) noexcept
    -> Quantity<decltype(l * r.value), kg, m, s>
    {
        return l * r.value;
    }

    constexpr Quantity& operator*=(auto r) noexcept {
        value *= r.value;
        return *this;
    }

    // Division ----------------------------------------------------------------

    template <typename R, int8_t Rkg, int8_t Rm, int8_t Rs>
    friend constexpr auto operator/(Quantity l, Quantity<R, Rkg, Rm, Rs> r) noexcept
    -> Quantity<decltype(l.value / r.value), kg - Rkg, m - Rm, s - Rs>
    {
        return l.value / r.value;
    }

    friend constexpr auto operator/(Quantity l, NotAQuantity auto r) noexcept
    -> Quantity<decltype(l.value / r), kg, m, s>
    {
        return l.value / r;
    }

    friend constexpr auto operator/(NotAQuantity auto l, Quantity r) noexcept
    -> Quantity<decltype(l / r.value), kg, m, s>
    {
        return l / r.value;
    }

    constexpr Quantity& operator/=(auto r) noexcept {
        value /= r.value;
        return *this;
    }

    // Exponentiation ----------------------------------------------------------

    template <int8_t n>
    friend constexpr Quantity<T, kg*n, m*n, s*n> pow(Quantity q) noexcept {
        return std::pow(q.value, n);
    }
};

template <typename T>
Quantity(T x) -> Quantity<T, 0, 0, 0>;

template <typename T>
using Dimensionless = Quantity<T, 0, 0, 0>;

template <typename T>
using Mass = Quantity<T, 1, 0, 0>;

template <typename T>
using Length = Quantity<T, 0, 1, 0>;

template <typename T>
using Time = Quantity<T, 0, 0, 1>;

template <typename T>
using Speed = Quantity<T, 0, 1, -1>;

template <typename T>
using Energy = Quantity<T, 1, 2, -2>;

} // namespace units
