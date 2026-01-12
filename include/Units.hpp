#pragma once

#include <cstdint>
#include <cmath>
#include <type_traits>

namespace units {

namespace detail {

template <unsigned N>
struct StringLiteral {
    char s[N]{};

    constexpr StringLiteral(const char(&arr)[N]) noexcept {
        for (unsigned i = 0; i < N; ++i)
            s[i] = arr[i];
    }
};

template <typename T>
struct ParsedLiteral {
    T factor = 1;
    int8_t kg = 0, m = 0, s = 0;

    template <unsigned N>
    constexpr ParsedLiteral(const StringLiteral<N>& str) {
        for (unsigned i = 0; i < N; ++i) {
            const char c = str.s[i];
            if ((c < 'A' || 'Z' < c) && (c < 'a' || 'z' < c)) {
                // TODO: write a custom comparison function (maybe lambda)
                const std::string_view v(str.s, i);
                if (v == "kg") {
                    kg += 1;
                } else if (v == "m") {
                    m += 1;
                } else if (v == "s") {
                    s += 1;
                }
                break; // TODO: complete parsing algorithm
            }
        }
    }
};

} // namespace detail

// TODO: Vector class to represent dimensions

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

namespace literals {

template <detail::StringLiteral s>
constexpr auto operator ""_u() noexcept {
    constexpr detail::ParsedLiteral<double> p(s);
    return Quantity<double, p.kg, p.m, p.s>(p.factor);
}

template <detail::StringLiteral s>
constexpr auto operator ""_uf() noexcept {
    constexpr detail::ParsedLiteral<float> p(s);
    return Quantity<float, p.kg, p.m, p.s>(p.factor);
}

} // namespace literals

using namespace literals;

template <typename T, int8_t kg, int8_t m, int8_t s>
class Quantity {
    T value;

    template <typename, int8_t, int8_t, int8_t>
    friend class Quantity;

public:
    // TODO: private constructors
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

#define DefineQuantity(M, L, T, NAME) \
    template <typename X> \
    using NAME = Quantity<X, M, L, T>;

DefineQuantity(0,  0,  0, Dimensionless)
DefineQuantity(0,  0,  1, Time)
DefineQuantity(0,  0, -1, AngularVelocity)
DefineQuantity(0,  0, -2, AngularAcceleration)
DefineQuantity(0,  1,  0, Length)
DefineQuantity(0,  1, -1, Speed)
DefineQuantity(0,  1, -2, Acceleration)
DefineQuantity(0,  2,  0, Area)
DefineQuantity(0,  3,  0, Volume)
DefineQuantity(1,  0,  0, Mass)
DefineQuantity(1,  1, -1, Momentum)
DefineQuantity(1,  1, -2, Force)
DefineQuantity(1,  2,  0, MomentOfInertia)
DefineQuantity(1,  2, -1, AngularMomentum)
DefineQuantity(1,  2, -2, Energy)
DefineQuantity(1,  2, -3, Power)
DefineQuantity(1, -1, -1, Viscosity)
DefineQuantity(1, -1, -2, Pressure)
DefineQuantity(1, -3,  0, Density)

#undef DefineQuantity
#define AliasQuantity(NAME, QUANTITY) \
    template <typename X> \
    using NAME = QUANTITY<X>;

AliasQuantity(Distance, Length)
AliasQuantity(Frequency, AngularVelocity)
AliasQuantity(Stress, Pressure)
AliasQuantity(Torque, Energy)
AliasQuantity(Velocity, Speed)
AliasQuantity(Weight, Force)
AliasQuantity(Work, Energy)

#undef AliasQuantity

} // namespace units
