#pragma once

#include <cmath>
#include <string_view>
#include <type_traits>

namespace units {

using dim_t = int;

struct Dimensions {
    dim_t mass = 0, length = 0, time = 0;

    constexpr operator bool() const noexcept {
        return mass || length || time;
    }

    constexpr Dimensions operator+(const Dimensions& r) const noexcept {
        // TODO: check for overflow
        return { mass + r.mass, length + r.length, time + r.time };
    }

    constexpr Dimensions& operator+=(const Dimensions& r) noexcept {
        // TODO: check for overflow
        mass += r.mass;
        length += r.length;
        time += r.time;
        return *this;
    }

    constexpr Dimensions operator-(const Dimensions& r) const noexcept {
        // TODO: check for underflow
        return { mass - r.mass, length - r.length, time - r.time };
    }

    constexpr Dimensions operator-() const noexcept {
        return { -mass, -length, -time };
    }

    constexpr Dimensions& operator-=(const Dimensions& r) noexcept {
        // TODO: check for underflow
        mass -= r.mass;
        length -= r.length;
        time -= r.time;
        return *this;
    }

    constexpr Dimensions operator*(dim_t n) const noexcept {
        // TODO: check for overflow
        return { mass * n, length * n, time * n };
    }

    constexpr Dimensions& operator*=(dim_t n) noexcept {
        // TODO: check for overflow
        mass *= n;
        length *= n;
        time *= n;
        return *this;
    }

    constexpr Dimensions operator/(dim_t n) const noexcept {
        // TODO: check for remainder == 0
        return { mass / n, length / n, time / n };
    }

    constexpr Dimensions& operator/=(dim_t n) noexcept {
        // TODO: check for remainder == 0
        mass /= n;
        length /= n;
        time /= n;
        return *this;
    }
};

namespace detail {

struct Identity { };

template <typename T>
constexpr bool isIdentity = false;

template <>
constexpr bool isIdentity<Identity> = true;

template <dim_t n, typename T, typename R = Identity>
constexpr auto Pow(const T& base, const R& result = {})
noexcept(std::is_arithmetic_v<T>)
{
    if constexpr (n == 0) {
        return 1;
    } else if constexpr (n < 0) {
        return 1. / Pow<-n>(base, result);
    } else {
        constexpr bool odd = n & 1;
        constexpr dim_t half = n >> 1;

        if constexpr (!half) {
            if constexpr (!odd) {
                return result;
            } else if constexpr (isIdentity<R>) {
                return base;
            } else {
                return result * base;
            }
        } else {
            if constexpr (!odd) {
                return Pow<half>(base * base, result);
            } else if constexpr (isIdentity<R>) {
                return Pow<half>(base * base, base);
            } else {
                return Pow<half>(base * base, result * base);
            }
        }
    }
}

template <dim_t n>
constexpr auto Root(const auto& x) noexcept {
    if constexpr (n == 0) {
        return 1;
    } else if constexpr (n < 0) {
        return 1. / Root<-n>(x);
    } else if constexpr (n == 1) {
        return x;
    } if constexpr (n == 2) {
        return std::sqrt(x);
    } if constexpr (n == 3) {
        return std::cbrt(x);
    } else {
        return std::pow(x, 1. / n);
    }
}

template <unsigned N>
struct StringLiteral {
    char s[N]{};

    constexpr StringLiteral(const char(&arr)[N]) noexcept {
        for (unsigned i = 0; i < N; ++i)
            s[i] = arr[i];
    }

    constexpr const char* begin() const noexcept {
        return s;
    }

    constexpr const char* end() const noexcept {
        return s + N;
    }
};

struct LiteralParser {
    Dimensions d { };
    double factor = 1;

    enum Category { alph, num, op, none };

    constexpr LiteralParser(const char* a, const char* const e) {
        const char* b = a;
        // Category cat = none, cat2 = none;

        for (; b < e; ++b) {
            const char c = *b;

            /*
            if (c == ' ' || c == '\t') {
                cat2 = none;
            } else if (('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z')) {
                cat2 = alph;
            } else {
                // TODO: error
            }

            if (cat != cat2) {
                switch (cat) {
                    case alph: {
                    } break;
                    case num: {
                    } break;
                }
                a = b;
                cat = cat2;
            }
            */

            if ((c < 'A' || 'Z' < c) && (c < 'a' || 'z' < c)) {
                // TODO: write a custom comparison function (maybe lambda)
                const std::string_view v(a, b);
                if (v == "kg") {
                    d.mass += 1;
                } else if (v == "m") {
                    d.length += 1;
                } else if (v == "s") {
                    d.time += 1;
                }
                break; // TODO: complete parsing algorithm
            }
        }
    }
};

} // namespace detail

template <Dimensions d, typename T>
class Quantity;

template <typename T>
constexpr bool isQuantity = false;

template <Dimensions d, typename T>
constexpr bool isQuantity<Quantity<d, T>> = true;

template <typename T>
concept AQuantity = isQuantity<T>;

template <typename T>
concept NotAQuantity = !isQuantity<T>;

template <Dimensions d, typename T>
constexpr Quantity<d, T> MakeQuantity(T value) noexcept {
    return value;
}

namespace literals {

template <detail::StringLiteral s>
constexpr auto operator ""_u() noexcept {
    constexpr detail::LiteralParser p(s.begin(), s.end());
    return MakeQuantity<p.d>(p.factor);
}

template <detail::StringLiteral s>
constexpr auto operator ""_uf() noexcept {
    constexpr detail::LiteralParser p(s.begin(), s.end());
    return MakeQuantity<p.d>(float(p.factor));
}

} // namespace literals

template <Dimensions d, NotAQuantity T>
class Quantity<d, T> {
    T value;

    template <Dimensions, typename>
    friend class Quantity;

    template <Dimensions d_, typename T_>
    friend constexpr Quantity<d_, T_> MakeQuantity(T_) noexcept;

    // Prevent default initialization
    Quantity() = delete;

    // Only friends can use implicit conversion for dimensionful quantities
    constexpr Quantity(T value) noexcept requires ((bool)d) : value(value) { }

public:
    // Implicit converting constructor for dimensionless quantities
    constexpr Quantity(T value) noexcept requires (!d) : value(value) { }

    template <typename R>
    constexpr Quantity(Quantity<d, R> r) noexcept : value(r.value) { }

    template <typename R>
    constexpr Quantity& operator=(Quantity<d, R> r) noexcept {
        value = r.value;
    }

    // Conversion operator for dimensionless quantities
    constexpr operator T() const noexcept requires (!d) {
        return value;
    }

    // Comparison --------------------------------------------------------------

    template <typename R>
    constexpr auto operator<=>(Quantity<d, R> r) const noexcept {
        return value <=> r.value;
    }

    template <typename R>
    constexpr bool operator==(Quantity<d, R> r) const noexcept {
        return value == r.value;
    }

    // Addition ----------------------------------------------------------------

    template <typename R>
    constexpr auto operator+(Quantity<d, R> r) const noexcept
    -> Quantity<d, decltype(value + r.value)>
    {
        return value + r.value;
    }

    template <typename R>
    constexpr Quantity& operator+=(Quantity<d, R> r) noexcept {
        value += r.value;
        return *this;
    }

    constexpr Quantity operator+() const noexcept {
        return +value;
    }

    // Subtraction -------------------------------------------------------------

    template <typename R>
    constexpr auto operator-(Quantity<d, R> r) const noexcept
    -> Quantity<d, decltype(value - r.value)>
    {
        return value - r.value;
    }

    template <typename R>
    constexpr Quantity& operator-=(Quantity<d, R> r) noexcept {
        value -= r.value;
        return *this;
    }

    constexpr Quantity operator-() const noexcept {
        return -value;
    }

    // Multiplication ----------------------------------------------------------

    template <Dimensions Rd, typename R>
    constexpr auto operator*(Quantity<Rd, R> r) const noexcept
    -> Quantity<d + Rd, decltype(value * r.value)>
    {
        return value * r.value;
    }

    template <NotAQuantity R>
    constexpr auto operator*(R r) const noexcept
    -> Quantity<d, decltype(value * r)>
    {
        return value * r;
    }

    template <NotAQuantity L, Dimensions Rd, typename R>
    friend constexpr auto operator*(L l, Quantity<Rd, R> r) noexcept
    -> Quantity<Rd, decltype(l * r.value)>;

    constexpr Quantity& operator*=(auto r) noexcept {
        value *= r;
        return *this;
    }

    // Division ----------------------------------------------------------------

    template <Dimensions Rd, typename R>
    constexpr auto operator/(Quantity<Rd, R> r) const noexcept
    -> Quantity<d - Rd, decltype(value / r.value)>
    {
        return value / r.value;
    }

    template <NotAQuantity R>
    constexpr auto operator/(R r) const noexcept
    -> Quantity<d, decltype(value / r)>
    {
        return value / r;
    }

    template <NotAQuantity L, Dimensions Rd, typename R>
    friend constexpr auto operator/(L l, Quantity<Rd, R> r) noexcept
    -> Quantity<-Rd, decltype(l / r.value)>;

    constexpr Quantity& operator/=(auto r) noexcept {
        value /= r;
        return *this;
    }

    // Exponentiation ----------------------------------------------------------

    template <dim_t n>
    friend constexpr auto Pow(Quantity q) noexcept(std::is_arithmetic_v<T>) {
        return MakeQuantity<(d * n)>(detail::Pow<n>(q.value));
    }

    template <dim_t n>
    friend constexpr auto Root(Quantity q) noexcept {
        return MakeQuantity<(d / n)>(detail::Root<n>(q.value));
    }
};

template <NotAQuantity L, Dimensions Rd, typename R>
constexpr auto operator*(L l, Quantity<Rd, R> r) noexcept
-> Quantity<Rd, decltype(l * r.value)>
{
    return l * r.value;
}

template <NotAQuantity L, Dimensions Rd, typename R>
constexpr auto operator/(L l, Quantity<Rd, R> r) noexcept
-> Quantity<-Rd, decltype(l / r.value)>
{
    return l / r.value;
}

template <typename T>
Quantity(T x) -> Quantity<{}, T>;

namespace quantities {

#define DefineQuantity(M, L, T, NAME) \
    template <typename X> \
    using NAME = Quantity<{dim_t(M), dim_t(L), dim_t(T)}, X>;

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

} // namespace quantities

} // namespace units
