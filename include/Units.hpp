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
constexpr auto Pow(const T& base, const R& result = Identity{})
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
    double factor = 1;
    Dimensions d { };

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

template <typename T, Dimensions d>
class Quantity;

template <typename T>
constexpr bool isQuantity = false;

template <typename T, Dimensions d>
constexpr bool isQuantity<Quantity<T, d>> = true;

template <typename T>
concept AQuantity = isQuantity<T>;

template <typename T>
concept NotAQuantity = !isQuantity<T>;

namespace literals {

template <detail::StringLiteral s>
constexpr auto operator ""_u() noexcept {
    constexpr detail::LiteralParser p(s.begin(), s.end());
    return Quantity<double, p.d>(p.factor);
}

template <detail::StringLiteral s>
constexpr auto operator ""_uf() noexcept {
    constexpr detail::LiteralParser p(s.begin(), s.end());
    return Quantity<float, p.d>(float(p.factor));
}

} // namespace literals

namespace detail {

template <Dimensions d, typename T>
constexpr Quantity<T, d> Make(T v) noexcept {
    return v;
}

} // namespace detail

template <NotAQuantity T, Dimensions d>
class Quantity<T, d> {
    T value;

    template <typename, Dimensions>
    friend class Quantity;

    template <Dimensions dd, typename V>
    friend constexpr Quantity<V, dd> detail::Make(V v) noexcept;

public:
    // TODO: private constructors
    constexpr Quantity() noexcept: value{} { }
    constexpr Quantity(T value) noexcept: value(value) { }

    constexpr operator T() const noexcept requires (!d) {
        return value;
    }

    // Comparison --------------------------------------------------------------

    template <typename R>
    constexpr auto operator<=>(Quantity<R, d> r) const noexcept {
        return value <=> r.value;
    }

    template <typename R>
    constexpr bool operator==(Quantity<R, d> r) const noexcept {
        return value == r.value;
    }

    // Addition ----------------------------------------------------------------

    template <typename R>
    constexpr auto operator+(Quantity<R, d> r) const noexcept
    -> Quantity<decltype(value + r.value), d>
    {
        return value + r.value;
    }

    template <typename R>
    constexpr Quantity& operator+=(Quantity<R, d> r) noexcept {
        value += r.value;
        return *this;
    }

    constexpr Quantity operator+() const noexcept {
        return +value;
    }

    // Subtraction -------------------------------------------------------------

    template <typename R>
    constexpr auto operator-(Quantity<R, d> r) const noexcept
    -> Quantity<decltype(value - r.value), d>
    {
        return value - r.value;
    }

    template <typename R>
    constexpr Quantity& operator-=(Quantity<R, d> r) noexcept {
        value -= r.value;
        return *this;
    }

    constexpr Quantity operator-() const noexcept {
        return -value;
    }

    // Multiplication ----------------------------------------------------------

    template <typename R, Dimensions Rd>
    constexpr auto operator*(Quantity<R, Rd> r) const noexcept
    -> Quantity<decltype(value * r.value), d + Rd>
    {
        return value * r.value;
    }

    template <NotAQuantity R>
    constexpr auto operator*(R r) const noexcept
    -> Quantity<decltype(value * r), d>
    {
        return value * r;
    }

    template <NotAQuantity L, typename R, Dimensions Rd>
    friend constexpr auto operator*(L l, Quantity<R, Rd> r) noexcept
    -> Quantity<decltype(l * r.value), Rd>;

    constexpr Quantity& operator*=(auto r) noexcept {
        value *= r;
        return *this;
    }

    // Division ----------------------------------------------------------------

    template <typename R, Dimensions Rd>
    constexpr auto operator/(Quantity<R, Rd> r) const noexcept
    -> Quantity<decltype(value / r.value), d - Rd>
    {
        return value / r.value;
    }

    template <NotAQuantity R>
    constexpr auto operator/(R r) const noexcept
    -> Quantity<decltype(value / r), d>
    {
        return value / r;
    }

    template <NotAQuantity L, typename R, Dimensions Rd>
    friend constexpr auto operator/(L l, Quantity<R, Rd> r) noexcept
    -> Quantity<decltype(l / r.value), -Rd>;

    constexpr Quantity& operator/=(auto r) noexcept {
        value /= r;
        return *this;
    }

    // Exponentiation ----------------------------------------------------------

    template <dim_t n>
    friend constexpr auto Pow(Quantity q) noexcept(std::is_arithmetic_v<T>) {
        return detail::Make<(d * n)>(detail::Pow<n>(q.value));
    }

    template <dim_t n>
    friend constexpr auto Root(Quantity q) noexcept {
        return detail::Make<(d / n)>(detail::Root<n>(q.value));
    }
};

template <NotAQuantity L, typename R, Dimensions Rd>
constexpr auto operator*(L l, Quantity<R, Rd> r) noexcept
-> Quantity<decltype(l * r.value), Rd>
{
    return l * r.value;
}

template <NotAQuantity L, typename R, Dimensions Rd>
constexpr auto operator/(L l, Quantity<R, Rd> r) noexcept
-> Quantity<decltype(l / r.value), -Rd>
{
    return l / r.value;
}

template <typename T>
Quantity(T x) -> Quantity<T, {}>;

#define DefineQuantity(M, L, T, NAME) \
    template <typename X> \
    using NAME = Quantity<X, {dim_t(M), dim_t(L), dim_t(T)}>;

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
