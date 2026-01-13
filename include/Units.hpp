#pragma once

#include <cmath>
#include <type_traits>

namespace units {

struct Dimensions {
    using type = int;
    type mass = 0, length = 0, time = 0;

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

    constexpr Dimensions operator*(type n) const noexcept {
        // TODO: check for overflow
        return { mass * n, length * n, time * n };
    }

    constexpr Dimensions& operator*=(type n) noexcept {
        // TODO: check for overflow
        mass *= n;
        length *= n;
        time *= n;
        return *this;
    }

    constexpr Dimensions operator/(type n) const noexcept {
        // TODO: check for remainder == 0
        return { mass / n, length / n, time / n };
    }

    constexpr Dimensions& operator/=(type n) noexcept {
        // TODO: check for remainder == 0
        mass /= n;
        length /= n;
        time /= n;
        return *this;
    }
};

namespace detail {

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

template <typename T>
struct ParsedLiteral {
    T factor = 1;
    Dimensions d { };

    enum Category { alph, num, op, none };

    constexpr ParsedLiteral(const char* a, const char* const e) {
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
struct IsQuantity : std::false_type { };

template <typename T, Dimensions d>
struct IsQuantity<Quantity<T, d>> : std::true_type { };

template <typename T>
concept AQuantity = IsQuantity<T>::value;

template <typename T>
concept NotAQuantity = !IsQuantity<T>::value;

namespace literals {

template <detail::StringLiteral s>
constexpr auto operator ""_u() noexcept {
    constexpr detail::ParsedLiteral<double> p(s.begin(), s.end());
    return Quantity<double, p.d>(p.factor);
}

template <detail::StringLiteral s>
constexpr auto operator ""_uf() noexcept {
    constexpr detail::ParsedLiteral<float> p(s.begin(), s.end());
    return Quantity<float, p.d>(p.factor);
}

} // namespace literals

template <typename T, Dimensions d>
class Quantity {
    T value;

    template <typename, Dimensions>
    friend class Quantity;

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
    friend constexpr auto operator+(Quantity l, Quantity<R, d> r) noexcept
    -> Quantity<decltype(l.value + r.value), d>
    {
        return l.value + r.value;
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
    friend constexpr auto operator-(Quantity l, Quantity<R, d> r) noexcept
    -> Quantity<decltype(l.value - r.value), d>
    {
        return l.value - r.value;
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
    friend constexpr auto operator*(Quantity l, Quantity<R, Rd> r) noexcept
    -> Quantity<decltype(l.value * r.value), d + Rd>
    {
        return l.value * r.value;
    }

    friend constexpr auto operator*(Quantity l, NotAQuantity auto r) noexcept
    -> Quantity<decltype(l.value * r), d>
    {
        return l.value * r;
    }

    friend constexpr auto operator*(NotAQuantity auto l, Quantity r) noexcept
    -> Quantity<decltype(l * r.value), d>
    {
        return l * r.value;
    }

    constexpr Quantity& operator*=(auto r) noexcept {
        value *= r.value;
        return *this;
    }

    // Division ----------------------------------------------------------------

    template <typename R, Dimensions Rd>
    friend constexpr auto operator/(Quantity l, Quantity<R, Rd> r) noexcept
    -> Quantity<decltype(l.value / r.value), d - Rd>
    {
        return l.value / r.value;
    }

    friend constexpr auto operator/(Quantity l, NotAQuantity auto r) noexcept
    -> Quantity<decltype(l.value / r), d>
    {
        return l.value / r;
    }

    friend constexpr auto operator/(NotAQuantity auto l, Quantity r) noexcept
    -> Quantity<decltype(l / r.value), -d>
    {
        return l / r.value;
    }

    constexpr Quantity& operator/=(auto r) noexcept {
        value /= r.value;
        return *this;
    }

    // Exponentiation ----------------------------------------------------------

    template <typename Dimensions::type n>
    friend constexpr Quantity<T, d * n> pow(Quantity q) noexcept {
        return std::pow(q.value, n);
    }
};

template <typename T>
Quantity(T x) -> Quantity<T, {}>;

#define DefineQuantity(M, L, T, NAME) \
    template <typename X> \
    using NAME = Quantity<X, {M, L, T}>;

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
