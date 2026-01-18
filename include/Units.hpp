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

template <dim_t N, typename T, typename R = Identity>
constexpr auto Pow(const T& base, const R& result = {})
noexcept(std::is_arithmetic_v<T>)
{
    if constexpr (N == 0) {
        return 1;
    } else if constexpr (N < 0) {
        return 1. / Pow<-N>(base, result);
    } else {
        constexpr bool odd = N & 1;
        constexpr dim_t half = N >> 1;

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

template <dim_t N>
constexpr auto Root(const auto& x) noexcept {
    if constexpr (N == 0) {
        return 1;
    } else if constexpr (N < 0) {
        return 1. / Root<-N>(x);
    } else if constexpr (N == 1) {
        return x;
    } if constexpr (N == 2) {
        return std::sqrt(x);
    } if constexpr (N == 3) {
        return std::cbrt(x);
    } else {
        return std::pow(x, 1. / N);
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

// kg
// kg2
// kg-2
// /kg
// /kg2
// /kg-2

struct LiteralParser {
    Dimensions d { };
    double factor = 1;

    constexpr LiteralParser(const char* a, const char* const e) {
        char cat = 0, catPrev = 0;
        // bool div = false;
        bool minus = false;
        bool num = false;
        dim_t* dim = nullptr;
        int n = 1;

        auto consume = [&] {
            // if (div)
            //     n = -n;
            if (minus)
                n = -n;
            *dim += n;
            // div = false;
            minus = false;
            num = false;
            dim = nullptr;
            n = 1;
        };

        const char* b = a;
        for (; b < e; ++b) {
            const char c = *b;
            if (('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z')) {
                cat = 'a';
            } else if ('0' <= c && c <= '9') {
                cat = '0';
            } else if (c == ' ' || c == '\t' || c == '\0') {
                cat = ' ';
            } else if (c == '-') {
                cat = c;
            // } else if (c == '/') {
            //     if (div || minus)
            //         throw "Unexpected division in unit literal";
            //     cat = c;
            //     div = true;
            } else {
                throw "Unexpected character in unit literal";
            }

            if (cat == catPrev)
                continue;

            switch (catPrev) {
                case 'a': {
                    const std::string_view token(a, b);
                    if (token == "kg") {
                        dim = &d.mass;
                    } else if (token == "m") {
                        dim = &d.length;
                    } else if (token == "s") {
                        dim = &d.time;
                    } else {
                        throw "Unexpected string in unit literal";
                    }
                } break;
                case '0': {
                    n = 0;
                    for (; a < b; ++a) {
                        (n *= 10) += *a - '0';
                    }
                    num = true;
                } break;
                default: ;
            }

            switch (cat) {
                case '0': {
                    if (!dim || num)
                        throw "Unexpected number in unit literal";
                } break;
                case '-': {
                    if (!dim || minus || num)
                        throw "Unexpected minus in unit literal";
                    minus = true;
                } break;
                default: ;
            }

            if (catPrev && cat == 'a')
                consume();

            catPrev = cat;
            a = b;
        }
        if (dim)
            consume();
    }
};

} // namespace detail

template <Dimensions D, typename T>
class Quantity;

template <typename T>
constexpr bool isQuantity = false;

template <Dimensions D, typename T>
constexpr bool isQuantity<Quantity<D, T>> = true;

template <typename T>
concept AQuantity = isQuantity<T>;

template <typename T>
concept NotAQuantity = !isQuantity<T>;

template <Dimensions D, typename T>
constexpr Quantity<D, T> MakeQuantity(T value) noexcept {
    return value;
}

template <Dimensions D, NotAQuantity T>
class Quantity<D, T> {
    T value;

    template <Dimensions, typename>
    friend class Quantity;

    template <Dimensions d_, typename T_>
    friend constexpr Quantity<d_, T_> MakeQuantity(T_) noexcept;

    // Prevent default initialization
    Quantity() = delete;

    // Only friends can use implicit conversion for dimensionful quantities
    constexpr Quantity(T value) noexcept requires ((bool)D) : value(value) { }

public:
    // Implicit converting constructor for dimensionless quantities
    constexpr Quantity(T value) noexcept requires (!D) : value(value) { }

    template <typename R>
    constexpr Quantity(Quantity<D, R> r) noexcept : value(r.value) { }

    template <typename R>
    constexpr Quantity& operator=(Quantity<D, R> r) noexcept {
        value = r.value;
    }

    // Conversion operator for dimensionless quantities
    constexpr operator T() const noexcept requires (!D) {
        return value;
    }

    // Comparison --------------------------------------------------------------

    template <typename R>
    constexpr auto operator<=>(Quantity<D, R> r) const noexcept {
        return value <=> r.value;
    }

    template <typename R>
    constexpr bool operator==(Quantity<D, R> r) const noexcept {
        return value == r.value;
    }

    // Addition ----------------------------------------------------------------

    template <typename R>
    constexpr auto operator+(Quantity<D, R> r) const noexcept
    -> Quantity<D, decltype(value + r.value)>
    {
        return value + r.value;
    }

    template <typename R>
    constexpr Quantity& operator+=(Quantity<D, R> r) noexcept {
        value += r.value;
        return *this;
    }

    constexpr Quantity operator+() const noexcept {
        return +value;
    }

    // Subtraction -------------------------------------------------------------

    template <typename R>
    constexpr auto operator-(Quantity<D, R> r) const noexcept
    -> Quantity<D, decltype(value - r.value)>
    {
        return value - r.value;
    }

    template <typename R>
    constexpr Quantity& operator-=(Quantity<D, R> r) noexcept {
        value -= r.value;
        return *this;
    }

    constexpr Quantity operator-() const noexcept {
        return -value;
    }

    // Multiplication ----------------------------------------------------------

    template <Dimensions RD, typename R>
    constexpr auto operator*(Quantity<RD, R> r) const noexcept
    -> Quantity<D + RD, decltype(value * r.value)>
    {
        return value * r.value;
    }

    template <NotAQuantity R>
    constexpr auto operator*(R r) const noexcept
    -> Quantity<D, decltype(value * r)>
    {
        return value * r;
    }

    template <NotAQuantity L, Dimensions RD, typename R>
    friend constexpr auto operator*(L l, Quantity<RD, R> r) noexcept
    -> Quantity<RD, decltype(l * r.value)>;

    constexpr Quantity& operator*=(auto r) noexcept {
        value *= r;
        return *this;
    }

    // Division ----------------------------------------------------------------

    template <Dimensions RD, typename R>
    constexpr auto operator/(Quantity<RD, R> r) const noexcept
    -> Quantity<D - RD, decltype(value / r.value)>
    {
        return value / r.value;
    }

    template <NotAQuantity R>
    constexpr auto operator/(R r) const noexcept
    -> Quantity<D, decltype(value / r)>
    {
        return value / r;
    }

    template <NotAQuantity L, Dimensions RD, typename R>
    friend constexpr auto operator/(L l, Quantity<RD, R> r) noexcept
    -> Quantity<-RD, decltype(l / r.value)>;

    constexpr Quantity& operator/=(auto r) noexcept {
        value /= r;
        return *this;
    }

    // Exponentiation ----------------------------------------------------------

    template <dim_t N>
    friend constexpr auto Pow(Quantity q) noexcept(std::is_arithmetic_v<T>) {
        return MakeQuantity<(D * N)>(detail::Pow<N>(q.value));
    }

    template <dim_t N>
    friend constexpr auto Root(Quantity q) noexcept {
        return MakeQuantity<(D / N)>(detail::Root<N>(q.value));
    }
};

template <NotAQuantity L, Dimensions RD, typename R>
constexpr auto operator*(L l, Quantity<RD, R> r) noexcept
-> Quantity<RD, decltype(l * r.value)>
{
    return l * r.value;
}

template <NotAQuantity L, Dimensions RD, typename R>
constexpr auto operator/(L l, Quantity<RD, R> r) noexcept
-> Quantity<-RD, decltype(l / r.value)>
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

} // namespace units
