#include "test.hpp"

#include "Units.hpp"

#include <array>
#include <type_traits>

struct vec2 {
    double x, y;
};

constexpr double operator*(const vec2& a, const vec2& b) noexcept {
    return a.x * b.x + a.y * b.y;
}

using units::Quantity, units::MakeQuantity;
using namespace units::literals;
using namespace units::quantities;

TEST(Pow) {
    using units::detail::Pow;
    {
        const auto x = "x";
        auto r = Pow<0>(x); // anything to the power 0 is 1
        static_assert(std::is_same_v<decltype(r), int>);
        TEST_EQ(r, 1);
    } {
        const auto x = 3;
        auto r = Pow<4>(x);
        static_assert(std::is_same_v<decltype(r), int>);
        TEST_EQ(r, 81);
    } {
        const auto x = 3;
        auto r = Pow<-3>(x);
        static_assert(std::is_same_v<decltype(r), double>);
        TEST_EQ(r, 1./27);
    } {
        const vec2 x { 3., 5. };
        auto r = Pow<2>(x);
        static_assert(std::is_same_v<decltype(r), double>);
        TEST_EQ(r, 34.);
    } {
        const vec2 x { 3., 5. };
        auto r = Pow<-2>(x);
        static_assert(std::is_same_v<decltype(r), double>);
        TEST_EQ(r, 1./34.);
    } {
        auto r = Pow<15>(2.f);
        static_assert(std::is_same_v<decltype(r), float>);
        TEST_EQ(r, 32'768);
    } {
        auto r = Pow<31>(2ul);
        static_assert(std::is_same_v<decltype(r), unsigned long>);
        TEST_EQ(r, 2'147'483'648);
    }
}

TEST(Size) {
    static_assert(sizeof(Quantity<{0,1,0},double>) == sizeof(double));
    static_assert(sizeof(Quantity<{0,1,0},float >) == sizeof(float ));
}

TEST(ImplicitConversion) {
    Quantity nounits = 5;
    static_assert(std::is_same_v<decltype(nounits), Quantity<{0,0,0},int>>);
    double value = nounits;
    TEST_EQ(value, 5);
}

TEST(Addition) {
    auto a = MakeQuantity<{1,0,0},double>(5.);
    auto b = MakeQuantity<{1,0,0},float>(2.);
    auto c = a + b;
    static_assert(std::is_same_v<decltype(c), Quantity<{1,0,0},double>>);
}

TEST(Addition2) {
    auto a = MakeQuantity<{1,0,0},double>(5.);
    auto b = MakeQuantity<{1,0,0},float>(2.);
    auto c = a + b;
    a += b;
    TEST_EQ(a, c);
}

TEST(Subtraction) {
    auto a = MakeQuantity<{1,0,0},double>(5.);
    auto b = MakeQuantity<{1,0,0},float>(2.);
    auto c = a - b;
    static_assert(std::is_same_v<decltype(c), Quantity<{1,0,0},double>>);
}

TEST(Subtraction2) {
    auto a = MakeQuantity<{1,0,0},double>(5.);
    auto b = MakeQuantity<{1,0,0},float>(2.);
    auto c = a - b;
    a -= b;
    TEST_EQ(a, c);
}

TEST(Comparison) {
    auto a = MakeQuantity<{1,0,0},double>(5.);
    auto b = MakeQuantity<{1,0,0},float>(2.);
    TEST_TRUE(a != b);
    TEST_TRUE(b < a);
    TEST_TRUE(a >= b);
    TEST_TRUE(b < a);
    TEST_TRUE(a >= b);
}

TEST(Multiplication) {
    auto length = MakeQuantity<{0,1,0},double>(5.);
    auto time = MakeQuantity<{0,0,1},float>(2.);
    auto product = length * time;
    static_assert(std::is_same_v<decltype(product), Quantity<{0,1,1},double>>);
}

TEST(Division) {
    auto length = MakeQuantity<{0,1,0},double>(5.);
    auto time = MakeQuantity<{0,0,1},float>(2.f);
    auto speed = length / time;
    static_assert(std::is_same_v<decltype(speed), Quantity<{0,1,-1},double>>);
    double value = speed / MakeQuantity<{0,1,-1}>(1);
    TEST_EQ(value, 2.5);
}

TEST(MultiplicationWithFactors) {
    Length<double> l1 = MakeQuantity<{0,1,0}>(5);
    auto l2 = 2. * l1;
    static_assert(std::is_same_v<decltype(l2), Length<double>>);
    auto l3 = l1 * 3.;
    static_assert(std::is_same_v<decltype(l3), Length<double>>);
    TEST_EQ(double(l2 / l1), 2);
    TEST_EQ(double(l3 / l1), 3);
}

TEST(DivisionWithFactors) {
    Length<double> l1 = MakeQuantity<{0,1,0}>(5);
    auto l2 = 2 / l1;
    static_assert(std::is_same_v<decltype(l2), Quantity<{0,-1,0},double>>);
    auto l3 = l1 / 4;
    static_assert(std::is_same_v<decltype(l3), Length<double>>);
    TEST_EQ(double(l2 * l1), 2);
    TEST_EQ(double(l3 / l1), 0.25);
}

TEST(AdditionInPlace) {
    Length<double> l1 = MakeQuantity<{0,1,0}>(5);
    l1 += MakeQuantity<{0,1,0},float>(3);
    auto l2 = MakeQuantity<{0,1,0}>(8);
    TEST_EQ(l1, l2);
}

TEST(SubtractionInPlace) {
    Length<double> l1 = MakeQuantity<{0,1,0}>(5);
    l1 -= MakeQuantity<{0,1,0},float>(3);
    auto l2 = MakeQuantity<{0,1,0}>(2);
    TEST_EQ(l1, l2);
}

TEST(MultiplicationInPlace) {
    Length<double> l1 = MakeQuantity<{0,1,0}>(5);
    const auto l2 = l1;
    l1 *= 2;
    TEST_EQ(double(l1 / l2), 2);
}

TEST(DivisionInPlace) {
    Length<double> l1 = MakeQuantity<{0,1,0}>(5);
    const auto l2 = l1;
    l1 /= 2;
    TEST_EQ(double(l1 / l2), 0.5);
}

TEST(ExponentiationPow) {
    {
        auto length = MakeQuantity<{0,1,0},float>(5.f);
        auto area = Pow<2>(length);
        static_assert(std::is_same_v<decltype(area), Quantity<{0,2,0},float>>);
        double value = area / MakeQuantity<{0,2,0}>(1);
        TEST_EQ(value, 25);
    } {
        auto v = MakeQuantity<{0,1,-1}>(vec2{ 3, 7 });
        auto speed2 = Pow<2>(v);
        static_assert(std::is_same_v<decltype(speed2), Quantity<{0,2,-2},double>>);
        double value = speed2 / MakeQuantity<{0,2,-2}>(1);
        TEST_EQ(value, 58.);
    }
}

TEST(ExponentiationRoot) {
    {
        auto volume = MakeQuantity<{0,3,0},float>(125.f);
        auto length = Root<3>(volume);
        static_assert(std::is_same_v<decltype(length), Quantity<{0,1,0},float>>);
        double value = length / MakeQuantity<{0,1,0}>(1);
        TEST_EQ(value, 5);
    } {
        auto v = MakeQuantity<{0,1,-1}>(vec2{ 3, 4 });
        auto speedInv = Root<-2>(Pow<2>(v));
        static_assert(std::is_same_v<decltype(speedInv), Quantity<{0,-1,1},double>>);
        double value = speedInv / MakeQuantity<{0,-1,1}>(1);
        TEST_EQ(value, 0.2);
    }
}

TEST(CompoundQuantity) {
    Mass<double> m = MakeQuantity<{1,0,0}>(7);
    Length<double> l = MakeQuantity<{0,1,0}>(3);
    Time<double> t = MakeQuantity<{0,0,1}>(5);
    auto v = l/t;
    static_assert(std::is_same_v<decltype(v), Speed<double>>);
    [](Energy<double> E){
        TEST_EQ(double(E / MakeQuantity<{1,2,-2}>(1)), 2.52);
    }(m * Pow<2>(v));
}

TEST(Inverse) {
    Time<double> time = MakeQuantity<{0,0,1}>(5);
    auto freq = 1 / time;
    static_assert(std::is_same_v<decltype(freq), Frequency<double>>);
}

TEST(LiteralDouble) {
    static_assert(std::is_same_v<decltype(1 * "kg"_u), Mass<double>>);
    static_assert(std::is_same_v<decltype(1 * "m"_u), Length<double>>);
    static_assert(std::is_same_v<decltype(1 * "s"_u), Time<double>>);
    static_assert(std::is_same_v<decltype("m"_u / "s"_u), Speed<double>>);
}

TEST(LiteralFloat) {
    static_assert(std::is_same_v<decltype(1 * "kg"_uf), Mass<float>>);
    static_assert(std::is_same_v<decltype(1 * "m"_uf), Length<float>>);
    static_assert(std::is_same_v<decltype(1 * "s"_uf), Time<float>>);
}

void type(auto);

TEST(Literals2) {
    type("m/s"_u);
    // static_assert(std::is_same_v<decltype("m/s"_u), Speed<double>>);
}
