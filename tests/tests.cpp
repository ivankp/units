#include "Units.hpp"

#include <type_traits>

#include "test.hpp"

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

TEST(Literals2) {
    static_assert(std::is_same_v<decltype("m"_u), Length<double>>);
    static_assert(std::is_same_v<decltype("m m"_u), Area<double>>);
    static_assert(std::is_same_v<decltype("m m m"_u), Volume<double>>);
    static_assert(std::is_same_v<decltype("m s"_u), Quantity<{0,1,1},double>>);
    static_assert(std::is_same_v<decltype("m1"_u), Length<double>>);
    static_assert(std::is_same_v<decltype("m2"_u), Area<double>>);
    static_assert(std::is_same_v<decltype("m3"_u), Volume<double>>);
    static_assert(std::is_same_v<decltype("m0"_u), Quantity<{0,0,0},double>>);

    static_assert("m0"_u == 1.);
    static_assert("kg0"_u == 1.);
    static_assert(""_u == 1.);

    static_assert(std::is_same_v<decltype("s-1"_u), Frequency<double>>);
    static_assert(std::is_same_v<decltype("s -1"_u), Frequency<double>>);
    static_assert(std::is_same_v<decltype("s- 1"_u), Frequency<double>>);
    static_assert(std::is_same_v<decltype("s - 1"_u), Frequency<double>>);
    static_assert(std::is_same_v<decltype("m s-1"_u), Speed<double>>);
    static_assert(std::is_same_v<decltype("s-1 m"_u), Speed<double>>);
    static_assert(std::is_same_v<decltype("s-1m"_u), Speed<double>>);
    static_assert(std::is_same_v<decltype("s- 1m"_u), Speed<double>>);

    static_assert(std::is_same_v<decltype("m/s"_u), Speed<double>>);
    static_assert("m/m"_u == 1);
    static_assert(std::is_same_v<decltype("m/m3"_u), Quantity<{0,-2,0},double>>);
    static_assert(std::is_same_v<decltype("/s"_u), Frequency<double>>);
    static_assert(std::is_same_v<decltype("/s-1"_u), Time<double>>);
    static_assert(std::is_same_v<decltype(Root<2>("m2"_u)), Length<double>>);

    static_assert(int("psi / Pa"_u + 0.5) == 6895);

    static_assert("lbf/in2"_u == "psi"_u);
}

TEST(LiteralsPrefixes) {
    static_assert("Qs"_u  == 1e30  * "s"_u);
    static_assert("Rs"_u  == 1e27  * "s"_u);
    static_assert("Ys"_u  == 1e24  * "s"_u);
    static_assert("Zs"_u  == 1e21  * "s"_u);
    static_assert("Es"_u  == 1e18  * "s"_u);
    static_assert("Ps"_u  == 1e15  * "s"_u);
    static_assert("Ts"_u  == 1e12  * "s"_u);
    static_assert("Gs"_u  == 1e9   * "s"_u);
    static_assert("Ms"_u  == 1e6   * "s"_u);
    static_assert("ks"_u  == 1e3   * "s"_u);
    static_assert("hs"_u  == 1e2   * "s"_u);
    static_assert("das"_u == 1e1   * "s"_u);
    static_assert("ds"_u  == 1e-1  * "s"_u);
    static_assert("cs"_u  == 1e-2  * "s"_u);
    static_assert("ms"_u  == 1e-3  * "s"_u);
    static_assert("us"_u  == 1e-6  * "s"_u);
    static_assert("ns"_u  == 1e-9  * "s"_u);
    static_assert("ps"_u  == 1e-12 * "s"_u);
    static_assert("fs"_u  == 1e-15 * "s"_u);
    static_assert("as"_u  == 1e-18 * "s"_u);
    static_assert("zs"_u  == 1e-21 * "s"_u);
    static_assert("ys"_u  == 1e-24 * "s"_u);
    static_assert("rs"_u  == 1e-27 * "s"_u);
    static_assert("qs"_u  == 1e-30 * "s"_u);

    static_assert("Qm / m"_u  == 1e30 );
    static_assert("Rm / m"_u  == 1e27 );
    static_assert("Ym / m"_u  == 1e24 );
    static_assert("Zm / m"_u  == 1e21 );
    static_assert("Em / m"_u  == 1e18 );
    static_assert("Pm / m"_u  == 1e15 );
    static_assert("Tm / m"_u  == 1e12 );
    static_assert("Gm / m"_u  == 1e9  );
    static_assert("Mm / m"_u  == 1e6  );
    static_assert("km / m"_u  == 1e3  );
    static_assert("hm / m"_u  == 1e2  );
    static_assert("dam / m"_u == 1e1  );
    static_assert("dm / m"_u  == 1e-1 );
    static_assert("cm / m"_u  == 1e-2 );
    static_assert("mm / m"_u  == 1e-3 );
    static_assert("um / m"_u  == 1e-6 );
    static_assert("nm / m"_u  == 1e-9 );
    static_assert("pm / m"_u  == 1e-12);
    static_assert("fm / m"_u  == 1e-15);
    static_assert("am / m"_u  == 1e-18);
    static_assert("zm / m"_u  == 1e-21);
    static_assert("ym / m"_u  == 1e-24);
    static_assert("rm / m"_u  == 1e-27);
    static_assert("qm / m"_u  == 1e-30);

    static_assert("m / Qm"_u  == double(1 / 1e30L ));
    static_assert("m / Rm"_u  == double(1 / 1e27L ));
    static_assert("m / Ym"_u  == double(1 / 1e24L ));
    static_assert("m / Zm"_u  == double(1 / 1e21L ));
    static_assert("m / Em"_u  == double(1 / 1e18L ));
    static_assert("m / Pm"_u  == double(1 / 1e15L ));
    static_assert("m / Tm"_u  == double(1 / 1e12L ));
    static_assert("m / Gm"_u  == double(1 / 1e9L  ));
    static_assert("m / Mm"_u  == double(1 / 1e6L  ));
    static_assert("m / km"_u  == double(1 / 1e3L  ));
    static_assert("m / hm"_u  == double(1 / 1e2L  ));
    static_assert("m / dam"_u == double(1 / 1e1L  ));
    static_assert("m / dm"_u  == double(1 / 1e-1L ));
    static_assert("m / cm"_u  == double(1 / 1e-2L ));
    static_assert("m / mm"_u  == double(1 / 1e-3L ));
    static_assert("m / um"_u  == double(1 / 1e-6L ));
    static_assert("m / nm"_u  == double(1 / 1e-9L ));
    static_assert("m / pm"_u  == double(1 / 1e-12L));
    static_assert("m / fm"_u  == double(1 / 1e-15L));
    static_assert("m / am"_u  == double(1 / 1e-18L));
    static_assert("m / zm"_u  == double(1 / 1e-21L));
    static_assert("m / ym"_u  == double(1 / 1e-24L));
    static_assert("m / rm"_u  == double(1 / 1e-27L));
    static_assert("m / qm"_u  == double(1 / 1e-30L));
}

TEST(LiteralsTime) {
    static_assert("min"_u == 60 * "s"_u);
    static_assert("h"_u == 60 * "min"_u);
    static_assert("d"_u == 24 * "h"_u);
}

TEST(Literals3) {
    static_assert(1e3 * "g"_u == "kg"_u);
    static_assert(1e6 * "mm"_u == "km"_u);
    static_assert(1e6 * "Pa"_u == "MPa"_u);

    static_assert("mL"_u / "cm3"_u == (double(1e-3L * 1e-3L) / double(1e-2L * 1e-2L * 1e-2L)));
    static_assert("L"_u / "dm3"_u == (double(1e-3L) / double(0.1L * 0.1L * 0.1L)));

    static_assert("Pa"_u == "N/m2"_u);

    static_assert("J"_u == "kg m2/s2"_u);
    static_assert("J"_u == "N m"_u);

    static_assert("W"_u == "kg m2 s-3"_u);
    static_assert("W"_u == "J/s"_u);

    static_assert("Hz"_u == "s-1"_u);

    static_assert("ha"_u == 1e4 * "m2"_u);
    static_assert("ha"_u == "hm2"_u);
    static_assert("ha"_u / "yd2"_u == (1e4 / double(0.9144L * 0.9144L)));

    static_assert("t"_u == 1e3 * "kg"_u);

    static_assert("EeV"_u / "J"_u / 0.1602 > 1);
    static_assert("EeV"_u / "J"_u / 0.1602 < 1 + 1e-3);

    static_assert(1e7 * "erg"_u == "J"_u);
    static_assert(1e-2 * "erg"_u == "nJ"_u);
    static_assert(10 * "erg"_u == "uJ"_u);

    static_assert("gTNT"_u == 1e3 * "cal"_u);
    static_assert("kgTNT"_u == "Mcal"_u);

    static_assert("BTU"_u == 1.0551 * "kJ"_u);

    static_assert("kW h"_u == 3.6 * "MJ"_u);
    static_assert("Wh"_u == "W h"_u);
}

TEST(LiteralsSpaces) {
    static_assert("s"_u == "s "_u);
    static_assert("s"_u == " s"_u);
    static_assert("s"_u == " s "_u);

    static_assert("s-1"_u == "s -1"_u);
    static_assert("s-1"_u == "s- 1"_u);
    static_assert("s-1"_u == "s - 1"_u);
    static_assert("s-1"_u == "/s"_u);
    static_assert("s-1"_u == "/ s"_u);

    static_assert("s"_u == "/s-1"_u);
    static_assert("s"_u == "/ s-1"_u);
    static_assert("s"_u == "/ s -1"_u);
    static_assert("s"_u == "/ s- 1"_u);
    static_assert("s"_u == "/ s - 1"_u);
}

TEST(SimpleLiterals) {
    using namespace units::simple_literals;

    static_assert(1._cm  == "cm"_u );
    static_assert(1._d   == "d"_u  );
    static_assert(1._ft  == "ft"_u );
    static_assert(1._g   == "g"_u  );
    static_assert(1._h   == "h"_u  );
    static_assert(1._ha  == "ha"_u );
    static_assert(1._Hz  == "Hz"_u );
    static_assert(1._in  == "in"_u );
    static_assert(1._J   == "J"_u  );
    static_assert(1._kg  == "kg"_u );
    static_assert(1._kHz == "kHz"_u);
    static_assert(1._kJ  == "kJ"_u );
    static_assert(1._km  == "km"_u );
    static_assert(1._kPa == "kPa"_u);
    static_assert(1._kt  == "kt"_u );
    static_assert(1._kW  == "kW"_u );
    static_assert(1._kWh == "kWh"_u);
    static_assert(1._L   == "L"_u  );
    static_assert(1._lb  == "lb"_u );
    static_assert(1._lbf == "lbf"_u);
    static_assert(1._lbm == "lbm"_u);
    static_assert(1._m   == "m"_u  );
    static_assert(1._MHz == "MHz"_u);
    static_assert(1._mi  == "mi"_u );
    static_assert(1._min == "min"_u);
    static_assert(1._MJ  == "MJ"_u );
    static_assert(1._mL  == "mL"_u );
    static_assert(1._mm  == "mm"_u );
    static_assert(1._MPa == "MPa"_u);
    static_assert(1._mph == "mph"_u);
    static_assert(1._N   == "N"_u  );
    static_assert(1._Pa  == "Pa"_u );
    static_assert(1._psi == "psi"_u);
    static_assert(1._s   == "s"_u  );
    static_assert(1._t   == "t"_u  );
    static_assert(1._W   == "W"_u  );
    static_assert(1._Wh  == "Wh"_u );
    static_assert(1._yd  == "yd"_u );

    static_assert(122.3_cm  == 1.223e2 * "cm"_u );
    static_assert(122.3_d   == 1.223e2 * "d"_u  );
    static_assert(122.3_ft  == 1.223e2 * "ft"_u );
    static_assert(122.3_g   == 1.223e2 * "g"_u  );
    static_assert(122.3_h   == 1.223e2 * "h"_u  );
    static_assert(122.3_ha  == 1.223e2 * "ha"_u );
    static_assert(122.3_Hz  == 1.223e2 * "Hz"_u );
    static_assert(122.3_in  == 1.223e2 * "in"_u );
    static_assert(122.3_J   == 1.223e2 * "J"_u  );
    static_assert(122.3_kg  == 1.223e2 * "kg"_u );
    static_assert(122.3_kHz == 1.223e2 * "kHz"_u);
    static_assert(122.3_kJ  == 1.223e2 * "kJ"_u );
    static_assert(122.3_km  == 1.223e2 * "km"_u );
    static_assert(122.3_kPa == 1.223e2 * "kPa"_u);
    static_assert(122.3_kt  == 1.223e2 * "kt"_u );
    static_assert(122.3_kW  == 1.223e2 * "kW"_u );
    static_assert(122.3_kWh == 1.223e2 * "kWh"_u);
    static_assert(122.3_L   == 1.223e2 * "L"_u  );
    static_assert(122.3_lb  == 1.223e2 * "lb"_u );
    static_assert(122.3_lbf == 1.223e2 * "lbf"_u);
    static_assert(122.3_lbm == 1.223e2 * "lbm"_u);
    static_assert(122.3_m   == 1.223e2 * "m"_u  );
    static_assert(122.3_MHz == 1.223e2 * "MHz"_u);
    static_assert(122.3_mi  == 1.223e2 * "mi"_u );
    static_assert(122.3_min == 1.223e2 * "min"_u);
    static_assert(122.3_MJ  == 1.223e2 * "MJ"_u );
    static_assert(122.3_mL  == 1.223e2 * "mL"_u );
    static_assert(122.3_mm  == 1.223e2 * "mm"_u );
    static_assert(122.3_MPa == 1.223e2 * "MPa"_u);
    static_assert(122.3_mph == 1.223e2 * "mph"_u);
    static_assert(122.3_N   == 1.223e2 * "N"_u  );
    static_assert(122.3_Pa  == 1.223e2 * "Pa"_u );
    static_assert(122.3_psi == 1.223e2 * "psi"_u);
    static_assert(122.3_s   == 1.223e2 * "s"_u  );
    static_assert(122.3_t   == 1.223e2 * "t"_u  );
    static_assert(122.3_W   == 1.223e2 * "W"_u  );
    static_assert(122.3_Wh  == 1.223e2 * "Wh"_u );
    static_assert(122.3_yd  == 1.223e2 * "yd"_u );
}
