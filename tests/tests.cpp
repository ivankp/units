#include "test.hpp"

#include "Units.hpp"

#include <type_traits>

using namespace units;

TEST(Size) {
    static_assert(sizeof(Quantity<double,0,1,0>) == sizeof(double));
    static_assert(sizeof(Quantity<float ,0,1,0>) == sizeof(float ));
}

TEST(ImplicitConversion) {
    Quantity quantity = 5.;
    double value = quantity;
    TEST_EQ(value, 5.);
}

TEST(Addition) {
    Quantity<double,1,0,0> a = 5.;
    Quantity<float,1,0,0> b = 2.;
    auto c = a + b;
    static_assert(std::is_same_v<decltype(c), Quantity<double,1,0,0>>);
}

TEST(Addition2) {
    Quantity<double,1,0,0> a = 5.;
    Quantity<float,1,0,0> b = 2.;
    a += b;
}

TEST(Subtraction) {
    Quantity<double,1,0,0> a = 5.;
    Quantity<float,1,0,0> b = 2.;
    auto c = a - b;
    static_assert(std::is_same_v<decltype(c), Quantity<double,1,0,0>>);
}

TEST(Subtraction2) {
    Quantity<double,1,0,0> a = 5.;
    Quantity<float,1,0,0> b = 2.;
    a -= b;
}

TEST(Multiplication) {
    Quantity<double,0,1,0> length = 5.;
    Quantity<float,0,0,1> time = 2.;
    auto product = length * time;
    static_assert(std::is_same_v<decltype(product), Quantity<double,0,1,1>>);
}

TEST(Division) {
    Quantity<double,0,1,0> length = 5.;
    Quantity<float,0,0,1> time = 2.f;
    auto speed = length / time;
    static_assert(std::is_same_v<decltype(speed), Quantity<double,0,1,-1>>);
    double value = speed / Quantity<int,0,1,-1>(1);
    TEST_EQ(value, 2.5);
}

TEST(MultiplicationByFactor) {
    Length<double> l1 = 5;
    auto l2 = 2. * l1;
    static_assert(std::is_same_v<decltype(l2), Length<double>>);
    auto l3 = l1 * 3.;
    static_assert(std::is_same_v<decltype(l3), Length<double>>);
    TEST_EQ(double(l2 / l1), 2);
    TEST_EQ(double(l3 / l1), 3);
}

TEST(Exponentiation) {
    Quantity<float,0,1,0> length = 5.f;
    auto area = pow<2>(length);
    static_assert(std::is_same_v<decltype(area), Quantity<float,0,2,0>>);
    double value = area / Quantity<int,0,2,0>(1);
    TEST_EQ(value, 25);
}

TEST(CompoundQuantity) {
    Mass<double> m = 7;
    Length<double> l = 3;
    Time<double> t = 5;
    auto v = l/t;
    static_assert(std::is_same_v<decltype(v), Speed<double>>);
    [](Energy<double> E){
        TEST_EQ(double(E / Energy<double>(1)), 2.52);
    }(m * pow<2>(v));
}
