#include "test.hpp"

#include "Units.hpp"

#include <type_traits>

using units::Quantity;

TEST(Size) {
    static_assert(sizeof(Quantity<0,1,0,double>) == sizeof(double));
    static_assert(sizeof(Quantity<0,1,0,float >) == sizeof(float ));
}

TEST(ImplicitConversion) {
    Quantity quantity = 5.;
    double value = quantity;
    TEST_EQ(value, 5.);
}

// TEST(Meter) {
//     Quantity<0,1,0> a = 5;
//     double b = a;
//     TEST_EQ(b, 5.);
// }

TEST(Multiplication) {
    Quantity<0,1,0,double> length = 5.;
    Quantity<0,0,1,float> time = 2.;
    auto product = length * time;
    static_assert(std::is_same_v<decltype(product), Quantity<0,1,1,double>>);
}

TEST(Division) {
    Quantity<0,1,0,double> length = 5.;
    Quantity<0,0,1,float> time = 2.f;
    auto speed = length / time;
    static_assert(std::is_same_v<decltype(speed), Quantity<0,1,-1,double>>);
    double value = speed / Quantity<0,1,-1,int>(1);
    TEST_EQ(value, 2.5);
}

TEST(Power) {
    Quantity<0,1,0,float> length = 5.f;
    auto area = pow<2>(length);
    static_assert(std::is_same_v<decltype(area), Quantity<0,2,0,float>>);
    double value = area / Quantity<0,2,0,int>(1);
    TEST_EQ(value, 25);
}
