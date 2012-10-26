#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <cmath>
#include "../src/optimizer.cpp"

TEST_CASE("vector/operations", "") {
    Vector<3> a(1,2,3);
    Vector<3> b(1,2,3);
    Vec3f d;

    d = a+b;
    Vector<3> c(2,4,6);
    REQUIRE(d == c);

    d = a*3;
    Vector<3> c2(3,6,9);
    REQUIRE(d == c2);

    d = c2/3;
    REQUIRE(d == a);

    d = c-b;
    REQUIRE(d == a);
}

TEST_CASE("nmo/create", "") {
    NelderMeadOptimizer<2> o;
}

TEST_CASE("nmo/syntax", "") {
    NelderMeadOptimizer<2> o;

    Vec2f a(0.5, 0.5);
    Vec2f new_values = o.step(a, 1);
}

// test all operations in a complete run
TEST_CASE("nmo/operations", "") {
    NelderMeadOptimizer<2> o;
    Vec2f a(0,0);
    Vec2f b(1,0);
    Vec2f c(0,1);
    Vec2f d(1,1);
    Vec2f e(1.5,1.5);
    Vec2f f(0,2);
    Vec2f g(0.75,0.5);
    Vec2f h(1,0.5);
    Vec2f i(0.5,1);

    o.step(b,0.1);
    o.step(a,0);
    Vec2f result;
    result = o.step(c,0.5);
    REQUIRE(result == d);
    result = o.step(result,1);
    REQUIRE(result == e);
    result = o.step(result,0);
    REQUIRE(result == f);
    result = o.step(result,0);
    REQUIRE(result == g);
    result = o.step(result,0);
    REQUIRE(result == i);
}

// a function with two maxima at (1.5|0) and (-1.5|0),
// looks like camel humps ;-)
float f(Vec2f v) {
    float x = v[0];
    float y = v[1];
    return ((-x*x*x*x+4.5*x*x+2)/pow(2.71828,2*y*y));
}

TEST_CASE("nmo/functionality", "") {
    NelderMeadOptimizer<2> o(0.001);

    // horrible start values
    Vec2f a(2, 1);
    Vec2f b(2.001, 0);
    Vec2f c(1000000, -200);

    o.step(a, f(a));
    o.step(b, f(b));
    Vec2f v = o.step(c, f(c));

    while (!o.done()) {
        v = o.step(v, f(v));
    }

    float tolerance = 0.001;
    REQUIRE(abs(v[0]) > 1.5-tolerance);
    REQUIRE(abs(v[0]) < 1.5+tolerance);
    REQUIRE(abs(v[1]) > -tolerance);
    REQUIRE(abs(v[1]) < +tolerance);
}
