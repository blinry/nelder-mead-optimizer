#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../src/optimizer.h"

TEST_CASE("vector/operations", "") {
    Vector a(1,2,3);
    Vector b(1,2,3);
    Vector d;

    d = a+b;
    Vector c(2,4,6);
    REQUIRE(d == c);

    d = a*3;
    Vector c2(3,6,9);
    REQUIRE(d == c2);

    d = c2/3;
    REQUIRE(d == a);

    d = c-b;
    REQUIRE(d == a);
}

TEST_CASE("nmo/create", "") {
    NelderMeadOptimizer o(2);
}

TEST_CASE("nmo/syntax", "") {
    NelderMeadOptimizer o(2);

    Vector a(0.5, 0.5);
    Vector new_values = o.step(a, 1);
}

// test all operations in a complete run
TEST_CASE("nmo/operations", "") {
    NelderMeadOptimizer o(2);
    Vector a(0,0);
    Vector b(1,0);
    Vector c(0,1);
    Vector d(1,1);
    Vector e(1.5,1.5);
    Vector f(0,2);
    Vector g(0.75,0.5);
    Vector h(1,0.5);
    Vector i(0.5,1);

    o.step(b,0.1);
    o.step(a,0);
    Vector result;
    result.prepare(2);
    o.step(c,0.5);
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
float f(Vector v) {
    float x = v[0];
    float y = v[1];
    return ((-x*x*x*x+4.5*x*x+2)/pow(2.71828,2*y*y));
}

TEST_CASE("nmo/functionality", "") {
    NelderMeadOptimizer o(2, 0.001);

    // horrible start values
    o.insert(Vector(2, 1));
    o.insert(Vector(2.001, 0));
    o.insert(Vector(1000000, -200));

    Vector v(2, 1);

    while (!o.done()) {
        float score = f(v);
        v = o.step(v, score);
    }

    float tolerance = 0.001;
    REQUIRE(abs(v[0]) > 1.5-tolerance);
    REQUIRE(abs(v[0]) < 1.5+tolerance);
    REQUIRE(abs(v[1]) > -tolerance);
    REQUIRE(abs(v[1]) < +tolerance);
}
