Description
-----------

Flexible Nelder-Mead-Optimizer with a simple interface.

- Comes with a small test suite that ensures correct behaviour
- Templatable to vectors of arbitrary dimension

Usage
-----

    float precision = 0.001;
    int dimension = 2;
    NelderMeadOptimizer o(dimension, precision);

    // request a simplex to start with
    Vector v(0.5, 0.5);
    o.insert(v);
    o.insert(Vector(0.1, 0.1));
    o.insert(Vector(0.2, 0.7));
    
    while (!o.done()) {
        v = o.step(v, f(v));
    }

See `t/test.cpp` for more examples. 

The optimizer works with float values. It is `done()` when each pair of candidates is at most `precision` away from each other.
