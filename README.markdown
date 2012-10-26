Description
-----------

Flexible Nelder-Mead-Optimizer with a simple interface.

- Comes with a small test suite that ensures correct behaviour
- Templatable to vectors of arbitrary dimension

Usage
-----

    float precision = 0.001;
    NelderMeadOptimizer<2> o(precision);

    Vector<2> v(0.5, 0.5);
    
    while (!o.done()) {
        v = o.step(v, f(v));
    }

See `t/test.cpp` for more examples. 

The optimizer works with float values. It is `done()` when each pair of coordinates is at most `precision` away from each other.
