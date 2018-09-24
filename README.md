Description
-----------

Flexible Nelder-Mead-Optimizer with a simple interface, that tries to maximize the target function.

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

License
-------

This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

See LICENSE file for a copy of the GNU General Public License.

Copyright (C) 2013  Sebastian Morr <sebastian@morr.cc>
