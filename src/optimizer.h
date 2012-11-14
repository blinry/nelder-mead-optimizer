#include <ctime>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
using namespace std;

// Float vector with standard operations
class Vector {
public:
    Vector() {
    }
    Vector(float c0, float c1) {
        coords.push_back(c0);
        coords.push_back(c1);
    }
    Vector(float c0, float c1, float c2) {
        coords.push_back(c0);
        coords.push_back(c1);
        coords.push_back(c2);
    }

    // add more constructors when N gets > 3

    float& operator[](int i) {
        return coords[i];
    }
    float at(int i) const {
        return coords[i];
    }
    int dimension() const {
        return coords.size();
    }
    void prepare(int size) {
        for (int i=0; i<size; i++) {
            coords.push_back(0);
        }
    }
    Vector operator+(Vector other) {
        Vector result;
        result.prepare(dimension());
        for (int i=0; i<dimension(); i++) {
            result[i] = coords[i] + other[i];
        }
        return result;
    }
    void operator+=(Vector other) {
        for (int i=0; i<dimension(); i++) {
            coords[i] += other[i];
        }
    }
    Vector operator-(Vector other) {
        Vector result;
        result.prepare(dimension());
        for (int i=0; i<dimension(); i++) {
            result[i] = coords[i] - other[i];
        }
        return result;
    }
    bool operator==(Vector other) {
        if (dimension() != other.dimension()) {
            return false;
        }
        for (int i=0; i<dimension(); i++) {
            if (other[i] != coords[i]) {
                return false;
            }
        }
        return true;
    }
    Vector operator*(float factor) {
        Vector result;
        result.prepare(dimension());
        for (int i=0; i<dimension(); i++) {
            result[i] = coords[i]*factor;
        }
        return result;
    }
    Vector operator/(float factor) {
        Vector result;
        result.prepare(dimension());
        for (int i=0; i<dimension(); i++) {
            result[i] = coords[i]/factor;
        }
        return result;
    }
    void operator/=(float factor) {
        for (int i=0; i<dimension(); i++) {
            coords[i] /= factor;
        }
    }
    bool operator<(const Vector other) const {
        for (int i=0; i<dimension(); i++) {
            if (at(i) < other.at(i))
                return false;
            else if (at(i) > other.at(i))
                return true;
        }
        return false;
    }
    float length() {
        float sum = 0;
        for (int i=0; i<dimension(); i++) {
            sum += coords[i]*coords[i];
        }
        return pow(sum, 0.5f);
    }
private:
    vector<float> coords;
};

// This class stores known values for vectors. It throws unknown vectors.
class ValueDB {
    public:
        ValueDB() {
        }
        float lookup(Vector vec) {
            if (!contains(vec)) {
                throw vec;
            } else {
                return values[vec];
            }
        }
        void insert(Vector vec, float value) {
            values[vec] = value;
        }
    private:
        bool contains(Vector vec) {
            map<Vector, float>::iterator it = values.find(vec); // TODO add tolerance
            return it != values.end();
        }
        map<Vector, float> values;
};

class NelderMeadOptimizer {
    public:
        NelderMeadOptimizer(int dimension, float termination_distance=0.001) {
            this->dimension = dimension;
            srand(time(NULL));
            alpha = 1;
            gamma = 2;
            rho = -0.5;
            sigma = 0.5;
            this->termination_distance = termination_distance;
        }
        // used in `step` to sort the vectors
        bool operator()(const Vector& a, const Vector& b) {
            return db.lookup(a) < db.lookup(b);
        }
        // termination criteria: each pair of vectors in the simplex has to
        // have a distance of at most `termination_distance`
        bool done() {
            if (vectors.size() < dimension) {
                return false;
            }
            for (int i=0; i<dimension+1; i++) {
                for (int j=0; j<dimension+1; j++) {
                    if (i==j) continue;
                    if ((vectors[i]-vectors[j]).length() > termination_distance) {
                        return false;
                    }
                }
            }
            return true;
        }
        void insert(Vector vec) {
            if (vectors.size() < dimension+1) {
                vectors.push_back(vec);
            }
        }
        Vector step(Vector vec, float score) {
            db.insert(vec, score);
            try {
                if (vectors.size() < dimension+1) {
                    vectors.push_back(vec);
                }

                // otherwise: optimize!
                if (vectors.size() == dimension+1) {
                    while(!done()) {
                        sort(vectors.begin(), vectors.end(), *this);
                        Vector cog; // center of gravity
                        cog.prepare(dimension);
                        for (int i = 1; i<=dimension; i++) {
                            cog += vectors[i];
                        }
                        cog /= dimension;
                        Vector best = vectors[dimension];
                        Vector worst = vectors[0];
                        Vector second_worst = vectors[1];
                        // reflect
                        Vector reflected = cog + (cog - worst)*alpha;
                        if (f(reflected) > f(second_worst) && f(reflected) < f(best)) {
                            vectors[0] = reflected;
                        } else if (f(reflected) > f(best)) {
                            // expand
                            Vector expanded = cog + (cog - worst)*gamma;
                            if (f(expanded) > f(reflected)) {
                                vectors[0] = expanded;
                            } else {
                                vectors[0] = reflected;
                            }
                        } else {
                            // contract
                            Vector contracted = cog + (cog - worst)*rho;
                            if (f(contracted) > f(worst)) {
                                vectors[0] = contracted;
                            } else {
                                for (int i=0; i<dimension; i++) {
                                    vectors[i] = best + (vectors[i] - best)*sigma;
                                }
                            }
                        }
                    }

                    // algorithm is terminating, output: simplex' center of gravity
                    Vector cog;
                    for (int i = 0; i<=dimension; i++) {
                        cog += vectors[i];
                    }
                    return cog/(dimension+1);
                } else {
                    // as long as we don't have enough vectors, request random ones,
                    // with coordinates between 0 and 1. If you want other start vectors,
                    // simply ignore these and use `step` on the vectors you want.
                    Vector result;
                    result.prepare(dimension);
                    for (int i = 0; i<dimension; ++i) {
                        result[i] = 0.001*(rand()%1000);
                    }
                    return result;
                }
            } catch (Vector v) {
                return v;
            }
        }
    private:
        float f(Vector vec) {
            return db.lookup(vec);
        }
        int dimension;
        float alpha, gamma, rho, sigma;
        float termination_distance;
        vector<Vector> vectors;
        ValueDB db;
};
