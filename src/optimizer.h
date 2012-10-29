#include <ctime>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
using namespace std;

// Float vector with standard operations
template <unsigned int N>
class Vector {
public:
    Vector() {
        for (int i=0; i<N; i++) {
            coords[i] = 0;
        }
    }
    Vector(float c0, float c1) {
        coords[0] = c0;
        coords[1] = c1;
    }
    Vector(float c0, float c1, float c2) {
        coords[0] = c0;
        coords[1] = c1;
        coords[2] = c2;
    }

    // add more constructors when N gets > 3

    float& operator[](int i) {
        return coords[i];
    }
    float at(int i) const {
        return coords[i];
    }
    Vector<N> operator+(Vector<N> other) {
        Vector<N> result;
        for (int i=0; i<N; i++) {
            result[i] = coords[i] + other[i];
        }
        return result;
    }
    void operator+=(Vector<N> other) {
        for (int i=0; i<N; i++) {
            coords[i] += other[i];
        }
    }
    Vector<N> operator-(Vector<N> other) {
        Vector<N> result;
        for (int i=0; i<N; i++) {
            result[i] = coords[i] - other[i];
        }
        return result;
    }
    bool operator==(Vector<N> other) {
        for (int i=0; i<N; i++) {
            if (other[i] != coords[i]) {
                return false;
            }
        }
        return true;
    }
    Vector<N> operator*(float factor) {
        Vector<N> result;
        for (int i=0; i<N; i++) {
            result[i] = coords[i]*factor;
        }
        return result;
    }
    Vector<N> operator/(float factor) {
        Vector<N> result;
        for (int i=0; i<N; i++) {
            result[i] = coords[i]/factor;
        }
        return result;
    }
    void operator/=(float factor) {
        for (int i=0; i<N; i++) {
            coords[i] /= factor;
        }
    }
    bool operator<(const Vector<N> other) const {
        for (int i=0; i<N; i++) {
            if (at(i) < other.at(i))
                return false;
            else if (at(i) > other.at(i))
                return true;
        }
        return false;
    }
    float length() {
        float sum = 0;
        for (int i=0; i<N; i++) {
            sum += coords[i]*coords[i];
        }
        return pow(sum, 0.5);
    }
private:
    float coords[N];
};

typedef Vector<2> Vec2f;
typedef Vector<3> Vec3f;

// This class stores known values for vectors. It throws unknown vectors.
template <unsigned int N>
class ValueDB {
    public:
        ValueDB() {
        }
        float lookup(Vector<N> vec) {
            if (!contains(vec)) {
                throw vec;
            } else {
                return values[vec];
            }
        }
        void insert(Vector<N> vec, float value) {
            values[vec] = value;
        }
    private:
        bool contains(Vector<N> vec) {
            typename map<Vector<N>, float>::iterator it = values.find(vec); // TODO add tolerance
            return it != values.end();
        }
        map<Vector<N>, float> values;
};

template <unsigned int N>
class NelderMeadOptimizer {
    public:
        NelderMeadOptimizer(float termination_distance=0.001) {
            srand(time(NULL));
            alpha = 1;
            gamma = 2;
            rho = -0.5;
            sigma = 0.5;
            this->termination_distance = termination_distance;
        }
        // used in `step` to sort the vectors
        bool operator()(const Vector<N>& a, const Vector<N>& b) {
            return db.lookup(a) < db.lookup(b);
        }
        // termination criteria: each pair of vectors in the simplex has to
        // have a distance of at most `termination_distance`
        bool done() {
            if (vectors.size() < N) {
                return false;
            }
            for (int i=0; i<N+1; i++) {
                for (int j=0; j<N+1; j++) {
                    if (i==j) continue;
                    if ((vectors[i]-vectors[j]).length() > termination_distance) {
                        return false;
                    }
                }
            }
            return true;
        }
        Vector<N> step(Vector<N> vec, float score) {
            db.insert(vec, score);
            try {
                if (vectors.size() < N+1) {
                    vectors.push_back(vec);
                }

                // otherwise: optimize!
                if (vectors.size() == N+1) {
                    while(!done()) {
                        sort(vectors.begin(), vectors.end(), *this);
                        Vector<N> cog; // center of gravity
                        for (int i = 1; i<=N; i++) {
                            cog += vectors[i];
                        }
                        cog /= N;
                        Vector<N> best = vectors[N];
                        Vector<N> worst = vectors[0];
                        Vector<N> second_worst = vectors[1];
                        // reflect
                        Vector<N> reflected = cog + (cog - worst)*alpha;
                        if (f(reflected) > f(second_worst) && f(reflected) < f(best)) {
                            vectors[0] = reflected;
                        } else if (f(reflected) > f(best)) {
                            // expand
                            Vector<N> expanded = cog + (cog - worst)*gamma;
                            if (f(expanded) > f(reflected)) {
                                vectors[0] = expanded;
                            } else {
                                vectors[0] = reflected;
                            }
                        } else {
                            // contract
                            Vector<N> contracted = cog + (cog - worst)*rho;
                            if (f(contracted) > f(worst)) {
                                vectors[0] = contracted;
                            } else {
                                for (int i=0; i<N; i++) {
                                    vectors[i] = best + (vectors[i] - best)*sigma;
                                }
                            }
                        }
                    }

                    // algorithm is terminating, output: simplex' center of gravity
                    Vector<N> cog;
                    for (int i = 0; i<=N; i++) {
                        cog += vectors[i];
                    }
                    return cog/(N+1);
                } else {
                    // as long as we don't have enough vectors, request random ones,
                    // with coordinates between 0 and 1. If you want other start vectors,
                    // simply ignore these and use `step` on the vectors you want.
                    Vector<N> result;
                    for (int i = 0; i<N; ++i) {
                        result[i] = 0.001*(rand()%1000);
                    }
                    return result;
                }
            } catch (Vector<N> v) {
                return v;
            }
        }
    private:
        float f(Vector<N> vec) {
            return db.lookup(vec);
        }
        float alpha, gamma, rho, sigma;
        float termination_distance;
        vector<Vector<N> > vectors;
        ValueDB<N> db;
};
