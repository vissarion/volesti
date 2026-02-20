// VolEsti (volume computation and sampling library)

// Copyright (c) 2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Original C++ code from https://github.com/rzrsk/vaidya-walk by Raaz Dwivedi.

// Modified by Alexandros Manochis to be integrated in volesti, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef PWALK_UTIL_MATH_FUNCTIONS_HPP_
#define PWALK_UTIL_MATH_FUNCTIONS_HPP_

#include <Eigen/Dense>
#include <random>
#include <cmath>

// Define random number generator type
typedef std::mt19937 rng_t;

template <typename Dtype>
void sample_gaussian(const int n, const Dtype a,
                     const Dtype sigma, Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& r) {
    static rng_t gen(1234567);
    static std::normal_distribution<Dtype> random_distribution(a, sigma);

    for (int i = 0; i < n; ++i) {
        r[i] = random_distribution(gen);
    }
}

template <typename Dtype>
Dtype rng_uniform(const Dtype a, const Dtype b) {
    static rng_t gen(1234567);
    static std::uniform_real_distribution<Dtype> random_distribution(a, b);

    return random_distribution(gen);
}

template <typename Dtype>
Dtype gaussian_density(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& x, 
                       const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& mu,
                       const Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>& sqrt_cov) {
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> c = sqrt_cov * (x - mu);
    return std::exp(-0.5*c.dot(c)) * sqrt_cov.determinant();
}


#endif // PWALK_UTIL_MATH_FUNTIONS_HPP_
