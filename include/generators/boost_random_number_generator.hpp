// VolEsti (volume computation and sampling library)

// Copyright (c) 2020 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.


#ifndef GENERATORS_BOOST_RANDOM_NUMBER_GENERATOR_HPP
#define GENERATORS_BOOST_RANDOM_NUMBER_GENERATOR_HPP

#include <chrono>
#include <random>
#include <type_traits>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>

namespace detail {

/// Trait: is T a Boost random engine?
/// We detect this by checking for the BOOST_PREVENT_MACRO_SUBSTITUTION-based
/// min()/max() — concretely, boost engines inherit from
/// boost::random::detail::ptr_helper or similar. The simplest reliable check
/// is whether the type lives in the boost:: namespace via a specialisation.
template <typename T>
struct is_boost_engine : std::false_type {};

// Specialise for the two Boost MT engines used in volesti
template <>
struct is_boost_engine<boost::random::mt19937> : std::true_type {};
template <>
struct is_boost_engine<boost::random::mt11213b> : std::true_type {};

/// Distribution selector: pick boost::random or std:: based on engine type.
template <typename RNGType, typename NT, bool IsBoost = is_boost_engine<RNGType>::value>
struct rng_dist_traits;

template <typename RNGType, typename NT>
struct rng_dist_traits<RNGType, NT, /*IsBoost=*/true> {
    using urdist_t   = boost::random::uniform_real_distribution<NT>;
    using uidist_t   = boost::random::uniform_int_distribution<>;
    using ndist_t    = boost::random::normal_distribution<NT>;
    using expdist_t  = boost::random::exponential_distribution<NT>;
};

template <typename RNGType, typename NT>
struct rng_dist_traits<RNGType, NT, /*IsBoost=*/false> {
    using urdist_t   = std::uniform_real_distribution<NT>;
    using uidist_t   = std::uniform_int_distribution<>;
    using ndist_t    = std::normal_distribution<NT>;
    using expdist_t  = std::exponential_distribution<NT>;
};

template <typename RNG, typename NT>
inline NT sample_trunc_expdist(RNG& rng,
    typename rng_dist_traits<RNG, NT>::expdist_t& expdist)
{
    NT z;
    do { z = expdist(rng); } while (z > NT(1));
    return z;
}

} // namespace detail

/////////////////// Random numbers generator
///
/// \tparam RNGType
/// \tparam NT
/// \tparam Ts
///
/// Note: boost::mt11213b (used as default in some volume functions) has been
/// replaced by std::mt19937, which is a safe, standard drop-in with a longer
/// period. boost::mt19937 maps directly to std::mt19937.

template <typename RNGType, typename NT, int ... Ts>
struct BoostRandomNumberGenerator;

template <typename RNGType, typename NT>
struct BoostRandomNumberGenerator<RNGType, NT>
{
    using _traits    = detail::rng_dist_traits<RNGType, NT>;
    using _urdist_t  = typename _traits::urdist_t;
    using _uidist_t  = typename _traits::uidist_t;
    using _ndist_t   = typename _traits::ndist_t;
    using _expdist_t = typename _traits::expdist_t;

    BoostRandomNumberGenerator(int d)
            :   _rng(std::chrono::system_clock::now().time_since_epoch().count())
            ,   _urdist(0, 1)
            ,   _uidist(0, d-1)
            ,   _ndist(0, 1)
    {}

    NT sample_urdist()
    {
        return _urdist(_rng);
    }

    NT sample_uidist()
    {
        return _uidist(_rng);
    }

    NT sample_ndist()
    {
        return _ndist(_rng);
    }

    NT sample_trunc_expdist()
    {
        return detail::sample_trunc_expdist<RNGType, NT>(_rng, _expdist);
    }

    void set_seed(unsigned rng_seed){
        _rng.seed(rng_seed);
    }

private :
    RNGType    _rng;
    _urdist_t  _urdist;
    _uidist_t  _uidist;
    _ndist_t   _ndist;
    _expdist_t _expdist;
};


template <typename RNGType, typename NT, int Seed>
struct BoostRandomNumberGenerator<RNGType, NT, Seed>
{
    using _traits    = detail::rng_dist_traits<RNGType, NT>;
    using _urdist_t  = typename _traits::urdist_t;
    using _uidist_t  = typename _traits::uidist_t;
    using _ndist_t   = typename _traits::ndist_t;
    using _expdist_t = typename _traits::expdist_t;

    BoostRandomNumberGenerator(int d=1)
            :   _rng(Seed)
            ,   _urdist(0, 1)
            ,   _uidist(0, d-1)
            ,   _ndist(0, 1)
    {}

    NT sample_urdist()
    {
        return _urdist(_rng);
    }

    NT sample_uidist()
    {
        return _uidist(_rng);
    }

    NT sample_ndist()
    {
        return _ndist(_rng);
    }

    NT sample_trunc_expdist()
    {
        return detail::sample_trunc_expdist<RNGType, NT>(_rng, _expdist);
    }

    void set_seed(unsigned rng_seed){
        _rng.seed(rng_seed);
    }

private :
    RNGType    _rng;
    _urdist_t  _urdist;
    _uidist_t  _uidist;
    _ndist_t   _ndist;
    _expdist_t _expdist;
};

#endif // GENERATORS_BOOST_RANDOM_NUMBER_GENERATOR_HPP
