// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NEW_BASIC_SAMPLING_FEATURES_H
#define NEW_BASIC_SAMPLING_FEATURES_H

#include <chrono>

/////////////////// Random numbers generator
///
/// \tparam RNGType
/// \tparam NT
/// \tparam Ts

template <typename RNGType, typename NT, int ... Ts>
struct BoostRandomNumberGenerator;

template <typename RNGType, typename NT>
struct BoostRandomNumberGenerator<RNGType, NT>
{
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

    void set_seed(unsigned rng_seed){
        _rng.seed(rng_seed);
    }

private :
    RNGType _rng;
    boost::random::uniform_real_distribution<NT> _urdist;
    boost::random::uniform_int_distribution<> _uidist;
    boost::random::normal_distribution<NT> _ndist;
};


template <typename RNGType, typename NT, int Seed>
struct BoostRandomNumberGenerator<RNGType, NT, Seed>
{
    BoostRandomNumberGenerator(int d)
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

    void set_seed(unsigned rng_seed){
        _rng.seed(rng_seed);
    }

private :
    RNGType _rng;
    boost::random::uniform_real_distribution<NT> _urdist;
    boost::random::uniform_int_distribution<> _uidist;
    boost::random::normal_distribution<NT> _ndist;
};



///
///////////////////// Random generators' policies

struct PushBackWalkPolicy
{
    template <typename PointList, typename Point>
    void apply(PointList &randPoints,
               Point &p) const
    {
        randPoints.push_back(p);
    }
};

template <typename BallPoly>
struct CountingWalkPolicy
{
    CountingWalkPolicy(unsigned int const& nump_PBSmall, BallPoly const& PBSmall)
            :   _nump_PBSmall(nump_PBSmall)
            ,   _PBSmall(PBSmall)
    {}

    template <typename PointList, typename Point>
    void apply(PointList &randPoints,
               Point &p)
    {
        if (_PBSmall.second().is_in(p) == -1)//is in
        {
            randPoints.push_back(p);
            ++_nump_PBSmall;
        }
    }

    unsigned int get_nump_PBSmall() const
    {
        return _nump_PBSmall;
    }

private :
    unsigned int _nump_PBSmall;
    BallPoly _PBSmall;
};


/////////////////// Random walk helpers
///


template <typename Point>
struct GetDirection
{
    typedef typename Point::FT NT;

    template <typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              RandomNumberGenerator &rng)
    {
        NT normal = NT(0);
        Point p(dim);
        NT* data = p.pointerToData();

        for (unsigned int i=0; i<dim; ++i)
        {
            *data = rng.sample_ndist();
            normal += *data * *data;
            data++;
        }

        normal = NT(1)/std::sqrt(normal);
        p *= normal;

        return p;
    }
};

template <typename Point>
struct GetPointInDsphere
{
    template <typename NT, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              NT const& radius,
                              RandomNumberGenerator &rng)
    {
        Point p = GetDirection<Point>::apply(dim, rng);
        NT U = rng.sample_urdist();
        U = std::pow(U, NT(1)/(NT(dim)));
        p *= radius * U;
        return p;
    }
};

template <typename Point>
struct GetPointOnDsphere
{
    template <typename NT, typename RandomNumberGenerator>
    inline static Point apply(unsigned int const& dim,
                              NT const& radius,
                              RandomNumberGenerator &rng)
    {
        Point p = GetDirection<Point>::apply(dim, rng);
        if (radius != 0) p *= radius;
        return p;
    }
};





#endif
