// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NEW_VOLUME_H
#define NEW_VOLUME_H


#include <iterator>
#include <vector>
#include <list>
#include <math.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include "vars.h"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "ball.h"
#include "ballintersectconvex.h"
#include "vpolyintersectvpoly.h"
#include "samplers.h"
#include "rounding.h"
#include "gaussian_samplers.h"
#include "gaussian_annealing.h"

#include "khach.h"

/////////////////// Random Walks

// ball walk with uniform target distribution
template
<
    typename Polytope,
    typename RNGType
>
struct BallWalk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    BallWalk(Polytope P)
    {
        P.ComputeInnerBall();
    }

    BallWalk(BallPolytope)
    {}

    template
    <
        typename BallPolytope,
        typename PointList,
        typename WalkPolicy
    >
    void apply(BallPolytope &P,
               Point &p,   // a point to start
               const unsigned int rnum,
               const unsigned int walk_length,
               PointList &randPoints,
               WalkPolicy &policy)
    {
        const NT delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());

        for (auto i=0; i<rnum; ++i)
        {
            for (auto j=0; j<walk_length; ++j)
            {
                Point y = get_point_in_Dsphere<RNGType, Point>(p.dimension(), delta);
                y = y + p;
                if (P.is_in(y)==-1) p = y;
            }
            policy.apply(randPoints, p);
        }
    }

    template
    <
        typename BallPolytope
    >
    void apply(BallPolytope &P,
               Point &p,   // a point to start
               const unsigned int walk_length)
    {
        const NT delta = ((P.InnerBall()).second * NT(4)) / NT(P.dimension());

        for (auto j=0; j<walk_length; ++j)
        {
            Point y = get_point_in_Dsphere<RNGType, Point>(p.dimension(), delta);
            y = y + p;
            if (P.is_in(y)==-1) p = y;
        }
    }
};

// random directions hit-and-run walk with uniform target distribution
template
<
    typename Polytope,
    typename RNGType
>
struct RDHRWalk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    RDHRWalk(Polytope P)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    RDHRWalk(BallPolytope P)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    template
    <
        typename BallPolytope,
        typename PointList,
        typename WalkPolicy
    >
    void apply(BallPolytope &P,
               Point &p,   // a point to start
               const unsigned int rnum,
               const unsigned int walk_length,
               PointList &randPoints,
               WalkPolicy &policy)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);

        initialize(P, p);

        for (auto i=0; i<rnum; ++i)
        {
            for (auto j=0; j<walk_length; ++j)
            {
                Point v = get_direction<RNGType, Point, NT>(p.dimension());
                std::pair<NT, NT> bpair = P.line_intersect(_p, v, _lamdas, _Av,
                                                           _lambda);
                _lambda = urdist(rng) * (bpair.first - bpair.second)
                        + bpair.second;
                _p = (_lambda * v) + _p;
            }
            policy.apply(randPoints, _p);
        }
    }

    template
    <
        typename BallPolytope
    >
    void apply(BallPolytope &P,
               Point &p,   // a point to start
               const unsigned int walk_length)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);

        for (auto j=0; j<walk_length; ++j)
        {
            Point v = get_direction<RNGType, Point, NT>(p.dimension());
            std::pair<NT, NT> bpair = P.line_intersect(_p, v, _lamdas, _Av,
                                                       _lambda);
            _lambda = urdist(rng) * (bpair.first - bpair.second)
                    + bpair.second;
            _p = (_lambda * v) + _p;
        }
        p = _p;
    }

private :

    template <typename BallPolytope>
    void initialize(BallPolytope &P, Point &p)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);

        _lamdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());

        Point v = get_direction<RNGType, Point, NT>(p.dimension());
        std::pair<NT, NT> bpair = P.line_intersect(p, v, _lamdas, _Av);
        _lambda = urdist(rng) * (bpair.first - bpair.second) + bpair.second;
        _p = (_lambda * v) + p;
    }

    Point _p;
    NT _lambda;
    typename Point::Coeff _lamdas;
    typename Point::Coeff _Av;
};


// random directions hit-and-run walk with uniform target distribution
template
<
    typename Polytope,
    typename RNGType
>
struct CDHRWalk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    CDHRWalk(Polytope P)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    CDHRWalk(BallPolytope P)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    template
    <
        typename BallPolytope,
        typename PointList,
        typename WalkPolicy
    >
    void apply(BallPolytope &P,
               Point &p,   // a point to start
               const unsigned int rnum,
               const unsigned int walk_length,
               PointList &randPoints,
               WalkPolicy &policy)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, p.dimension()-1);

        initialize(P, p);

        for (auto i=0; i<rnum; ++i)
        {
            for (auto j=0; j<walk_length; ++j)
            {
                auto rand_coord_prev = _rand_coord;
                _rand_coord = uidist(rng);
                NT kapa = urdist(rng);
                std::pair<NT, NT> bpair = P.line_intersect_coord(_p,
                                                                 _p_prev,
                                                                 _rand_coord,
                                                                 rand_coord_prev,
                                                                 _lamdas);
                _p_prev = _p;
                _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                             * (bpair.second - bpair.first));
            }
            policy.apply(randPoints, _p);
        }
    }

    template
    <
        typename BallPolytope
    >
    void apply(BallPolytope &P,
               Point &p,   // a point to start
               const unsigned int walk_length)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, p.dimension()-1);

        for (auto j=0; j<walk_length; ++j)
        {
            auto rand_coord_prev = _rand_coord;
            _rand_coord = uidist(rng);
            NT kapa = urdist(rng);
            std::pair<NT, NT> bpair = P.line_intersect_coord(_p,
                                                             _p_prev,
                                                             _rand_coord,
                                                             rand_coord_prev,
                                                             _lamdas);
            _p_prev = _p;
            _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                         * (bpair.second - bpair.first));
        }
        p = _p;
    }

private :

    template <typename BallPolytope>
    void initialize(BallPolytope &P, Point &p)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, p.dimension()-1);

        _lamdas.setZero(P.num_of_hyperplanes());
        _rand_coord = uidist(rng);
        NT kapa = urdist(rng);
        _p=p;
        std::pair<NT, NT> bpair = P.line_intersect_coord(_p, _rand_coord, _lamdas);
        _p_prev = _p;
        _p.set_coord(_rand_coord, _p[_rand_coord] + bpair.first + kapa
                    * (bpair.second - bpair.first));
    }

    unsigned int _rand_coord;
    Point _p;
    Point _p_prev;
    typename Point::Coeff _lamdas;
};


// billiard walk for uniform distribution
template
<
    typename Polytope,
    typename RNGType
>
struct BilliardWalk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    BilliardWalk(Polytope P)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    BilliardWalk(BallPolytope P)
    {
        Point center = P.InnerBall().first;
        initialize(P, center);
    }

    template
    <
        typename GenericPolytope
    >
    void apply(GenericPolytope &P,
               Point &p,   // a point to start
               const unsigned int walk_length)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);

        unsigned int n = P.dimension();
        NT T = urdist(rng) * P.ComputeDiameter();
        const NT dl = 0.995;
        NT diameter = P.ComputeDiameter();

        for (auto j=0; j<walk_length; ++j)
        {
            T = urdist(rng) * diameter;
            _v = get_direction<RNGType, Point, NT>(n);
            Point p0 = _p;

            auto it = 0;
            while (it < 10*n)
            {
                std::pair<NT, int> pbpair;
                //if (i==0 && j==0)
                if (j==0)
                {
                    pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av);
                } else {
                    pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev);
                }
                if (T <= pbpair.first) {
                    _p = (T * _v) + _p;
                    _lambda_prev = T;
                    break;
                }

                _lambda_prev = dl * pbpair.first;
                _p = (_lambda_prev * _v) + _p;
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
            if (it == 10*n) _p = p0;
        }

        p = _p;
    }

private :

    template
    <
        typename GenericPolytope
    >
    void initialize(GenericPolytope &P,
                    Point &p)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_real_distribution<> urdist(0, 1);

        unsigned int n = P.dimension();
        NT T;// = urdist(rng) * P.ComputeDiameter();
        const NT dl = 0.995;
        NT diameter = P.ComputeDiameter();

        _lambdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());

        _p = p;

        T = urdist(rng) * diameter;
        _v = get_direction<RNGType, Point, NT>(n);
        Point p0 = _p;

        auto it = 0;
        while (it < 10*n)
        {
            std::pair<NT, int> pbpair;
            pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av);
            if (T <= pbpair.first) {
                _p = (T * _v) + _p;
                _lambda_prev = T;
                break;
            }

            _lambda_prev = dl * pbpair.first;
            _p = (_lambda_prev * _v) + _p;
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);
            it++;
        }
        if (it == 10*n) _p = p0;

    }

    Point _p;
    Point _v;
    NT _lambda_prev;
    typename Point::Coeff _lambdas;
    typename Point::Coeff _Av;
};


///
/// Random generators' policies

struct PushBackWalkPolicy
{
    template <typename PointList, typename Point>
    void apply(PointList &randPoints,
               Point &p)
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

    unsigned int get_nump_PBSmall()
    {
        return _nump_PBSmall;
    }

private :
    unsigned int _nump_PBSmall;
    BallPoly _PBSmall;
};


////////////////////////////// Random Point Generators
///

template
<
    typename Walk,
    typename RNGType
>
struct RandomPointGenerator
{
    template
    <
        typename Polytope,
        typename Point,
        typename PointList,
        typename WalkPolicy
    >
    static void apply(Polytope &P,
                      Point &p,   // a point to start
                      const unsigned int rnum,
                      const unsigned int walk_length,
                      PointList &randPoints,
                      WalkPolicy &policy)
    {
        Walk walk(P);
        for (auto i=0; i<rnum; ++i)
        {
            walk.template apply(P, p, walk_length);
            policy.apply(randPoints, p);
        }
    }
};



////////////////////////////// Algorithms


// ----- ROUNDING ------ //
/*
 *     auto round = false;

    //1. Rounding of the polytope if round=true
    NT round_value=1;
    if(round){
#ifdef VOLESTI_DEBUG
        std::cout<<"\nRounding.."<<std::endl;
        double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
#endif
        std::pair<NT,NT> res_round = rounding_min_ellipsoid<NT, RNGType>(P, 1, walk);
        round_value=res_round.first;
        std::pair<Point,NT> res=P.ComputeInnerBall();
        c=res.first;
        radius=res.second;
        P.comp_diam(var.diameter, radius);
        if (var.ball_walk){
            var.delta = 4.0 * radius / NT(n);
        }

#ifdef VOLESTI_DEBUG
        double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
#endif
    }
    */
/*
        //apply rounding
        double tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
        //res_round = rounding_min_ellipsoid(P, CheBall, var);
        std::pair<NT,NT> res_round = rounding_min_ellipsoid<NT, RNGType>(P, 1, walk);
        round_value = round_value * res_round.first;
        auto ratio2 = res_round.second;
        auto ratio1 = 0.0;
        int count=1;
        //apply rounding until conditios are satisfied
        while(ratio2>ratio1 && count<=1) {
            P.ComputeInnerBall(); //compute the new chebychev center
            //res_round = rounding_min_ellipsoid(P, CheBall, var);
            res_round = rounding_min_ellipsoid<NT, RNGType>(P, 1, walk);
            round_value = round_value * res_round.first;
            ratio1=ratio2;
            ratio2 = res_round.second;
            count++;
        }
        double tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout<<"\nround value is: "<<round_value<<std::endl;
        std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
        P.ComputeInnerBall(); //compute the new chebychev center
*/

// main rounding function
template <typename NT, typename RNGType, typename Polytope, typename WalkType>
std::pair<NT, NT> rounding_min_ellipsoid(Polytope &P,
                                          unsigned int walk_length,
                                          WalkType walk) {

    typedef typename Polytope::PolytopePoint Point;
    //typedef typename Polytope::FT NT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::VT VT;

    P.ComputeInnerBall();
    const std::pair<Point,NT> InnerBall = P.InnerBall();
    auto n = P.dimension();
    auto i = 0;
    auto j = 0;
    Point c = InnerBall.first;
    NT radius = InnerBall.second;
    std::list<Point> randPoints; //ds for storing rand points

    //std::cout << radius << ",";
    //c.print();

    if (!P.get_points_for_rounding(randPoints)) {  // If P is a V-polytope then it will store its vertices in randPoints
        // If P is not a V-Polytope or number_of_vertices>20*domension
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
        Point p = get_point_in_Dsphere<RNGType, Point>(n, radius);
        p = p + c;

        //use a large walk length e.g. 1000
        //rand_point_generator(P, p, 1, 10*n, randPoints, var);
        PushBackWalkPolicy policy;
        walk.template apply<RNGType>(P, p, 1, 10*n, randPoints, policy);

        // 3. Sample points from P
        unsigned int num_of_samples = 10*n;//this is the number of sample points will used to compute min_ellipoid
        randPoints.clear();
        walk.template apply<RNGType>(P, p, num_of_samples, 10 + n / 10, randPoints, policy);

        //if (var.bill_walk) {
        //    rand_point_generator(P, p, num_of_samples, 5, randPoints, var);
        //} else {
        //    rand_point_generator(P, p, num_of_samples, 10 + n / 10, randPoints, var);
        //}
    }
    //TODO: Update diameter for Billiard and radius for Ball walk

    // Store points in a matrix to call Khachiyan algorithm for the minimum volume enclosing ellipsoid
    boost::numeric::ublas::matrix<double> Ap(n,randPoints.size());
    typename std::list<Point>::iterator rpit=randPoints.begin();

    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        for (i=0 ; i<rpit->dimension(); i++){
            Ap(i,j)=double((*rpit)[i]);
        }
    }
    boost::numeric::ublas::matrix<double> Q(n,n); //TODO: remove dependence on ublas and copy to eigen
    boost::numeric::ublas::vector<double> c2(n);
    size_t w=1000;
    KhachiyanAlgo(Ap,0.01,w,Q,c2); // call Khachiyan algorithm

    MT E(n,n);
    VT e(n);

    //Get ellipsoid matrix and center as Eigen objects
    for(unsigned int i=0; i<n; i++){
        e(i)=NT(c2(i));
        for (unsigned int j=0; j<n; j++){
            E(i,j)=NT(Q(i,j));
        }
    }


    //Find the smallest and the largest axes of the elliposoid
    Eigen::EigenSolver<MT> eigensolver(E);
    NT rel = std::real(eigensolver.eigenvalues()[0]);
    NT Rel = std::real(eigensolver.eigenvalues()[0]);
    for(unsigned int i=1; i<n; i++){
        if(std::real(eigensolver.eigenvalues()[i])<rel) rel=std::real(eigensolver.eigenvalues()[i]);
        if(std::real(eigensolver.eigenvalues()[i])>Rel) Rel=std::real(eigensolver.eigenvalues()[i]);
    }

    Eigen::LLT<MT> lltOfA(E); // compute the Cholesky decomposition of E
    MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition

    //Shift polytope in order to contain the origin (center of the ellipsoid)
    P.shift(e);

    MT L_1 = L.inverse();
    P.linear_transformIt(L_1.transpose());

    return std::pair<NT, NT> (L_1.determinant(),rel/Rel);
}


// ----- VOLUME ------ //

template
<
    typename Polytope,
    typename RNGType = boost::mt19937,
    typename WalkType = BallWalk<Polytope,RNGType>
>
double volume(Polytope &P,
              double error = 1.0,
              unsigned int walk_length = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef Ball<Point> Ball;
    typedef BallIntersectPolytope<Polytope,Ball> BallPoly;

    typedef RandomPointGenerator<WalkType, RNGType> RandomPointGenerator;

    unsigned int n = P.dimension();
    unsigned int rnum = std::pow(error, -2) * 400 * n * std::log(n);
    unsigned int n_threads = 1;

    //0. Get the Chebychev ball (largest inscribed ball) with center and radius
    P.ComputeInnerBall();
    auto InnerBall = P.InnerBall();
    Point c = InnerBall.first;
    NT radius = InnerBall.second;

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    P.shift(c.getCoefficients());
    c = Point(n);

    rnum=rnum/n_threads;
    NT vol=0;

    // Perform the procedure for a number of threads and then take the average
    for(unsigned int t=0; t<n_threads; t++){
        // 2. Generate the first random point in P
        // Perform random walk on random point in the Chebychev ball
#ifdef VOLESTI_DEBUG
        std::cout<<"\nGenerate the first random point in P"<<std::endl;
#endif
        Point p = get_point_on_Dsphere<RNGType , Point>(n, radius);
        std::list<Point> randPoints; //ds for storing rand points

        //use a large walk length e.g. 1000
        //rand_point_generator(P, p, 1, 50*n, randPoints, var);
        PushBackWalkPolicy policy;
        //walk.template apply(P, p, 1, 50*n, randPoints, policy);
        RandomPointGenerator::apply(P, p, 1, 50*n, randPoints, policy);

        // 3. Sample "rnum" points from P
#ifdef VOLESTI_DEBUG
        double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout<<"\nCompute "<<rnum<<" random points in P"<<std::endl;
#endif
        //rand_point_generator(P, p, rnum-1, walk_len, randPoints, var);
        //walk.template apply<RNGType>(P, p, rnum-1, walk_length, randPoints, policy);
        RandomPointGenerator::apply(P, p, rnum-1, walk_length, randPoints, policy);

#ifdef VOLESTI_DEBUG
        double tstop2 = (double)clock()/(double)CLOCKS_PER_SEC;
        std::cout << "First random points construction time = "
                  << tstop2 - tstart2 << std::endl;
#endif

        // 4.  Construct the sequence of balls
        // 4a. compute the radius of the largest ball
        NT current_dist, max_dist=NT(0);
        for(typename  std::list<Point>::iterator pit=randPoints.begin();
            pit!=randPoints.end(); ++pit){
            current_dist=(*pit).squared_length();
            if(current_dist>max_dist){
                max_dist=current_dist;
            }
        }
        max_dist=std::sqrt(max_dist);
#ifdef VOLESTI_DEBUG
        std::cout<<"\nFurthest distance from Chebychev point= "<<max_dist
                <<std::endl;
        std::cout<<"\nConstructing the sequence of balls"<<std::endl;
#endif

        //
        // 4b. Number of balls
        int nb1 = n * (std::log(radius)/std::log(2.0));
        int nb2 = std::ceil(n * (std::log(max_dist)/std::log(2.0)));

        std::vector<Ball> balls;

        for (int i=nb1; i<=nb2; ++i)
        {
            if (i==nb1)
            {
                balls.push_back(Ball(c,radius*radius));
                vol = (std::pow(M_PI,n/2.0)*(std::pow(balls[0].radius(), n) ) )
                       / (tgamma(n/2.0+1));
            } else {
                balls.push_back(Ball(c,std::pow(std::pow(2.0,NT(i)/NT(n)),2)));
            }
        }
        assert(!balls.empty());

#ifdef VOLESTI_DEBUG
        std::cout<<"---------"<<std::endl;
#endif

        // 5. Estimate Vol(P)

        typename std::vector<Ball>::iterator bit2=balls.end();
        bit2--;

        while(bit2!=balls.begin()){

            //each step starts with some random points in PBLarge stored
            //in list "randPoints", these points have been generated in a
            //previous step

            BallPoly PBLarge(P,*bit2);
            --bit2;
            BallPoly PBSmall(P,*bit2);

#ifdef VOLESTI_DEBUG
            std::cout<<"("<<balls.end()-bit2<<"/"<<balls.end()-balls.begin()
                     <<")Ball ratio radius="
                     <<PBLarge.second().radius()<<","
                     <<PBSmall.second().radius()<<std::endl;
            std::cout<<"Points in PBLarge="<<randPoints.size()<<std::endl;
#endif

            // choose a point in PBLarge to be used to generate more rand points
            Point p_gen = *randPoints.begin();

            // num of points in PBSmall and PBLarge
            unsigned int nump_PBSmall = 0;
            unsigned int nump_PBLarge = randPoints.size();

            //keep the points in randPoints that fall in PBSmall
            typename std::list<Point>::iterator rpit=randPoints.begin();
            while (rpit!=randPoints.end())
            {
                if (PBSmall.second().is_in(*rpit) == 0)//not in
                {
                    rpit=randPoints.erase(rpit);
                } else {
                    ++nump_PBSmall;
                    ++rpit;
                }
            }

#ifdef VOLESTI_DEBUG
            std::cout<<"Points in PBSmall="<<randPoints.size()
                     <<"\nRatio= "<<NT(nump_PBLarge)/NT(nump_PBSmall)
                     <<std::endl;
            std::cout<<"Generate "<<rnum-nump_PBLarge<<" more "<<std::endl;
#endif

            //generate more random points in PBLarge to have "rnum" in total
            //rand_point_generator(PBLarge,p_gen,rnum-nump_PBLarge,walk_len,
            //                     randPoints,PBSmall,nump_PBSmall,var,walk);

            CountingWalkPolicy<BallPoly> counting_policy(nump_PBSmall, PBSmall);

            //BilliardWalkOld<Point> walk_old;
            //walk_old.template apply<RNGType>(PBLarge, p_gen, rnum-nump_PBLarge,
            //                            walk_length, randPoints,
            //                             counting_policy);

            RandomPointGenerator::apply(PBLarge, p_gen, rnum-nump_PBLarge,
                             walk_length, randPoints,
                             counting_policy);

            nump_PBSmall = counting_policy.get_nump_PBSmall();

            vol *= NT(rnum)/NT(nump_PBSmall);

#ifdef VOLESTI_DEBUG
            std::cout<<nump_PBSmall<<"/"<<rnum<<" = "<<NT(rnum)/nump_PBSmall
                     <<"\ncurrent_vol = "<<vol
                     <<"\n--------------------------"<<std::endl;
#endif

            //don't continue in pairs of balls that are almost inside P, i.e. ratio ~= 2
        }
    }
    //vol=round_value*vol;
#ifdef VOLESTI_DEBUG
    std::cout<<"rand points = "<<rnum<<std::endl;
    std::cout<<"walk len = "<<walk_length<<std::endl;
    //std::cout<<"round_value: "<<round_value<<std::endl;
    std::cout<<"volume computed: "<<vol<<std::endl;
#endif

    P.free_them_all();
    return vol;
}

#endif
