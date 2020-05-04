// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "vars.h"
#include "hpolytope.h"
#include "vpolytope.h"
#include "zpolytope.h"
#include "samplers.h"
#include "rounding.h"
#include "vpolyintersectvpoly.h"
#include "extractMatPoly.h"

//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param P A convex polytope (H- or V-representation or zonotope).
//'
//' @section warning:
//' Do not use this function.
//'
//' @return A numerical matrix that describes the rounded polytope and contains the round value.
// [[Rcpp::export]]
Rcpp::List rounding (Rcpp::Reference P){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef VPolytope<Point, RNGType > Vpolytope;
    typedef Zonotope<Point> zonotope;
    typedef IntersectionOfVpoly<Vpolytope> InterVP;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    Hpolytope HP;
    Vpolytope VP;
    zonotope ZP;
    InterVP VPcVP;

    bool cdhr=false, rdhr = false, ball_walk = false, billiard = false;
    unsigned int n = P.field("dimension"), walkL, type = P.field("type");;

    std::pair <Point, NT> InnerBall;
    Rcpp::NumericMatrix Mat;

    if (type == 1) {
        walkL = 10 + 10*n;
        cdhr = true;
    } else {
        walkL = 5;
        billiard = true;
    }

    switch (type) {
        case 1: {
            // Hpolytope
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            InnerBall = HP.ComputeInnerBall();
            //if (billiard && diam < 0.0) HP.comp_diam(diam, InnerBall.second);
            break;
        }
        case 2: {
            VP.init(n, Rcpp::as<MT>(P.field("V")), VT::Ones(Rcpp::as<MT>(P.field("V")).rows()));
            InnerBall = VP.ComputeInnerBall();
            break;
        }
        case 3: {
            // Zonotope
            ZP.init(n, Rcpp::as<MT>(P.field("G")), VT::Ones(Rcpp::as<MT>(P.field("G")).rows()));
            InnerBall = ZP.ComputeInnerBall();
            break;
        }
        case 4: {
            throw Rcpp::exception("volesti does not support rounding for this representation currently.");
        }
    }

    std::pair< std::pair<MT, VT>, NT > round_res;

    switch (type) {
        case 1: {
            round_polytope(Polytope &P, std::pair<Point,NT> &InnerBall,
            const unsigned int &walk_length, RNG &rng)
            round_res = rounding_min_ellipsoid<MT, VT>(HP, InnerBall, var);
            Mat = extractMatPoly(HP);
            break;
        }
        case 2: {
            round_res = rounding_min_ellipsoid<MT, VT>(VP, InnerBall, var);
            Mat = extractMatPoly(VP);
            break;
        }
        case 3: {
            round_res = rounding_min_ellipsoid<MT, VT>(ZP, InnerBall, var);
            Mat = extractMatPoly(ZP);
            break;
        }
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Mat, Rcpp::Named("T") = Rcpp::wrap(round_res.first.first),
                              Rcpp::Named("shift") = Rcpp::wrap(round_res.first.second),
                              Rcpp::Named("round_value") = round_res.second);

}
