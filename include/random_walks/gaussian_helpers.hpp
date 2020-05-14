#ifndef GAUSSIAN_HELPERS_HPP
#define GAUSSIAN_HELPERS_HPP


// evaluate the pdf of point p
template <typename Point, typename NT>
NT eval_exp(Point &p, const NT &a) {
    return std::exp(-a * p.squared_length());
}


template <typename Point, typename NT>
NT get_max(Point &l, Point &u, const NT &a_i) {
    NT res;
    Point a = -1.0 * l, bef = u - l;
    Point b = (1.0 / std::sqrt((bef).squared_length())) * bef;
    Point z = (a.dot(b) * b) + l;
    NT low_bd = (l[0] - z[0]) / b[0], up_bd = (u[0] - z[0]) / b[0];
    if (low_bd * up_bd > 0) {
        //if(std::signbit(low_bd)==std::signbit(up_bd)){
        res = std::max(eval_exp(u, a_i), eval_exp(l, a_i));
    } else {
        res = eval_exp(z, a_i);
    }

    return res;
}


template <typename NT>
NT get_max_coord(const NT &l, const NT &u, const NT &a_i) {
    if (l < 0.0 && u > 0.0) {
        return 1.0;
    }
    return std::max(std::exp(-a_i * l * l), std::exp(-a_i * u * u));
}

#endif // GAUSSIAN_HELPERS_HPP
