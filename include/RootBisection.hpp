/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include "Exception.hpp"

#ifndef ANPI_ROOT_BISECTION_HPP
#define ANPI_ROOT_BISECTION_HPP

namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], using the bisection method.
     *
     * @param funct a std::function of the form "T funct(T x)"
     * @param xl lower interval limit
     * @param xu upper interval limit
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootBisection(const std::function<T(T)> &funct, T xl, T xu, const T eps) {

        if (xl >= xu) {//check for correct interval
            throw anpi::Exception("inverted interval");
        }
        if (funct(xl) * funct(xu) > T(0)) {//check if both extremes have same sign
            throw anpi::Exception("unenclosed root");
        }

        int maxi = 10000;
        T xr = xl; // estimate the root as the lower limit
        T fl = funct(xl);
        T ea = T(); //relative error
        for (int i = 0; i < maxi; ++i) {
            T xrold(xr);
            xr = (xl + xu) / T(2);
            T fr = funct(xr);

            if (std::abs(xr) > std::numeric_limits<T>::epsilon()) {//check division for 0
                ea = std::abs((xrold - xr) / xr) * T(100); // calculate the approximate relative error
            }
            T cond = fl * fr;
            if (cond < T(0)) {//continuew with left side
                xu = xr;
            } else if (cond > T(0)) {//continuew with right side
                xl = xr;
                fl = fr;
            } else {// there was no error!
                ea = T(0);
                xr = (std::abs(fl) < std::numeric_limits<T>::epsilon()) ? xl : xr;
            }
            if (ea < eps) {
                return xr;
            }
        }
        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }
}

#endif

