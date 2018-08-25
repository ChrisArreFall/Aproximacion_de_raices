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

#ifndef ANPI_ROOT_SECANT_HPP
#define ANPI_ROOT_SECANT_HPP

namespace anpi {

    /**
     * Find a root of the function funct looking for it starting at xi
     * by means of the secant method.
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xi initial position
     * @param xii second initial position
     *
     * @return root found, or NaN if no root could be found
     */
    template<typename T>
    T rootSecant(const std::function<T(T)> &funct, T xi, T xii, const T eps) {

        int maxi = 10000;
        T dx;
        T tmp;
        for (int i = 0; i < maxi; ++i) {
            tmp = xii;
            xii = xii - funct(xii) * (xi - xii) / (funct(xi) - funct(xii));
            xi = tmp;
            dx = std::abs(xii-xi);
            if (dx < eps){
                return std::abs(funct(xii))<=eps? xii : std::numeric_limits<T>::quiet_NaN();
            }
        }
        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }

}

#endif

