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

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {
    /**
     * Calculate the derivative of the function by means of centered differences
     * @tparam T
     * @param funct function to be derived
     * @param xi location to find the derivative
     * @return the derivative at xi
     */
    template<typename T>
    T derivative(const std::function<T(T)> &funct, T xi) {


      T h = T(.001);
      T fp = (funct(xi + h) - funct(xi - h)) / (T(2)*h) ; //diferencias centradas
      return fp;
    }

    /**
     * Find the roots of the function funct looking by means of the
     * Newton-Raphson method
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xi initial root guess
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootNewtonRaphson(const std::function<T(T)>& funct,T xi,const T eps) {

      T x2;
      T ea =T();
      int maxi=10000;
      for (int i  = 0; i< maxi; ++i) {
          x2 = xi - funct(xi) / derivative(funct, xi); //calculate next x

          if (std::abs(x2)  !=T(0)) { //check division for 0
              ea = std::abs((x2-xi)/x2);
          }

          if (funct(x2) == T(0) || (std::abs(ea) < eps) ) return x2; //check if precision is reached

          xi = x2;
      }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();

    }

}

#endif
