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
    //newton
    template<typename T>
    T derivate(const std::function<T(T)>& funct, T xi) {


      T h = T(.01);
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

      T x2 = T();
      T ea =T();

      int maxi=10000;
      for (int i  = 0; i< maxi; ++i){
        x2 = xi - funct(xi)/derivate(funct, xi);
        if(x2 != T(0)) {
          ea = std::abs((x2 - xi) / x2) * T(100);
        }else{
          ea = std::abs(x2-xi);
        }
        if(ea < eps){
          if(funct(x2)<=std::numeric_limits<T>::epsilon()){
            return x2;
          }
        }
        xi = x2;
      }

      return std::numeric_limits<T>::quiet_NaN();
    }

}

#endif
