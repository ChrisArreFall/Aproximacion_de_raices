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

#ifndef ANPI_ROOT_INTERPOLATION_HPP
#define ANPI_ROOT_INTERPOLATION_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], by means of the interpolation method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xl lower interval limit
   * @param xu upper interval limit
   *
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */

  template<typename T>
  T rootInterpolation(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

    //interpolation
    T error = std::abs((xu - xl) / xu);
    T y_l = funct(xl);
    T y_u = funct(xu);
    T f_root;
    T raiz = xu - ((y_u * (xl - xu)) / (y_l - y_u));
    T raizProx = 0;


    while (error > eps) {
      y_l = funct(xl);
      f_root = funct(raiz);

      if (y_l * f_root < 0) {
        xu = raiz;
        raizProx  = xu - ((y_u * (xl - xu)) / (y_l - y_u));
        error = std::abs((raizProx - raiz) / raizProx);
        raiz = raizProx;
      }
      else if (y_l * f_root > 0) {
        xl = raiz;
        raizProx  = xu - ((y_u * (xl - xu)) / (y_l - y_u));
        error = std::abs((raizProx - raiz) / raizProx);
        raiz = raizProx;
      }
      else{
        break;
      }
    }

    //if y(raiz) is really different than 0, it did not find the root
    if (abs(funct(raiz)) > 0.5) {
      return std::numeric_limits<T>::quiet_NaN();
    }
    else {
      return raiz;
    }
  }
}


  
#endif

