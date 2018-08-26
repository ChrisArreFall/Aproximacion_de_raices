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
  T rootBisection(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

    //bisection
    if (funct(xl) * funct(xu) >= 0)
    {
      return std::numeric_limits<T>::quiet_NaN();
    }

    T c = xl;
    while ((xu-xl) >= eps)
    {
      // Find middle point
      c = (xl+xu)/2;

      // Check if middle point is root
      if (funct(c) == 0.0)
          break;

          // Decide the side to repeat the steps
      else if (funct(c)*funct(xl) < 0)
          xu = c;
      else
          xl = c;
    }
    return c;

  }

}
  
#endif

