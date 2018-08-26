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
  T rootSecant(const std::function<T(T)>& funct,T xi,T xii,const T eps) {

    auto maxi = 10000; // max number of iterations
    T ea=T(); //delta of the change in X
    T tmp; // temp value to store the current xii
    for (int i = 0; i < maxi; ++i) { //loop until the max iterations
      tmp = xii; //store the current xi
      xii = xii - funct(xii) * (xi - xii) / (funct(xi) - funct(xii)); //calculate next x
      xi = tmp; //set the x i-1 as the last xii
      if(std::abs(xii)>std::numeric_limits<T>::epsilon()) {
        ea = std::abs((xii - xi) / xii) * T(100); //calculate delta
      }
      if (ea < eps){ //if the delta is less than the epsilon
        // we check that the X we get is actually a root
        return std::abs(funct(xii))<=eps? xii : std::numeric_limits<T>::quiet_NaN();
      }
    }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

