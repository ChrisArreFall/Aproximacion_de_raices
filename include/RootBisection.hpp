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

    T intentos = 10000000;
    T i = 0;
    T medio = 0;
    T y = 1;
    T raiz = 1;
    while((abs(y) > eps) && i < intentos ){
      medio = (xl+xu)/2;
      y = funct(medio);
      if(y>0){
        xu = medio;
      }
      else if(y<0){
        xl = medio;
      }
      else{
        raiz = medio;
        return raiz;
      }
      i++;
    }
    if(intentos == i){
        return std::numeric_limits<T>::quiet_NaN();
    }
    else{
        raiz = medio;
        return raiz;
    }

  }

}
  
#endif

