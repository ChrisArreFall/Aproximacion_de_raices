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
  T derivate(const std::function<T(T)>& funct, T xi, const T eps) {

    T aprox = 0;
    T anterior = 0;
    T error = 0;
    T h = 10^(-15);
    T lambda = 0.5;


    do{

      h = h/lambda;
      anterior = aprox;
      aprox = (funct(xi) - funct(xi - h)  ) / h ;
      error = abs( (aprox - anterior) / aprox );

    }while(error > eps);

    return aprox;
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

    T x2 = 0;
    T f = 0;
    T error = 1;
    T i = 0;

    while(error > eps){

      x2 = xi - funct(xi)/derivate(funct, xi, eps);
      f = funct(x2);
      error = abs( (f - funct(xi)) / f );
      xi = x2;
      ++i;

      if(i == 10000){
        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
      }
    }

    return f;

  }

}
  
#endif
