/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 04.08.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_RIDDER_HPP
#define ANPI_ROOT_RIDDER_HPP

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
    T rootRidder(const std::function<T(T)>& funct,T xi,T xii,const T eps) {

        auto maxit = std::numeric_limits<T>::digits;
        T fl = funct(xi);
        T fh = funct(xii);
        T xl = xi;
        T xh = xii;
        T ea=T();

        for (int j = 0; j < maxit; ++j){
            T xm = 0.5 * ( xl + xh );
            T fm = funct (xm);
            T s = sqrt( fm * fm -fl * fh);

            T xnew = xm + (xm-xl)*( ( fl>=fh? T(1) : T(-1) ) * fm/s );

            if(xm!=T(0)) {
                ea= std::abs((xnew-xm)/xm)*T(100);
            }
            if( ea < eps){
                return xnew;
            }
            T fnew = funct(xnew);
            if (fnew == T(0)){
                return xnew;
            }
            if( fnew*fm < T(0) ){
                xl = xm;
                fl = fm;
                xh =xnew;
                fh = fnew;
            }else if( fl*fnew  < T(0) ){
                xh = xnew;
                fh = fnew;
            }else if( fh*fnew < T(0) ){
                xl =xnew;
                fl =fnew;
            }

        }
        return std::numeric_limits<T>::quiet_NaN();

    }

}

#endif

