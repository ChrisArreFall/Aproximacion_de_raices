/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <cstdlib>
#include <iostream>
#include "RootSecant.hpp"

double x2(double x){
    return std::abs(x)-std::exp(-x);
}

int main() {

  // Put your main code in here
  auto xr = anpi::rootSecant<double>(x2,0.0,1.0,0.00001);
  std::cout << xr << std::endl; // REMOVE-ME!
  
  return EXIT_FAILURE;
}
  
