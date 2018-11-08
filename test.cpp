#include <iostream>
#include<math.h>
#include<conio.h>
#include<stdio.h>
#include<fstream>

using namespace std;

main()
{
  double m_p, T_mean, mach_no, vel_disp, pdf[], n_H[], x[], lambda_jeans_cm[];
  double tau[], sigma_H2, n_LW[], n_H2[], X_H2[], K_b, G, c_s;
  float x_mean, n_H_range, n_H_mean, local_rad_field, Z;

  n_H_range = 100

  for(int i=1; i<=n_H_range; i++)
  {

  }
}


m_p = 1.672621777e-24  # g
  T_mean = 10.           #K
  mach_no = 5.
  n_H_range = 1000000
  n_H_mean = 100    # [H] cm^-3
  vel_disp = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #vel_disp in pdf
  x_mean = 1
  pdf = []
  n_H = []    # [H] cm^-3
  x = []
  lambda_jeans_cm = []     # cm
  tau = []    #optical depth
  sigma_H2 = 1000 * m_p       # assuming solar metallicity #cm^2
  n_LW = []     #number of Lyman-Werner photons
  local_rad_field = 1.0  #solar units
  n_H2 = []       #[H]/cc
  Z = 1. # solar metallicity
  X_H2 = []
