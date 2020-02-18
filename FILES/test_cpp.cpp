#include<iostream>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<conio.h>

using namespace std;

class H2_fraction
{
public:
  const float pi = 3.14159265
  float make_pdf(float, float, float); //if this shows an error, replace variables in paranthesis with just 'float'
  float make_lambda_jeans(float, float);
  float make_X_H2(float, float, float);
  float make_integral1();
  float varying_M();
  float varying_Z();
  float varying_G_o();
};

float H2_fraction :: make_pdf(float s, float s_bar, float sigma_s)
{
  float pdf;
  pdf = (1/np.sqrt(2*pi*pow(sigma_s, 2))) * (exp(-0.5*pow(((s - s_bar)/sigma_s), 2)));
  return (pdf);
}

float H2_fraction :: make_lambda_jeans(float n_H, float s)
{
  float m_p, K_b, G, lambda_jeans;
  int T_mean;
  m_p = 1.672621777*pow(10, -24);   // g
  T_mean = 10;                      //K
  K_b = 1.38064852*pow(10, -16);    // ergs K-1
  G = 6.67408*pow(10, -8);          // dyne cm^2 g^-2
  lambda_jeans = ((sqrt(K_b * T_mean / m_p)) / sqrt(4 * pi * G * n_H * m_p));
  return (lambda_jeans);
}

float H2_fraction :: make_X_H2(float n_H, float Z, float G_o)
{
  float m_p, K_b, G, kappa, DC, CC, rad_field_outside, exp_tau, numerator, denominator, X_H2;
  int T_mean;
  m_p = 1.672621777*pow(10, -24);   // g
  T_mean = 10;                      //K
  K_b = 1.38064852*pow(10, -16);    // ergs K-1
  G = 6.67408*pow(10, -8);          // dyne cm^2 g^-2
  kappa = 1000 * m_p;
  DC = 1.7*pow(10,-11);             //Destruction coefficient
  CC = 2.5*pow(10, -17);            //cm3 s-1  //Construction coefficient
  rad_field_outside = G_o           //in solar units
  exp_tau = exp(-kappa * n_H * ((sqrt(K_b * T_mean / m_p)) / sqrt(4* pi * G * n_H * m_p)));
  numerator = DC * rad_field_outside * exp_tau;
  denominator = CC * Z * n_H;
  X_H2 = 1 / (2 + (numerator/denominator) );
  return (X_H2);
}

float H2_fraction :: make_integral1()
{
  float s[1000], pdf[1000], n_H[1000], lambda_jeans[1000], X_H2[1000], X_H2_bar[1000];
  float sigma_s, smin, s_bar, ds, integral1, n_H_mean, Z, G_o;
  int mach_no;
  sigma_s = sqrt(log(1 + pow((0.3 * mach_no), 2)));
  s_bar = -0.5*pow(sigma_s, 2);
  smin = -4*sigma_s + s_bar;
  smax = 4*sigma_s + s_bar;
  ds = (smax - smin)/1000;
  n_H_mean = pow(10, 4);
  Z = 1;
  G_o = 1;
  for(int i=0; i<1000; i++)
  {
    s[i] = smin + i*ds
    n_H[i] = n_H_mean * np.exp(s)
    lambda_jeans[i] = make_lambda_jeans(n_H[i], s[i])
    X_H2[i] = make_X_H2(n_H[i], Z, G_o)
    pdf[i] = make_pdf(s[i], s_bar, sigma_s)
    integral1 += np.exp(s[i]) * pdf[i] * ds   //this should be ~1
  }
  //plotting(n_H, pdf, lambda_jeans, X_H2);
  return (integral1);
}

float H2_fraction ::make_integral()
{
  float s[1000], pdf[1000], n_H[1000], lambda_jeans[1000], 
}

main()
{
  // order of variables:
  // s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, integral1, X_H2_bar
  g1, g2, g3, g4, g5, g6, g7, g8, g9, g10 = varying_G_o()  //varying G_o
  m1, m2, m3, m4, m5, m6, m7, m8, m9, m10 = varying_M()  //varying mach_no
  z1, z2, z3, z4, z5, z6, z7, z8, z9, z10 = varying_Z()  //varying matallicity

}
