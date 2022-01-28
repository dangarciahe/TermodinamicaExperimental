#include <iostream>
#include <cmath>

double integrando1(double x); //x^3 / e^x-1
double integrando2(double theta, double T);
double D(double x);
double C_V (double theta, double T);
double integral1(double a, double b, int n_pasos); //Integral de Riemann
double integral2(double a, double b, int n_pasos, double theta);
double biseccion (double tol, double theta_0, double theta_f); //Metodo de la biseccion  
double termino(double theta); //Pasar el delta Q al lado derecho, y esto es igual a ese monstruo


//Inputs
const double m = 6.9e-3; //Cambio de masa 
const double masa = 25.2e-3; //Masa de la muestra
const double molecular = 63.546e-3; //Masa molar de la muestra


//Constantes
const double N = 6.023e23; //Numero de Avogadro
const double k_B = 1.380649e-23; //k Boltzmann
const double R = 8.314462; //Constante de los gases ideales
const double L = 266.77e3; //Calor latente
const double Q = m*L; // Calor
const double T_i = 292.15; //Temperatura inicial
const double T_f = 77; //Temperatura final
const int n_steps = 4000; //N de la integral



int main(){
  std::cout.precision(4);
  std::cout << "La temperatura de Debye es:  "
	    << biseccion (1e-1, 1, 1000) << " K \n"; 
  return 0;
}



double integrando1(double x)
{
  double valor = (x*x*x*x*std::exp(x))/((std::exp(x)-0.99999999)*(std::exp(x)-0.99999999));
  return valor;
}

double integral1(double a, double b, int n_pasos)
{
  double delta = (b-a)/n_pasos;
  double suma = 0;
  double x = a;
  for (int ii = 0; ii<=n_pasos; ii++)
    {
      double y = integrando1(x);
      suma += y;
      x+= delta;
    }
  double valor = suma*delta; 
  return valor;
}


double D(double x){  
  double valor;
  valor = (3/(x*x*x))*integral1(0, x, n_steps);
  return valor;
}


double C_V (double theta, double T){
  double valor = 3*R*D(theta/T);
  return valor;
}



double integrando2 (double theta, double T){
  double valor = (masa/molecular)*C_V(theta, T);
  return valor;
}


double integral2(double a, double b, int n_pasos, double theta){
  double delta = (b-a)/n_pasos;
  double suma = 0;
  double x = a;
  for (int ii = 0; ii<=n_pasos; ii++)
    {
      double y = integrando2(theta, x);
      suma += y;
      x+= delta;
    }
  double valor = suma*delta; 
  return valor;
}


double termino(double theta){
  double valor = m*L - integral2(T_f, T_i, n_steps, theta);
  return valor;
}



double biseccion (double tol, double theta_0, double theta_f)
{
  double xmedio = (theta_0+theta_f)/2;
  double y0 = termino (theta_0);
  double yf = termino (theta_f);
  double y = termino (xmedio);
  if (y0*yf>0)
    {
      std::cout << "Errooooooooooooooo";
    }
  
  while (tol <= std::fabs(y)){
    if (y*yf > 0){
      theta_f = xmedio;
    }
    if (y*y0 > 0){
      theta_0 = xmedio;
    }
    
    xmedio = (theta_0+theta_f)/2;
    y0 = termino (theta_0);
    yf = termino (theta_f);
    y = termino (xmedio);
  }
  
  return xmedio;
}

