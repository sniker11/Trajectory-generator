#include <umrob/polynomial.h>
#include <math.h>
#include <iostream>

namespace umrob {

Polynomial::Constraints& Polynomial::constraints() {
    return constraints_;
}

const Polynomial::Constraints& Polynomial::constraints() const {
    return constraints_;
}


// clang-format off
/* Simplified coefficient for xi = 0 and dx = xf-xi
 * a = -(12*yi - 12*yf + 6*dx*dyf + 6*dx*dyi - d2yf*dx^2 + d2yi*dx^2)/(2*dx^5)
 * b = (30*yi - 30*yf + 14*dx*dyf + 16*dx*dyi - 2*d2yf*dx^2 + 3*d2yi*dx^2)/(2*dx^4) 
 * c = -(20*yi - 20*yf + 8*dx*dyf + 12*dx*dyi - d2yf*dx^2 + 3*d2yi*dx^2)/(2*dx^3) 
 * d = d2yi/2 Important avant de vous connecter
Par mesure de sécurité, votre ordinateur doit être à jour tant au niveau système, qu'antivirus. De plus, vous devez avoir pris connaissance de la charte informatique de l'Université de Montpellier et de RENATER dont vous 
 * e = dyi 
 * f = yi
 */
void Polynomial::computeCoefficients() {

   /* double dx=polyY.xf-polyY.xi;
    polyY.xi=0;
coef.a= -(12*polyY.yi - 12*polyY.yf+6*dx*polyY.dyf+6*dx*polyY.dyi-polyY.d2yf*dx*dx+polyY.d2yi*dx*dx)/(2*pow(dx,5));
coef.b= (30*polyY.yi - 30*polyY.yf+14*dx*polyY.dyf+16*dx*polyY.dyi-2*polyY.d2yf*dx*dx+3*polyY.d2yi*dx*dx)/(2*pow(dx,4));
coef.c= -(20*polyY.yi - 20*polyY.yf+8*dx*polyY.dyf+12*dx*polyY.dyi-polyY.d2yf*dx*dx+3*polyY.d2yi*dx*dx)/(2*pow(dx,3));
coef.d=polyY.d2yi/2
coef.e=polyY.dyi;
coef.f=polyY.yi;*/
double dx=constraints_.xf-constraints_.xi;
    constraints_.xi=0;
coefficients_.a= -(12*constraints_.yi - 12*constraints_.yf+6*dx*constraints_.dyf+6*dx*constraints_.dyi-constraints_.d2yf*dx*dx+constraints_.d2yi*dx*dx)/(2*pow(dx,5));
coefficients_.b= (30*constraints_.yi - 30*constraints_.yf+14*dx*constraints_.dyf+16*dx*constraints_.dyi-2*constraints_.d2yf*dx*dx+3*constraints_.d2yi*dx*dx)/(2*pow(dx,4));
coefficients_.c= -(20*constraints_.yi - 20*constraints_.yf+8*dx*constraints_.dyf+12*dx*constraints_.dyi-constraints_.d2yf*dx*dx+3*constraints_.d2yi*dx*dx)/(2*pow(dx,3));
coefficients_.d=constraints_.d2yi/2;
coefficients_.e=constraints_.dyi;
coefficients_.f=constraints_.yi;
}
// clang-format on

//! y = ax^5 + bx^4 + cx^3 + dx^2 + ex + f
double Polynomial::evaluate(double x) {

double y;
y=coefficients_.a*pow(x,5)+coefficients_.b*pow(x,4)+coefficients_.c*pow(x,3)+coefficients_.d*pow(x,2)+coefficients_.e*x+coefficients_.f;
if (constraints_.xf> constraints_.xi)
{if(x>constraints_.xf)
  y=constraints_.yf;  
 if(x<constraints_.xi)
 y=constraints_.yi;}
else
{if (constraints_.xf< constraints_.xi)
 {if(x<constraints_.xf)
  y=constraints_.yf;  
 if(x>constraints_.xi)
 y=constraints_.yi;}
 else
 std::cout<<"xi=xf";
 
}
    return y;
}
//! dy = 5ax^4 + 4bx^3 + 3cx^2 + 2dx + e
double Polynomial::evaluateFirstDerivative(double x) {
double dy;
 dy=5*coefficients_.a*pow(x,4)+4*coefficients_.b*pow(x,3)+3*coefficients_.c*pow(x,2)+2*coefficients_.d*x+coefficients_.e;
 if (constraints_.xf> constraints_.xi)
{
 if(x>constraints_.xf)
  dy=constraints_.dyf;  

 if(x<constraints_.xi)
 dy=constraints_.dyi;
}
else
{if (constraints_.xf< constraints_.xi)
 {
     if(x<constraints_.xf)
  dy=constraints_.dyf;  

 if(x>constraints_.xi)
 dy=constraints_.dyi;
 }
 else
 {
    std::cout<<"xi=xf";
 }
}
    return dy;
}

//! d2y = 20ax^3 + 12bx^2 + 6cx + 2d
double Polynomial::evaluateSecondDerivative(double x) {
    double d2y;
 d2y=20*coefficients_.a*pow(x,3)+12*coefficients_.b*pow(x,2)+12*coefficients_.c*x+2*coefficients_.d;
 if (constraints_.xf> constraints_.xi)
{
 if(x>constraints_.xf)
  d2y=constraints_.d2yf;  

 if(x<constraints_.xi)
 d2y=constraints_.d2yi;
}
else
{if (constraints_.xf< constraints_.xi)
 {
     if(x<constraints_.xf)
  d2y=constraints_.d2yf;  

 if(x>constraints_.xi)
 d2y=constraints_.d2yi;
 }
 else
 {
    std::cout<<"xi=xf";
 }
}
    return d2y;
}


} // namespace umrob