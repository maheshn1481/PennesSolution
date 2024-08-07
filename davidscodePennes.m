clc
clear
close all

u0=21;
k0= .57;
w0= 6.;
P= 10.;
r = .001;

temp = pennesmht(u0,k0,w0,P,r)

function temperature = pennesmht(u0,k,w,P,r)
      mua   = 0.45e+2
      mus   = 47.0e+2
      anfact= .9
      mutr  = mua+mus*(1.0-anfact)
      mueff = sqrt(3.0*mua*mutr)
      ua    =  310
      R1    =  .001
      R2    =  .03
      cblood = 3840.0
      s1 = 3.0/4.0/pi*P*mua*mutr/(w-k*mueff*mueff)*exp(-mueff*r)/r+ua;
      s2 = s1;
      s5 = 1/r*exp(sqrt(w/k)*r)*(-4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*u0*pi*R1*w+4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff+3.0*sqrt(w/k)*R2*P*mua*mutr*exp(-sqrt(w/k)*R2-mueff*R1)+4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*ua*pi*R1*w-4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff-3.0*P*mua*mutr*mueff*R2*exp(-mueff*R2-sqrt(w/k)*R1)-3.0*P*mua*mutr*exp(-mueff*R2-sqrt(w/k)*R1)+4.0*exp(-sqrt(w/k)*R2)*ua*pi*R1*w-4.0*exp(-sqrt(w/k)*R2)*u0*pi*R1*w+4.0*exp(-sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff+3.0*P*mua*mutr*exp(-sqrt(w/k)*R2-mueff*R1)-4.0*exp(-sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff)/4.0;     
      s6 = exp(-sqrt(w/k)*(-R1+R2))/(-w+k*mueff*mueff)/pi/(exp(-2.0*sqrt(w/k)*(-R1+R2))+sqrt(w/k)*R2*exp(-2.0*sqrt(w/k)*(-R1+R2))-1.0+sqrt(w/k)*R2);      
      s4 = s5*s6;      
      s6 = 1/r*exp(-sqrt(w/k)*r)*exp(-sqrt(w/k)*(-R1+R2))/4.0;     
      s9 = 4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*w-4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*w*sqrt(w/k)*R2-4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff+4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff*sqrt(w/k)*R2-3.0*P*mua*mutr*exp(sqrt(w/k)*R2-mueff*R1)+3.0*P*mua*mutr*sqrt(w/k)*R2*exp(sqrt(w/k)*R2-mueff*R1)-4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*w+4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*w*sqrt(w/k)*R2+4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff-4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff*sqrt(w/k)*R2+3.0*P*mua*mutr*mueff*R2*exp(sqrt(w/k)*R1-mueff*R2)+3.0*P*mua*mutr*exp(sqrt(w/k)*R1-mueff*R2);      
      s10 = 1/(exp(-2.0*sqrt(w/k)*(-R1+R2))+sqrt(w/k)*R2*exp(-2.0*sqrt(w/k)*(-R1+R2))-1.0+sqrt(w/k)*R2);      
      s8 = s9*s10;      
      s9 = 1/pi/(-w+k*mueff*mueff);      
      s7 = s8*s9;      
      s5 = s6*s7;      
      s3 = s4+s5;      
      temperature  = s2+s3;
end