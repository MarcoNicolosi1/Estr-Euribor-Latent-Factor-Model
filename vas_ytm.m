function [A,B,y,P]=vas_ytm(a,b,sigma,r0,tau)
% [P,y]=vas_zc(a,b,sigma,r,tau)
%  Price and yield of a zero coupon bond in the Vasicek model
% a,b,sigma, parameters of the model.
% r0 is the initial short rate
% tau is the time to maturity

B=(1-exp(-a.*tau))/a;
A=((B-tau).*(a^2*b-sigma^2/2)/a^2)-sigma^2.*B.^2/(4*a);
P=exp(A-B.*r0);
y=(B.*r0-A)./tau;

A = -A/tau;
B = B/tau;