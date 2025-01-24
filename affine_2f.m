function [phi_2f,psi_2f] = affine_2f(T,k_theta,k_xi,rho,sig_xi,sig_theta,theta_bar,u)
% [phi_2f,psi_2f] = affine_2f(T,k_theta,k_xi,rho,sig_xi,sig_theta,theta_bar,u)
%   Detailed explanation goes here
u_1=u(1);
u_2=u(2);

t2 = T.*k_theta;
t3 = T.*k_xi;
t4 = k_theta+k_xi;
t5 = k_xi.^2;
t6 = sig_theta.^2;
t7 = u_1.^2;
t10 = -k_xi;
t11 = 1.0./k_theta;
t8 = t2.*2.0;
t9 = t3.*2.0;
t12 = -t2;
t14 = -t3;
t19 = k_theta+t10;
t20 = 1.0./t4;
t13 = -t8;
t15 = -t9;
t16 = exp(t12);
t24 = 1.0./t19;
t26 = t12+t14;
t17 = exp(t13);
t18 = exp(t15);
t21 = t16-1.0;
t25 = t24.^2;
t27 = exp(t26);
t22 = t17-1.0;
t23 = t18-1.0;
t28 = t27-1.0;
phi_2f = -t21.*theta_bar.*u_2-(t6.*t11.*t22.*u_2.^2)./4.0-(sig_xi.^2.*t7.*t23)./(k_xi.*4.0)-k_theta.*t24.*theta_bar.*u_1.*(exp(t14)-1.0)-(k_xi.*t6.*t7.*t23.*t25)./4.0+k_xi.*t21.*t24.*theta_bar.*u_1-(rho.*sig_xi.*sig_theta.*t7.*t23.*t24)./2.0-(t5.*t6.*t7.*t11.*t22.*t25)./4.0+t5.*t6.*t7.*t20.*t25.*t28+(k_xi.*t6.*t11.*t22.*t24.*u_1.*u_2)./2.0-rho.*sig_xi.*sig_theta.*t20.*t28.*u_1.*u_2+t6.*t10.*t20.*t24.*t28.*u_1.*u_2+k_xi.*rho.*sig_xi.*sig_theta.*t7.*t20.*t24.*t28;

t2 = T.*k_theta;
t3 = T.*k_xi;
t4 = -t2;
t5 = -t3;
t6 = exp(t4);
t7 = exp(t5);
psi_2f = [t7.*u_1;t6.*u_2-(u_1.*(k_xi.*t6-k_xi.*t7))./(k_theta-k_xi)];
end