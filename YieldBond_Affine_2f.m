function [a, b] = YieldBond_Affine_2f(t,T, t_J, par1,par2,rho,parJ, r_t)

% P(t,T) = exp(-(T-t)*yield) = exp(alpha+beta*X_t)
% yield = - alpha/(T-t) - beta/(T-t)*X_t
% a=- alpha/(T-t);
% b=-beta/(T-t)
% yield = a + b*X_t 

% dr_t=J_t dD_t
% J_t \sim N(a_J+g_J*X_t,w_J^2)
% g_J=[1;0];
% E_t[exp(u'*X_T]=psi(T-t,u)+phi(T-t,u)'*X_t


% P(t,T) = E_t exp(-\int_t^T r_s ds) =
% = exp(-r_t*(T-t))*E_t \exp(-sum_k(J_k (T-T_k)))=



div=360;


k_xi=par1{1};  
sig_xi = par1{2};
aJ=parJ{1};
gJ=parJ{2};
wJ=parJ{3};

k_theta=par2{1};
theta_bar=par2{2};
sig_theta=par2{3};



t_J = t_J((t <= t_J)&(t_J < T)); %a column vector

if isempty(t_J)
    a = r_t;
    b= [0;0];
else
    len_tJ=length(t_J); 
    T_tJ = days(T-t_J)./div;% maturities in years act/360
    
    diff_tJ=days(diff(t_J))./div;
    
    
    N_fact=2;
    
    F=zeros(N_fact,len_tJ);
    F(:,end)=-gJ*(T_tJ(end));
    G=zeros(1,len_tJ);
    G(end)=-aJ*(T_tJ(end))+0.5*wJ^2*T_tJ(end)^2;
    for k=len_tJ-1:-1:1
        %[phi,psi]=affine_1f(diff_tJ(k+1),bar_xi,k_xi,sig_xi,F(k+1));
        [phi,psi]=affine_2f(diff_tJ(k),k_theta,k_xi,rho,sig_xi,sig_theta,theta_bar,F(:,k+1));
        F(:,k)=psi-gJ*(T_tJ(k));
        G(k)=-aJ*(T_tJ(k))+0.5*wJ^2*T_tJ(k)^2+phi+G(k+1);
    end
    % ultimo passo

    
    [phi0,psi0]=affine_2f(days(t_J(1)-t)/div,k_theta,k_xi,rho,sig_xi,sig_theta,theta_bar,F(:,1));
    mat=days(T-t)/div;
    alpha=-r_t*mat+G(1)+phi0;
    beta=psi0;
    a=- alpha/mat;
    b=-beta/mat;
end






    function [A,B,V]=fun_ABV(tau)
        % formule 37-39-39 pag. 34 Schlogl et al.

        B=exp(-tau*K); 
        A=(1-B)*x_bar;
        V=sig_x^2/(2*K)*(1-exp(-2*K*tau));
    end

end
