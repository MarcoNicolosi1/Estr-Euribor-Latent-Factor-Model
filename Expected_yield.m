function [Ey,Er_T1] = Expected_yield(t,T1,T2, t_J, par1,par2,parJ,rho,lambda_xi,lambda_theta,r_t,X_t)
% Calcola il valore atteso E_t[y(T1,T2)]
% restituisce anche E_tr_T1
div=365;
k_xi=par1{1};  
sig_xi = par1{2};
aJ=parJ{1};
gJ=parJ{2};
wJ=parJ{3};

k_theta=par2{1};
theta_bar=par2{2};
sig_theta=par2{3};


Er_T1=E_tr_T();
[a, b] = YieldBond_Affine_2f(T1,T2, t_J, par1,par2,rho,parJ, Er_T1);
EX=E_tX_T(T1);
%Ey=a-b'*EX;
Ey=a+b'*EX; %MARCO

function E=E_tX_T(T)
tau=days(T-t)/div;
A_Obj = A_Obj_2f(tau,k_theta,k_xi,lambda_xi,lambda_theta,rho,sig_xi,sig_theta,theta_bar);
B_Obj = B_Obj_2f(tau,k_theta,k_xi);
E = A_Obj + B_Obj*X_t; 
end

function Er=E_tr_T()
% sommare tutti i salti...
% E_tr_T=r_t+sum(E_tJ_k)=r_t+sum(E_tX_k)
tempo_di_salto = t_J(and(t_J > t,t_J<=T1));
Er=r_t;
    %valore atteso del salto k-simo
    for k=1:length(tempo_di_salto)
    tau = days(tempo_di_salto(k) - t)/div;
    A_Obj = A_Obj_2f(tau,k_theta,k_xi,lambda_xi,lambda_theta,rho,sig_xi,sig_theta,theta_bar);
    B_Obj = B_Obj_2f(tau,k_theta,k_xi);
    E_tX_k=E_tX_T(tempo_di_salto(k));
    Er = Er+A_Obj(1) + B_Obj(1,:)*E_tX_k; 
    end

end

end






 