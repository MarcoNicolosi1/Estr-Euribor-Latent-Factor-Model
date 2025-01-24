function [Z, d, GG, T_mtx, c, HH, a, P] = SetStateSpaceModel_4F(hyper_par, ESTR, t, T, t_J, N_Estr, N_Euribor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_obs, N] = size(T); %maturities, N = number of series for estimation
m_states = 4;     % number of states
g = N+m_states;   % number of errors


% Measurement equation
% y_{t} = Z_t \alpha_{t} +  d_t + G_t \epsilon_t

% Transition equation
% \alpha_{t+1} = T_mtx_t \alpha_{t} + c_t+ H_t \epsilon_t

mZ_Euribor = nan(N_Euribor,m_states,n_obs);
vd_Euribor = nan(N_Euribor,m_states,n_obs);
mZ_Ester = nan(N_Estr,m_states,n_obs);
vd_Ester = nan(N_Estr,n_obs);

Z = nan(N,m_states,n_obs);
d = nan(N,n_obs);

%a = nan(n,N);
%b = nan(n,N);
T_mtx = nan(m_states,m_states,n_obs);
c = nan(m_states,n_obs);
HH = nan(m_states,m_states,n_obs);



% gJ=parJ{2};
% wJ=parJ{3};
% 

k_xi = exp(hyper_par(1));
sig_xi = exp(hyper_par(2));
k_theta = exp(hyper_par(3));
theta_bar = hyper_par(4);
sig_theta = exp(hyper_par(5));




par1{1} = k_xi;
par1{2} = sig_xi;
par2{1} = k_theta;
par2{2} = theta_bar;
par2{3} = sig_theta;
parJ{1} = 0;%aJ
parJ{2} = [1; 0];% gJ;
parJ{3} = exp(hyper_par(6));% w_J % std dev of jumps
w = parJ{3}; 
rho = 0.9*(1 - 2*exp(hyper_par(7))/(1+exp(hyper_par(7))));% correlazione tra le due variabili di stato xi e theta


k_spread = exp(hyper_par(8));
sig_spread = exp(hyper_par(9));
spread_bar = hyper_par(10);


sigma = exp(hyper_par(11));  % std dev measurement error 

if length(hyper_par) == 11
    lambda_xi = 0; % market price of risk of xi
    lambda_theta = 0; % market price of risk of theta
    gamma_J = 0;
    lambda_spread = 0;
else
    lambda_xi = exp(hyper_par(12)); % market price of risk of xi
    lambda_theta = exp(hyper_par(13)); % market price of risk of theta
    gamma_J = exp(hyper_par(14)); % market price of risk for jumps
    lambda_spread = exp(hyper_par(15)); % market price of risk of spread
end

div = 360;

for i = 1:n_obs
    for j = 1:N_Estr
      %  target_t = ESTR(i);
        target_t = 0;        
        [a_Estr, b_Estr] = YieldBond_Affine_2f(t(i),T(i,j), t_J, par1,par2,rho,parJ, target_t);
        mZ_Ester(j,:,i) = [b_Estr' 1 0];
        vd_Ester(j,i) = a_Estr;
    end
    for j = 1:N_Euribor
        tau = days(T(i,N_Estr+j)-t(i))/div;
        r0 = 0;
        [a_spread,b_spread]=vas_ytm(k_spread,spread_bar,sig_spread,r0,tau);
        target_t = 0;
        [a_Estr, b_Estr] = YieldBond_Affine_2f(t(i),T(i,N_Estr+j), t_J, par1,par2,rho,parJ, target_t);
        mZ_Euribor(j,:,i) = [b_Estr' 1 b_spread];
        vd_Euribor(j,i) = a_Estr+a_spread;
    end
    Z(:,:,i) = [mZ_Ester(:,:,i); mZ_Euribor(:,:,i)];
    d(:,i) = [vd_Ester(:,i); vd_Euribor(:,i)];

    if i == n_obs
        dt = days(t(n_obs)-t(n_obs-1))/div;
    else
        dt = days(t(i+1)-t(i))/div;
    end
    
    X_bar_Obj = -[-1/k_xi -1/k_theta; 0 -1/k_theta]*[sig_xi*lambda_xi; k_theta*theta_bar + sig_theta*(rho*lambda_xi + sqrt(1-rho^2)*lambda_theta)];
   
    A_Obj = A_Obj_2f(dt,k_theta,k_xi,lambda_xi,lambda_theta,rho,sig_xi,sig_theta,theta_bar);
   
    B_Obj = B_Obj_2f(dt,k_theta,k_xi);
    V_Obj = V_Obj_2f(dt,k_theta,k_xi,rho,sig_xi,sig_theta);

    d_J = sum(t(i) == t_J);

    B_ester = [B_Obj, [0; 0]; [gamma_J*d_J 0 1]];    
    T_mtx(:,:,i)  = [B_ester, zeros(m_states-1,1); zeros(1,m_states-1), exp(-k_spread*dt)];    
        
    spread_bar_Obj = spread_bar +sig_spread*lambda_spread/k_spread;
    c(:,i)  = [A_Obj; 0; spread_bar_Obj*(1-exp(-k_spread*dt))];

    cholV = chol(V_Obj,"lower");
    Var_ester = [zeros(m_states-1,N), [cholV, [0;0]; [0 0 w*d_J]] ]*[zeros(m_states-1,N), [cholV, [0;0]; [0 0 w*d_J]]]'; % H_t * H_t'  
    HH(:,:,i) = [Var_ester, zeros(m_states-1,1); zeros(1,m_states-1), sig_spread^2/2/k_spread*(1-exp(-2*k_spread*dt))];  

end

GG = [sigma*eye(N,N) zeros(N,m_states)]*[sigma*eye(N,N) zeros(N,m_states)]'; % G_t * G_t' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial conditions 
a = [X_bar_Obj; ESTR(n_obs); spread_bar_Obj];     % a_{1|0}

tau = 1000;
mP_old = [V_Obj_2f(tau,k_theta,k_xi,rho,sig_xi,sig_theta) [0; 0]; [0 0 w^2]];
P = [mP_old zeros(m_states-1,1);zeros(1,m_states-1),sig_spread^2/2/k_spread]; % P_{1|0}     


 

end


