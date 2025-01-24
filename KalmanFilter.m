function [Innovations, VarInnovations, StatePred, CovStatePred,LogLik] = ...
    KalmanFilter(y, Z, d, GG, T_mtx, c, HH, a, P)
%----------------------------------------------------------------------
m_states = length(a); % n. of state elements 
[n_obs,N]= size(y);   

y = y';
%----------------------------------------------------------------------
% Initialisation of matrices and scalars
LogF = 0; SumSquares = 0;     
Innovations = NaN(N, n_obs);  % stores the KF innovations
VarInnovations = NaN(N,N, n_obs); 
StatePred = NaN(m_states, n_obs); % stores the states predictions
CovStatePred = NaN(m_states, m_states, n_obs);

for i = 1:n_obs    
	v = y(:,i) - Z(:,:,i) * a - d(:,i); 		
    F = Z(:,:,i) * P * Z(:,:,i)' + GG;                                            
    vK = T_mtx(:,:,i) * P * Z(:,:,i)' / F;
	a = T_mtx(:,:,i) * a + c(:,i)+ vK * v;		
    P = T_mtx(:,:,i) * P * T_mtx(:,:,i)' + HH(:,:,i) - vK * F * vK';
    LogF = LogF + log(det(F));   
    SumSquares = SumSquares + v'*(F\v);

    % store quantities
    Innovations(:,i) = v;       
    VarInnovations(:,:,i) = F;
    StatePred(:,i) = a;       
    CovStatePred(:,:,i) = P;

end 
LogLik = -0.5 * (LogF + SumSquares );
end
 