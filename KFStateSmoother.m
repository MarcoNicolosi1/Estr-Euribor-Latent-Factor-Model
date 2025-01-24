function [mAs, mPs] = KFStateSmoother(y, Z, d, GG, T_mtx, c, HH, a, P)
%----------------------------------------------------------------------
 
m_states = length(a);
[n_obs,N]= size(y);   
y = y';

%----------------------------------------------------------------------
% Initialisation of arrays of matrices and scalars
av = NaN([N 1 n_obs]); 
aF = NaN([N N n_obs]);
aK = NaN([m_states N n_obs]);
aaf = NaN([m_states 1 n_obs]);
aP = NaN([m_states m_states n_obs]);
for i = 1:n_obs
	aaf(:,:,i) = a; 	
    aP(:,:,i) = P;	 % store a_{t|t-1}, P_{t|t-1}    
	v = y(:,i) - Z(:,:,i) * a - d(:,i); 		
    F = Z(:,:,i) * P * Z(:,:,i)' + GG;
    K = T_mtx(:,:,i) * P * Z(:,:,i)' / F;
    av(:,:,i) = v;             
    aF(:,:,i) = F;
	aK(:,:,i) = K;
    a = T_mtx(:,:,i) * a + c(:,i)+ K * v;		
    P = T_mtx(:,:,i) * P * T_mtx(:,:,i)' + HH(:,:,i) - K * F * K';
end

vr = zeros(m_states,1);  
mN = zeros(m_states,m_states);
mAs = NaN(m_states, n_obs); 	
mPs = NaN(m_states, n_obs); 

for i = n_obs:-1.0:1
 	Finv =  inv(aF(:,:,i));
	mL = T_mtx(:,:,i) - aK(:,:,i) * Z(:,:,i);			
	vr = Z(:,:,i)' * Finv * av(:,:,i) + mL' * vr;	
	mN = Z(:,:,i)' * Finv * Z(:,:,i) + mL' * mN * mL;
	a = aaf(:,:,i); P = aP(:,:,i);
    mAs(:,i) =  a + P * vr;
	mPs(:,i) = diag(P - P * mN * P);        
end
end