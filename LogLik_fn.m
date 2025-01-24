function LogLik = LogLik_fn(hyper_par, y, ESTR, t, T, t_J, model,N_Estr,N_Euribor)


    switch model            
        case '4F'
        [Z, d, GG, T_mtx, c, HH, a, P] = SetStateSpaceModel_4F(hyper_par, ESTR, t, T, t_J,N_Estr,N_Euribor);
    end

    LogF = 0; SumSquares = 0;     
    [n_obs,~]= size(y);   

    y = y';
    for i = 1:n_obs
	    v = y(:,i) - Z(:,:,i) * a - d(:,i); 		
        F = Z(:,:,i) * P * Z(:,:,i)' + GG;
        K = T_mtx(:,:,i) * P * Z(:,:,i)' / F;
	    a = T_mtx(:,:,i) * a + c(:,i)+ K * v;		    
        P = T_mtx(:,:,i) * P * T_mtx(:,:,i)' + HH(:,:,i) - K * F * K';
        LogF = LogF + log(det(F));   
        SumSquares = SumSquares + v'*(F\v);
    end
    LogLik = 0.5 * (LogF + SumSquares ); 



end
 