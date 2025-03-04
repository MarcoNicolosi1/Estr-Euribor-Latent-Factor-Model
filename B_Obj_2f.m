function B_Obj_2f = B_Obj_2f(T,k_theta,k_xi)
%B_Obj_2f
%    B_Obj_2f = B_Obj_2f(T,K_THETA,K_XI)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    27-Jul-2024 16:25:17

t2 = T.*k_theta;
t3 = T.*k_xi;
t4 = -t2;
t5 = -t3;
t6 = exp(t4);
t7 = exp(t5);
B_Obj_2f = reshape([t7,0.0,-(k_xi.*t6-k_xi.*t7)./(k_theta-k_xi),t6],[2,2]);
end
