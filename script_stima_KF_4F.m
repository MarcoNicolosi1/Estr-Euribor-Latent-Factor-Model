clear all 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to estimate a 4-dimensional latent factor model that simultaneously 
% captures the dynamics of the Estr OIS term structure and the EURIBOR term structure.
%
% Let y_Estr(t,T) be the yield to maturity T observed at time t derived by the
% term structure of Estr OIS.
%
% y_Estr(t,T) = -1/(T-t)ln(P_Estr(t,T)) where
% P_Estr(t,T) = E_t exp(-\int_t^T r_s ds)
%
% r_t is a short latent rate that is constant up to the next BCE meeting
% date when it can jump. The jump size is normal with expectation \xi_t and
% variance w_J^2. 
% \xi_t is a latent factor that follows an Ornstein Oulenbeck process whose 
% long term drift \theta is another latent factor, following an Ornstein Oulenbeck
% process.
%
% r_t = r_0 + \sum_{i \in times of jumps} J_i
% J_t \sim N(\xi_t,w_J^2)
% d\xi_t = k_x*(\theta_t-\xi_t)*dt + sig_xi*dW1_t
% d\theta_t = k_x*(bar_theta-\theta_t)*dt + sig_theta*dW2_t
%
% The model is affine. Hence
% y_Estr = a + b*[\xi_t;\theta_t;r_t] 
%
% Let s_t be a latent spread process following a Vasicek dynamics.
%
% r_t+s_t is used to model the term structure  of yields implied by Euribor
%
% The model is estimated via maximum Log likelihood and Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load DATA.mat

% Convert tables to timetable for data synchronization
ESTR = table2timetable(ESTR);
OISZC = table2timetable(OISZC);
EURZ = table2timetable(EURZ);

% Synchronize and remove missing values
OISZC = synchronize(ESTR, OISZC);          
EstrEuribZC = synchronize(OISZC,EURZ);
EstrEuribZC = rmmissing(EstrEuribZC);

% Extract and normalize interest rates
EstrZY=[EstrEuribZC.EUESTONZY,EstrEuribZC.EUESTTNZY,EstrEuribZC.EUEST1WZY,EstrEuribZC.EUEST1MZY,EstrEuribZC.EUEST3MZY,EstrEuribZC.EUEST6MZY,EstrEuribZC.EUEST9MZY,...
    EstrEuribZC.EUEST1YZY,EstrEuribZC.EUEST2YZY,EstrEuribZC.EUEST3YZY,EstrEuribZC.EUEST4YZY,EstrEuribZC.EUEST5YZY,EstrEuribZC.EUEST6YZY,...
    EstrEuribZC.EUEST7YZY,EstrEuribZC.EUEST8YZY,EstrEuribZC.EUEST9YZY,EstrEuribZC.EUES10YZY]/100;
ESTR = EstrEuribZC.EUROSTR/100;

EURIBZY=[EstrEuribZC.EURONZRZY,EstrEuribZC.EURTNZRZY,EstrEuribZC.EUR1WZRZY,EstrEuribZC.EUR3MZRZY,EstrEuribZC.EUR3MZRZY,EstrEuribZC.EUR6MZRZY,EstrEuribZC.EUR9MZRZY,...
    EstrEuribZC.EUR1YZRZY,EstrEuribZC.EUR2YZRZY,EstrEuribZC.EUR3YZRZY,EstrEuribZC.EUR4YZRZY,EstrEuribZC.EUR5YZRZY,EstrEuribZC.EUR6YZRZY,...
    EstrEuribZC.EUR7YZRZY,EstrEuribZC.EUR8YZRZY,EstrEuribZC.EUR9YZRZY,EstrEuribZC.EUR10YZRZY]/100;

[~, N_mat_Estr] = size(EstrZY);
[~, N_mat_Euribor] = size(EURIBZY);

% Configure maturities (same maturities for each obervation time)
lag=2; %temporal lag
t = EstrEuribZC.Date;
start=t+lag;
matZC=repmat(start,1,N_mat_Estr);
matZC(:,1:2)=start+days(1:2);
matZC(:,3)=start+days(7);
m_shift=[1,3,6,9];
matZC(:,4)=datetime(year(start),month(start)+m_shift(1),day(start));
matZC(:,5)=datetime(year(start),month(start)+m_shift(2),day(start));
matZC(:,6)=datetime(year(start),month(start)+m_shift(3),day(start));
matZC(:,7)=datetime(year(start),month(start)+m_shift(4),day(start));
y_shift=1:10;
for i = 1:10
    matZC(:,8+(i-1)) = datetime(year(start)+y_shift(i),month(start),day(start));
end

div=360;
mat=days(matZC-t)./div;% maturities in years act/360
   

%choose the time windows for estimation

t0 = 768; %no negative rates
%t0 = 2; %all period

%All observed maturities 
mat_Estr = mat(t0:end, 1:N_mat_Estr);
matZC_Estr = matZC(t0:end, 1:N_mat_Estr);
yEster =  log(1+EstrZY(t0:end, 1:N_mat_Estr));
mat_Euribor = mat_Estr;
matZC_Euribor = matZC_Estr;
yEuribor =  log(1+EURIBZY(t0:end,1:N_mat_Euribor));

%maturities on Estr OIS used to calibrate the model
mat_in_Estr_cal = 8;
mat_end_Estr_cal = N_mat_Estr;
mat_Estr_cal = mat(t0:end, mat_in_Estr_cal:mat_end_Estr_cal);
matZC_Estr_cal = matZC(t0:end,mat_in_Estr_cal:mat_end_Estr_cal);
yEster_cal =  log(1+EstrZY(t0:end,mat_in_Estr_cal:mat_end_Estr_cal));

ESTR = ESTR(t0:end);
t_y = EstrEuribZC.Date(t0:end);

[~, N_mat_Estr_cal] = size(yEster_cal);
[Nobs, N_mat_Euribor] = size(yEuribor);

dt = 0; % ahead prediction
T_est = Nobs-dt; % length of estimation window

%Dates of BCE meetings. Jumps may occur at BCE meetings plus 6 days
t_J = BCE_meetings+days(6);
t_J = [t_J; t_J(end)+ [42:42:42*80]']; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimation

% transformed hyperparameter:   
% k_xi = exp(hyper_par(1));
% sig_xi = exp(hyper_par(2));
% k_theta = exp(hyper_par(3));
% theta_bar = hyper_par(4);
% sig_theta = exp(hyper_par(5));
% w_J = exp(hyper_par(6));% w_J % std dev of jumps
% rho = 1 - 2*exp(hyper_par(7))/(1+exp(hyper_par(7)));% correlazione tra le due variabili di stato xi e theta
%
% k_spread = exp(hyper_par(8));
% sig_spread = exp(hyper_par(9));
% spread_bar = hyper_par(10);
%
% sigma = exp(hyper_par(11));  % std dev measurement error 
%
% lambda_xi = exp(hyper_par(12)); % market price of risk of xi
% lambda_theta = exp(hyper_par(13)); % market price of risk of theta
% gamma = exp(hyper_par(14)); % change of measure fo jumps
%
% lambda_spread = exp(hyper_par(15)); % market price of risk of spread

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ProbMeas = 'Change of meausure';
%ProbMeas = 'RN';

% Change of measure
switch ProbMeas
    case 'RN'
        hyper_par0 = [0.4258; -4.7781; -0.1263; 0.0009;  -5.7466; -4.9906;  10.8253;  0.4891; -3.4068; 0.0031;  -7.7348]; 
    case 'Change of meausure'
        hyper_par0 = [0.1325; -5.2107; 0.2559; 0.0005;  -5.7630; -6.1829;  10.8253; 0.4720; -3.4726; 0.0017; -7.4858; -4.6052; -4.6052; -0.0086 ; -4.6052];  
end

% Estimate parameters using maximum likelihood
f  = @(hyper_par)LogLik_fn(hyper_par, [yEster_cal, yEuribor], ESTR, t_y, [matZC_Estr_cal, matZC_Euribor], t_J,'4F',N_mat_Estr_cal,N_mat_Euribor); 
opts = optimset('Display','iter','TolX',1e-9,'TolFun',1e-9,...
                'Diagnostics','on', 'MaxIter',1000, 'MaxFunEvals', 10000,...
                'LargeScale', 'off', 'PlotFcns', @optimplotfval);
[hyper_par, fval, exitflag, output, grad, hessian] = fminunc(f, hyper_par0, opts);
 
% Display estimated results
disp(['MLE of transformed parameters:  ', num2str(hyper_par')]);
disp(['Likelihood:  ', num2str(-1.0*fval)]);
k_xi = exp(hyper_par(1));
sig_xi = exp(hyper_par(2));
k_theta = exp(hyper_par(3));
theta_bar = hyper_par(4);
sig_theta = exp(hyper_par(5));
w_J = exp(hyper_par(6));% w_J % std dev of jumps
rho = 0.9*(1 - 2*exp(hyper_par(7))/(1+exp(hyper_par(7))));% correlazione tra le due variabili di stato xi e theta
k_spread = exp(hyper_par(8));
sig_spread = exp(hyper_par(9));
spread_bar = hyper_par(10);
sigma = exp(hyper_par(11));  % std dev measurement error 

cov_par = inv(hessian);
sig_CI = real(sqrt(diag(cov_par)));
CI_k_xi = k_xi*[1-1.96*sig_CI(1) 1+1.96*sig_CI(1)];
CI_sig_xi = sig_xi*[1-1.96*sig_CI(2) 1+1.96*sig_CI(2)];
CI_k_theta = k_theta*[1-1.96*sig_CI(3) 1+1.96*sig_CI(3)];
CI_theta_bar = [theta_bar-1.96/sqrt(hessian(4,4)) theta_bar+1.96/sqrt(hessian(4,4))];
CI_sig_theta = sig_theta*[1-1.96*sig_CI(5) 1+1.96*sig_CI(5)];
CI_w_J = w_J*[1-1.96*sig_CI(6) 1+1.96*sig_CI(6)];
CI_rho = rho*[1-1.96*sig_CI(7) 1+1.96*sig_CI(7)];
CI_k_spread = k_spread*[1-1.96*sig_CI(8) 1+1.96*sig_CI(8)];
CI_spread_bar = [spread_bar-1.96/sqrt(hessian(9,9)) spread_bar+1.96/sqrt(hessian(9,9))];
CI_sig_spread = sig_spread*[1-1.96*sig_CI(10) 1+1.96*sig_CI(10)];
CI_sigma = sigma*[1-1.96*sig_CI(11) 1+1.96*sig_CI(11)];


if length(hyper_par) == 11
    lambda_xi = 0; % market price of risk of xi
    lambda_theta = 0; % market price of risk of theta
    gamma_J = 0; % market price of risk for jumps
    lambda_spread = 0;
else
    lambda_xi = exp(hyper_par(12)); % market price of risk of xi
    lambda_theta = exp(hyper_par(13)); % market price of risk of theta
    gamma_J = exp(hyper_par(14)); % market price of risk for jumps
    lambda_spread = exp(hyper_par(15)); % market price of risk for spread
    CI_lambda_xi = lambda_xi*[1-1.96*sig_CI(12) 1+1.96*sig_CI(12)];
    CI_lambda_theta = lambda_theta*[1-1.96*sig_CI(13) 1+1.96*sig_CI(13)];
    CI_gamma_J = gamma_J*[1-1.96*sig_CI(14) 1+1.96*sig_CI(14)];
    CI_lambda_spread = lambda_spread*[1-1.96*sig_CI(15) 1+1.96*sig_CI(15)];
end

disp('ESTIMATION RESULTS         ');
disp(['k_xi = ', num2str(k_xi), '     CI: ', num2str(CI_k_xi)]);
disp(['sigma_xi = ', num2str(sig_xi), '     CI: ', num2str(CI_sig_xi)]);
disp(['k_theta = ', num2str(k_theta), '     CI: ', num2str(CI_k_theta)]);
disp(['theta_bar = ', num2str(theta_bar), '     CI: ', num2str(CI_theta_bar)]);
disp(['sigma_theta = ', num2str(sig_theta), '     CI: ', num2str(CI_sig_theta)]);
disp(['w_J = ', num2str(w_J), '     CI: ', num2str(CI_w_J)]);
disp(['rho = ', num2str(rho), '     CI: ', num2str(CI_rho)]);
disp(['k_spread = ', num2str(k_spread), '     CI: ', num2str(CI_k_spread)]);
disp(['spread_bar = ', num2str(spread_bar), '     CI: ', num2str(CI_spread_bar)]);
disp(['sigma_spread = ', num2str(sig_spread), '     CI: ', num2str(CI_sig_spread)]);
disp(['sigma = ', num2str(sigma), '     CI: ', num2str(CI_sigma)]);

if length(hyper_par) == 15
    disp(['lambda_xi ', num2str(lambda_xi), '     CI: ', num2str(CI_lambda_xi)]);
    disp(['lambda_theta ', num2str(lambda_theta), '     CI: ', num2str(CI_lambda_theta)]);
    disp(['gamma_J', num2str(gamma_J), '     CI: ', num2str(CI_gamma_J)]);
    disp(['lambda_xspread ', num2str(lambda_spread), '     CI: ', num2str(CI_lambda_spread)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create and estimate the state-space model
[Z, d, GG, T, c, HH, a, P] = SetStateSpaceModel_4F(hyper_par, ESTR, t_y, [matZC_Estr matZC_Euribor], t_J, N_mat_Estr, N_mat_Euribor);
[Inn, F, StatePred, CovStatePred, LogLik] = KalmanFilter([yEster yEuribor], Z, d, GG, T, c, HH, a, P);

%% Smoothed estimates of latent states
[State, VarState] = KFStateSmoother([yEster yEuribor], Z, d, GG, T, c, HH, a, P);

x_est1 = State(1,:); %xi
Varx_est1 = VarState(1,:);
x_est2 = State(2,:); % theta
Varx_est2 = VarState(2,:);
x_est3 = State(3,:); % r
Varx_est3 = VarState(3,:);
x_est4 = State(4,:); % spread
Varx_est4 = VarState(4,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model term structure

yEster_model = yEster;
yEuribor_model = yEuribor;


for i = 1:Nobs
    yEster_model(i,1:N_mat_Estr) = d(1:N_mat_Estr,i) + Z(1:N_mat_Estr,:,i)*State(:,i); 
    yEuribor_model(i,1:N_mat_Euribor) = d(N_mat_Estr+1:end,i) + Z(N_mat_Estr+1:end,:,i)*State(:,i);         
end

SE_Ester=(yEster_model-yEster).^2;
RSE_Ester=sqrt(mean(SE_Ester,2));
RMSE_Ester=sqrt(mean(SE_Ester))*10000;

SE_Euribor=(yEuribor_model-yEuribor).^2;
RSE_Euribor=sqrt(mean(SE_Euribor,2));
RMSE_Euribor=sqrt(mean(SE_Euribor))*10000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
d_y = datestr(now, 'dd-mm-yyyy HH-MM-SS');
filename = ['4F ' d_y '.mat'];
save(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure
h1 = figure(1);
subplot(2,1,1)
plot(t_y, x_est1, 'b', ...
     t_y, x_est1+1.96 * sqrt(Varx_est1), '-g',...
     t_y, x_est1-1.96 * sqrt(Varx_est1), '-g'); 
title('Estimate of the latent variable \xi');
subplot(2,1,2)
plot(t_y, x_est2, 'b', ...
     t_y, x_est2+1.96 * sqrt(Varx_est2), '-g',...
     t_y, x_est2-1.96 * sqrt(Varx_est2), '-g'); 
title('Estimate of the latent variable \theta');
%SaveFigureFullScreenPDF(h1,'Estimate of the latent variable')

h2 = figure(2);
plot(t_y, ESTR)
hold on;
plot(t_y, x_est3)
plot(t_y, yEster(:,1))
hold off;
% xline(BCE_meetings,':')
legend('ESTR','filtered short rate','y(1) yield','Location','SouthEast')
% xlim([t_y(1),t_y(end)])


%t_val = 1;
% t_val = 120;
% t_val = 100;
t_val = length(t_y);
h3 = figure(3);
subplot(2,1,1)
plot(mat_Estr(t_val,:), 100*yEster(t_val,:),'o');
hold on
plot(mat_Estr(t_val,:), 100*yEster_model(t_val,:))
hold off
legend('market','model')
xlabel('maturities (years)')
title(['Ester yield fit at ', datestr(t_y(t_val))])
subplot(2,1,2)
plot(mat_Euribor(t_val,:), 100*yEuribor(t_val,:),'o');
hold on
plot(mat_Euribor(t_val,:), 100*yEuribor_model(t_val,:))
hold off
legend('market','model')
xlabel('maturities (years)')
title(['Euribor yield fit at ', datestr(t_y(t_val))])
%SaveFigureFullScreenPDF(h3,['fit at ', datestr(t_y(t_val))])


h4 = figure(4);
sgtitle('Estimate for Estr')
subplot(3,3,1)
i = 1;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
plot(t_y,ESTR,'.g');
%legend('Model','observed','ESTR','Location','NorthWest');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i)*div,'%.0f'),' days'])
subplot(3,3,2)
i = 3;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i)*div,'%.0f'),' days'])
subplot(3,3,3)
i = 5;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.2f'),' years'])
subplot(3,3,4)
i = 6;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.1f'),' years'])
subplot(3,3,5)
i = 8;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' year'])
subplot(3,3,6)
i = 9;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' years'])
subplot(3,3,7)
i = 12;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' years'])
subplot(3,3,8)
i = 14;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' years'])
subplot(3,3,9)
i = 17;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' years'])
%SaveFigureFullScreenPDF(h4,'Estr Yields to maturity')


h5 = figure(5);
sgtitle('Estimate for Euribor')
subplot(3,3,1)
i = 1;
plot(t_y,yEuribor_model(:,i))
hold on;
plot(t_y,yEuribor(:,i),'.r');
plot(t_y,ESTR,'.g');
%legend('Model','observed','ESTR','Location','NorthWest');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i)*div,'%.0f'),' days'])
subplot(3,3,2)
i = 3;
plot(t_y,yEuribor_model(:,i))
hold on;
plot(t_y,yEuribor(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i)*div,'%.0f'),' days'])
subplot(3,3,3)
i = 5;
plot(t_y,yEuribor_model(:,i))
hold on;
plot(t_y,yEuribor(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i),'%.2f'),' years'])
subplot(3,3,4)
i = 6;
plot(t_y,yEuribor_model(:,i))
hold on;
plot(t_y,yEuribor(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i),'%.1f'),' years'])
subplot(3,3,5)
i = 8;
plot(t_y,yEuribor_model(:,i))
hold on;
plot(t_y,yEuribor(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i),'%.0f'),' year'])
subplot(3,3,6)
i = 9;
plot(t_y,yEuribor_model(:,i))
hold on;
plot(t_y,yEuribor(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i),'%.0f'),' years'])
subplot(3,3,7)
i = 12;
plot(t_y,yEuribor_model(:,i))
hold on;
plot(t_y,yEuribor(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i),'%.0f'),' years'])
subplot(3,3,8)
i = 14;
plot(t_y,yEuribor_model(:,i))
hold on;
plot(t_y,yEuribor(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i),'%.0f'),' years'])
subplot(3,3,9)
i = 17;
plot(t_y,yEuribor_model(:,i))
hold on;
plot(t_y,yEuribor(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Euribor(1,i),'%.0f'),' years'])
%SaveFigureFullScreenPDF(h5,'Euribor Yields to maturity')



figure(6);
plot(t_y,RSE_Euribor*10^4,'.')
title('Root mean square error for Euribor')
ylabel('bsp')


figure(7);
RSE_Ester_2=sqrt(mean(SE_Ester(:,mat_in_Estr_cal:mat_end_Estr_cal),2));
plot(t_y,RSE_Ester_2*10^4)
title('Root mean square error for Estr')
ylabel('bsp')

figure(8);
SE_total=[SE_Ester(:,mat_in_Estr_cal:mat_end_Estr_cal),SE_Euribor];
RSE_total=sqrt(mean(SE_total,2));
plot(t_y,RSE_total*10^4)
title('Root mean square error')
ylabel('bsp')


figure(9)
sgtitle('Latent variables estimates')
subplot(2,2,1)
plot(t_y, x_est1, '.') 
ylim([-8,15]*10^-3)
title('\xi');
subplot(2,2,2)
plot(t_y, x_est2, '.')
ylim([-8,15]*10^-3)
title('\theta');
subplot(2,2,3)
plot(t_y, x_est3, '.')
title('r');
subplot(2,2,4)
plot(t_y, x_est4, '.')
ylim([-8,15]*10^-3)
title('s');

%% plot of the expected value E_tR_T

t=t_y(end);
ind=find(t_J>t);
N_forecast=5;
ind=ind(1:N_forecast);
next_jumps=t_J(ind);
forecasts=zeros(size(next_jumps));
T2=t+year(10);
par1={k_xi,sig_xi};
par2={k_theta,theta_bar,sig_theta};
aJ=0;
gJ=1;
wJ=0;
parJ={aJ,gJ,wJ};
r_t=x_est3(end);
X_t=[x_est1(end);x_est2(end)];

for i=1:N_forecast
    [Ey,forecasts(i)]=Expected_yield(t,next_jumps(i),T2, t_J, par1,par2,parJ,rho,lambda_xi,lambda_theta,r_t,X_t);
end

figure(10)
plot(next_jumps,100*forecasts,':o')
xticks(next_jumps);
datetick('x', 'dd-mmm', 'keepticks');
title(['Forecast at time ' datestr(t) ' of the next jumps'])
xlabel('Times of jumps')
ylabel('%')