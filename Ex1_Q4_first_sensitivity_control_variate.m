%Question 4 Part 2 ------------------------------------------------------------------------------------------------------------------------------
clear all
close all
clc

S0 = 100; K = [90 100 110]; sigma = 0.3; r = 0.04;T=1; dt=[T/4 T/12 T/50]; 
M = 1E5; alpha = 0.95;
Mb = 1E4;

DeltaArithLBasCV = [];
MCStdArithLBasCV =[];
MCConfintArithLBasCV = [];
EfficientRatioA = [];
col = 0;
for i = 1:3
    for j = 1:3
    col = (j - 1) * 2;
    tic;
    [correlcoef, bstarhat] = AsianMCCtrlVarEstimateb(S0, K(i), sigma, r, T, dt(j), Mb);
    [MCstd, MCAsianArithDelta, MCConfInt] = AsianArithMCCtrlVar(S0, K(i), sigma, r, T, dt(j), M, bstarhat, alpha);
    TimeA = toc;
    DeltaArithLBasCV(i,j) = MCAsianArithDelta;
    MCStdArithLBasCV(i,j) = MCstd;
    MCConfintArithLBasCV(i,col+1:col+2)=MCConfInt;
    SquaredSEArithLBasCV = MCstd^2/M;
    EfficientRatioA(i,j) = TimeA * SquaredSEArithLBasCV;
    end
end
display('LB_LR as the control variate')


DeltaArithGnAsCV = [];
MCStdArithGnAsCV =[];
MCConfintArithGnAsCV = [];
EfficientRatioB = [];
col = 0;
for i = 1:3
    for j = 1:3
    col = (j - 1) * 2;
    tic;
    [correlcoef, bstarhat] = AsianMCCtrlVarEstimateb1(S0, K(i), sigma, r, T, dt(j), Mb);
    [MCstd, MCAsianArithDelta, MCConfInt] = AsianArithMCCtrlVar1(S0, K(i), sigma, r, T, dt(j), M, bstarhat, alpha);
    TimeB = toc;
    DeltaArithGnAsCV(i,j) = MCAsianArithDelta;
    MCStdArithGnAsCV(i,j) = MCstd;
    MCConfintArithGnAsCV(i,col+1:col+2)=MCConfInt;
    SquaredSEArithGnAsCV = MCstd^2/M;
    EfficientRatioB(i,j) = TimeB * SquaredSEArithGnAsCV;

    end
end
% Create vectors for K and dt pairing
[K_mesh, dt_mesh] = meshgrid(K, dt);
K_vals = K_mesh(:);  % Flattened vector of K values
dt_vals = dt_mesh(:);  % Flattened vector of dt values

% Correct table for LB as control variate
LBTable = table(K_vals, dt_vals, ...
    reshape(DeltaArithLBasCV', [], 1), reshape(MCStdArithLBasCV', [], 1), ...
    reshape(MCConfintArithLBasCV(:,1:2:end)', [], 1), reshape(MCConfintArithLBasCV(:,2:2:end)', [], 1), ...
    reshape(EfficientRatioA', [], 1), ...
    'VariableNames', {'K', 'dt', ...
                      'Delta_LB', 'MCStd_LB', 'ConfInt_LB_Lower', 'ConfInt_LB_Upper', ...
                      'EfficientRatioA'});

% Correct table for GM as control variate
GNTable = table(K_vals, dt_vals, ...
    reshape(DeltaArithGnAsCV', [], 1), reshape(MCStdArithGnAsCV', [], 1), ...
    reshape(MCConfintArithGnAsCV(:,1:2:end)', [], 1), reshape(MCConfintArithGnAsCV(:,2:2:end)', [], 1), ...
    reshape(EfficientRatioB', [], 1), ...
    'VariableNames', {'K', 'dt', ...
                      'Delta_GM', 'MCStd_GM', 'ConfInt_GM_Lower', 'ConfInt_GM_Upper', ...
                      'EfficientRatioB'});

% Displaying the corrected tables
disp('Table for LB as control variate:')
disp(LBTable)

disp('Table for Gn as control variate:')
disp(GNTable)

% Efficiency Ratio Comparison Table
EfficiencyComparison = EfficientRatioA ./ EfficientRatioB;
EfficiencyTable = table(K_vals, dt_vals, reshape(EfficiencyComparison', [], 1), ...
    'VariableNames', {'K', 'dt', 'EfficiencyRatioComparison'});

disp('Efficiency Ratio Comparison:')
disp(EfficiencyTable)


%------------------------------------------------------------------------------------------------------------------------------
function [MCstd, MCAsianArithDelta, MCConfInt] = AsianArithMCCtrlVar(S0, K, sigma, r, T, dt, M, bstarhat, alpha)
b = r-sigma^2/2;
n = T/dt; % no. of monitoring dates of the underlying
LRDeltaAn = zeros(1,M);
LRDeltaLB = zeros(1,M);
for j = 1:M
    S = [S0 zeros(1,n)];
    for i = 1:n
        S(i+1) = S(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*randn); % simulation of stock price trajectory (exact method)
    end
    zetaS1 = (log(S(2)/S(1))-b*dt)/(sigma*sqrt(dt)); % computes the zeta function at S1 (i.e., S(2)) conditional on S0 (i.e., S(1))
    LRDeltaAn(j) = exp(-r*T)*max(mean(S(2:n+1))-K,0)*(zetaS1/(S0*sigma*sqrt(dt))); % simulates samples of the LR delta of An
    LRDeltaLB(j) = exp(-r*T)*(mean(S(2:n+1))-K)*((prod(S(2:n+1))^(1/n))>K)*(zetaS1/(S0*sigma*sqrt(dt))); % simulates samples of the LR delta of LB
end

LRDeltaLBTrue = LBfirstDiff(S0, K, r, T, dt, sigma); % calls the true pricing formula of the geometric Asian option
LRDeltaAnb = LRDeltaAn - bstarhat*(LRDeltaLB-LRDeltaLBTrue); % computes the sample of CV discounted payoffs of the arithmetic Asian option
MCAsianArithDelta = mean(LRDeltaAnb); % computes the CV price estimate of the arithmetic Asian option
MCstd = std(LRDeltaAnb)/sqrt(M); % computes the standard error of the CV price estimate of the arithmetic Asian option
format long g
MCConfInt = MCAsianArithDelta + norminv(0.5+alpha/2)*MCstd*[-1 1]; % 100*alpha% confidence interval
end

% The following estimates the optimal coefficient b* in the control variate
% method using a separate set of Monte Carlo simulations. It returns the
% correlation coefficient between samples of discount payoffs of the arithmetic
% and geometric Asian options, and the estimate b*.

function [correlcoef, bstarhat] = AsianMCCtrlVarEstimateb(S0, K, sigma, r, T, dt, Mb)

n = T/dt; % no. of monitoring dates of the underlying
b = r-sigma^2/2;
LRDeltaAn = zeros(1,Mb);
LRDeltaLB = zeros(1,Mb);
for j = 1:Mb
    S = [S0 zeros(1,n)];
    for i = 1:n
        S(i+1) = S(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*randn); % simulation of stock price trajectory (exact method)
    end
    zetaS1 = (log(S(2)/S(1))-b*dt)/(sigma*sqrt(dt)); % computes the zeta function at S1 (i.e., S(2)) conditional on S0 (i.e., S(1))
    LRDeltaAn(j) = exp(-r*T)*max(mean(S(2:n+1))-K,0)*(zetaS1/(S0*sigma*sqrt(dt))); % simulates samples of the LR delta of An
    LRDeltaLB(j) = exp(-r*T)*(mean(S(2:n+1))-K)*((prod(S(2:n+1))^(1/n))>K)*(zetaS1/(S0*sigma*sqrt(dt))); % simulates samples of the LR delta of LB
end


regress_coefs = regress(LRDeltaAn', [ones(Mb,1) LRDeltaLB']); % regresses C samples on P samples
bstarhat = sum((LRDeltaAn-mean(LRDeltaAn)).*(LRDeltaLB-mean(LRDeltaLB)))/sum((LRDeltaLB-mean(LRDeltaLB)).^2); % computes the estimate of the optimal coefficient b*
correlcoef = sum((LRDeltaAn-mean(LRDeltaAn)).*(LRDeltaLB-mean(LRDeltaLB)))/sqrt(sum((LRDeltaAn-mean(LRDeltaAn)).^2)*sum((LRDeltaLB-mean(LRDeltaLB)).^2)); % computes the correlation coef
% between samples of discounted payoffs of the arithmetic and geometric Asian options
end

function LBdeltaTrue = LBfirstDiff(S0, K, r, T, dt, sigma)
n = T/dt; %Delta

sigBar = sigma*sqrt((2*n+1)/(3*n)); %Sigma bar
tBar = (n+1)*dt/2; %T bar

b = (log(S0/K)+(r-sigma^2/2)*tBar)/(sigBar*sqrt(tBar)); %b

total1 = 0;
total2 = 0;
for m = 1:1:n %Loop to calculate the sums within the final equation
    muK = (r - sigma^2/2)*m*dt;
    sigK = sigma*sqrt(m*dt);
    aK = sigma*dt^0.5*(m*(n+1-(m+1)/2))/(sqrt(n*(n+1)*(2*n+1)/6));

    total1 = total1 + exp(muK+sigK^2/2)*normcdf(b+aK);
    total2 = total2 + (1/(sqrt(2*pi)))*exp(-(b+aK)^2/2)*(1/(S0*sigBar*sqrt(tBar)))*exp(muK+sigK^2/2);
end

LBdeltaTrue = (exp(-r*T)/n)*total1 + ((S0*exp(-r*T))/n)*total2 - K*exp(-r*T)*(1/sqrt(2*pi))*exp(-b^2/2)*(1/(S0*sigBar*sqrt(tBar)));

end

%------------------------------------------------------------------------------------------------------------------------------
function [MCstd, MCAsianArithDelta, MCConfInt] = AsianArithMCCtrlVar1(S0, K, sigma, r, T, dt, M, bstarhat, alpha)
b = r-sigma^2/2;
n = T/dt; % no. of monitoring dates of the underlying
LRDeltaAn = zeros(1,M);
LRDeltaGn = zeros(1,M);
for j = 1:M
    S = [S0 zeros(1,n)];
    for i = 1:n
        S(i+1) = S(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*randn); % simulation of stock price trajectory (exact method)
    end
    zetaS1 = (log(S(2)/S(1))-b*dt)/(sigma*sqrt(dt)); % computes the zeta function at S1 (i.e., S(2)) conditional on S0 (i.e., S(1))
    LRDeltaAn(j) = exp(-r*T)*max(mean(S(2:n+1))-K,0)*(zetaS1/(S0*sigma*sqrt(dt))); % simulates samples of the LR delta of An
    LRDeltaGn(j) = exp(-r*T)*(max((prod(S(2:n+1))^(1/n)-K),0))*(zetaS1/(S0*sigma*sqrt(dt))); % simulates samples of the LR delta of LB
end

LRDeltaGnTrue = GnfirstDiff(S0, K, r, T, dt, sigma); % calls the true pricing formula of the geometric Asian option
LRDeltaAnb = LRDeltaAn - bstarhat*(LRDeltaGn-LRDeltaGnTrue); % computes the sample of CV discounted payoffs of the arithmetic Asian option
MCAsianArithDelta = mean(LRDeltaAnb); % computes the CV price estimate of the arithmetic Asian option
MCstd = std(LRDeltaAnb)/sqrt(M); % computes the standard error of the CV price estimate of the arithmetic Asian option
format long g
MCConfInt = MCAsianArithDelta + norminv(0.5+alpha/2)*MCstd*[-1 1]; % 100*alpha% confidence interval
end

% The following estimates the optimal coefficient b* in the control variate
% method using a separate set of Monte Carlo simulations. It returns the
% correlation coefficient between samples of discount payoffs of the arithmetic
% and geometric Asian options, and the estimate b*.

function [correlcoef, bstarhat] = AsianMCCtrlVarEstimateb1(S0, K, sigma, r, T, dt, Mb)

n = T/dt; % no. of monitoring dates of the underlying
b = r-sigma^2/2;
LRDeltaAn = zeros(1,Mb);
LRDeltaGn = zeros(1,Mb);
for j = 1:Mb
    S = [S0 zeros(1,n)];
    for i = 1:n
        S(i+1) = S(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*randn); % simulation of stock price trajectory (exact method)
    end
    zetaS1 = (log(S(2)/S(1))-b*dt)/(sigma*sqrt(dt)); % computes the zeta function at S1 (i.e., S(2)) conditional on S0 (i.e., S(1))
    LRDeltaAn(j) = exp(-r*T)*max(mean(S(2:n+1))-K,0)*(zetaS1/(S0*sigma*sqrt(dt))); % simulates samples of the LR delta of An
    LRDeltaGn(j) = exp(-r*T)*(max((prod(S(2:n+1))^(1/n)-K),0))*(zetaS1/(S0*sigma*sqrt(dt))); % simulates samples of the LR delta of LB
end


regress_coefs = regress(LRDeltaAn', [ones(Mb,1) LRDeltaGn']); 
bstarhat = sum((LRDeltaAn-mean(LRDeltaAn)).*(LRDeltaGn-mean(LRDeltaGn)))/sum((LRDeltaGn-mean(LRDeltaGn)).^2); % computes the estimate of the optimal coefficient b*
correlcoef = sum((LRDeltaAn-mean(LRDeltaAn)).*(LRDeltaGn-mean(LRDeltaGn)))/sqrt(sum((LRDeltaAn-mean(LRDeltaAn)).^2)*sum((LRDeltaGn-mean(LRDeltaGn)).^2)); % computes the correlation coef
% between samples of discounted payoffs of the arithmetic and geometric Asian options
end

function GndeltaTrue = GnfirstDiff(S0, K, r, T, dt, sigma)
n = T/dt;

sigBar = sigma*sqrt((2*n+1)/(3*n));
tBar = (n+1)*dt/2;

d = (log(S0/K)+(r-sigma^2/2+sigBar^2)*tBar)/(sigBar*sqrt(tBar));

GndeltaTrue= exp(-r*T)* ( exp((r-sigma^2/2+sigBar^2/2)*tBar)*normcdf(d) + S0* exp((r-sigma^2/2+sigBar^2/2)*tBar)*normpdf(d)/(S0*sigBar*sqrt(tBar)) - K*normpdf(d-sigBar*sqrt(tBar))/(S0*sigBar*sqrt(tBar)) );

end



