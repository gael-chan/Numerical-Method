clear all
close all
clc


S0 = 100; 
K = [90 100 110]; 
sigma = 0.3; 
r = 0.04;
T = 1; 
dt = [T/4 T/12 T/50]; 
M = 1E5; 
alpha = 0.95;
Mb = 1E4;

% Preallocating matrices for LB as control variate
PriceArithLBasCV = zeros(3,3);
MCStdArithLBasCV = zeros(3,3);
ConfInt_LB_Lower = zeros(3,3);
ConfInt_LB_Upper = zeros(3,3);
EfficientRatioA = zeros(3,3);

% Preallocating matrices for Gn as control variate
PriceArithGMasCV = zeros(3,3);
MCStdArithGMasCV = zeros(3,3);
ConfInt_GM_Lower = zeros(3,3);
ConfInt_GM_Upper = zeros(3,3);
EfficientRatioB = zeros(3,3);

% Loop for LB as control variate
for i = 1:3
    for j = 1:3
        tic;
        [correlcoefLB, bstarhat] = AsianMCCtrlVarEstimateb(S0, K(i), sigma, r, T, dt(j), Mb);
        [MCstd, MCAsianArithPrice, MCConfInt] = AsianArithMCCtrlVar(S0, K(i), sigma, r, T, dt(j), M, bstarhat, alpha);
        TimeA = toc;

        PriceArithLBasCV(i,j) = MCAsianArithPrice;
        MCStdArithLBasCV(i,j) = MCstd;
        ConfInt_LB_Lower(i,j) = MCConfInt(1);
        ConfInt_LB_Upper(i,j) = MCConfInt(2);
        SquaredSEArithLBasCV = MCstd^2/M;
        EfficientRatioA(i,j) = TimeA * SquaredSEArithLBasCV;
        CorrelationCoeffLB(i,j) = correlcoefLB; % Store correlation coefficient for LB
    end
end

% Loop for Gn as control variate
for i = 1:3
    for j = 1:3
        tic;
        [correlcoefGM, bstarhat] = AsianMCCtrlVarEstimateb2(S0, K(i), sigma, r, T, dt(j), Mb);
        [MCstd, MCAsianArithPrice, MCConfInt] = AsianArithMCCtrlVar2(S0, K(i), sigma, r, T, dt(j), M, bstarhat, alpha);
        TimeB = toc;

        PriceArithGMasCV(i,j) = MCAsianArithPrice;
        MCStdArithGMasCV(i,j) = MCstd;
        ConfInt_GM_Lower(i,j) = MCConfInt(1);
        ConfInt_GM_Upper(i,j) = MCConfInt(2);
        SquaredSEArithGMasCV = MCstd^2/M;
        EfficientRatioB(i,j) = TimeB * SquaredSEArithGMasCV;
        CorrelationCoeffGM(i,j) = correlcoefGM; % Store correlation coefficient for GM
    end
end



% Create vectors for K and dt pairing
[K_mesh, dt_mesh] = meshgrid(K, dt); 
K_vals = K_mesh(:);  
dt_vals = dt_mesh(:);  

% Correct table for LB as control variate
LBTable = table(K_vals, dt_vals, ...
    reshape(PriceArithLBasCV', [], 1), reshape(MCStdArithLBasCV', [], 1), ...
    reshape(ConfInt_LB_Lower', [], 1), reshape(ConfInt_LB_Upper', [], 1), ...
    reshape(EfficientRatioA', [], 1), ...
    'VariableNames', {'K', 'dt', ...
                      'Price_LB', 'MCStd_LB', 'ConfInt_LB_Lower', 'ConfInt_LB_Upper', ...
                      'EfficientRatioA'});

% Correct table for GM as control variate
GNTable = table(K_vals, dt_vals, ...
    reshape(PriceArithGMasCV', [], 1), reshape(MCStdArithGMasCV', [], 1), ...
    reshape(ConfInt_GM_Lower', [], 1), reshape(ConfInt_GM_Upper', [], 1), ...
    reshape(EfficientRatioB', [], 1), ...
    'VariableNames', {'K', 'dt', ...
                      'Price_GM', 'MCStd_GM', 'ConfInt_GM_Lower', 'ConfInt_GM_Upper', ...
                      'EfficientRatioB'});

% Efficiency Ratio Comparison Table
EfficiencyComparison = EfficientRatioA ./ EfficientRatioB;
EfficiencyTable = table(K_vals, dt_vals, reshape(EfficiencyComparison, [], 1), ...
    'VariableNames', {'K', 'dt', 'EfficiencyRatioComparison'});

disp('Efficiency Ratio Comparison:')
disp(EfficiencyTable)

disp('Table for LB as control variate:')
disp(LBTable)

disp('Table for Gn as control variate:')
disp(GNTable)



%----------------------------------------------------------------------------------------------------------------------------------

function [MCstd, MCAsianArithPrice, MCConfInt] = AsianArithMCCtrlVar(S0, K, sigma, r, T, dt, M, bstarhat, alpha)

n = T/dt;
C = zeros(1,M); 
P = zeros(1,M); 
for j = 1:M
    S = [S0 zeros(1,n)];
    for i = 1:n
        S(i+1) = S(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*randn); 
    end
    C(j) = exp(-r*T)*max(mean(S(2:n+1))-K,0); 
    P(j) = exp(-r*T)*(mean(S(2:n+1))-K)*((prod(S(2:n+1))^(1/n))>K); 
end
PTrue = ArithLB(S0, K, r, T, dt, sigma);
Cb = C - bstarhat*(P-PTrue); 
MCAsianArithPrice = mean(Cb); 
MCstd = std(Cb)/sqrt(M); 
format short g
MCConfInt = MCAsianArithPrice + norminv(0.5+alpha/2)*MCstd*[-1 1]; 
end



function [correlcoef, bstarhat] = AsianMCCtrlVarEstimateb(S0, K, sigma, r, T, dt, Mb)

n = T/dt; 
C = zeros(1,Mb); 
P = zeros(1,Mb); 
for j = 1:Mb
    S = [S0 zeros(1,n)];
    for i = 1:n
        S(i+1) = S(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*randn); 
    end
    C(j) = exp(-r*T)*max(mean(S(2:n+1))-K,0); 
    P(j) = exp(-r*T)*(mean(S(2:n+1))-K)*((prod(S(2:n+1))^(1/n))>K); 
end
regress_coefs = regress(C', [ones(Mb,1) P']); 
bstarhat = sum((C-mean(C)).*(P-mean(P)))/sum((P-mean(P)).^2); 
correlcoef = sum((C-mean(C)).*(P-mean(P)))/sqrt(sum((C-mean(C)).^2)*sum((P-mean(P)).^2));

end

function LBtrue = ArithLB(S0, K, r, T, dt, sigma)
n = T/dt;
sum = 0;

for k = 1:n
    mu_k = (r-sigma^2/2)*k*dt;
    sigma_k = sigma*sqrt(k*dt);
    a_k = sigma*sqrt(dt)* (k*(n+1-(k+1)/2))/(sqrt(n*(n+1)*(2*n+1)/6));
    sigmaBar = sigma*sqrt((2*n+1)/(3*n));
    Tbar = (n+1)*dt/2;
    b = (log(S0/K)+(r-sigma^2/2)*Tbar)/(sigmaBar*sqrt(Tbar));
    sum = sum + (exp(mu_k+sigma_k^2/2)*normcdf(b+a_k));
end
LBtrue = S0*exp(-r*T)/n * sum -K*exp(-r*T)*normcdf(b);

end


%----------------------------------------------------------------------------------------------------------------------------------
function [MCstd, MCAsianArithPrice, MCConfInt] = AsianArithMCCtrlVar2(S0, K, sigma, r, T, dt, M, bstarhat, alpha)

n = T/dt; 
C = zeros(1,M); 
P = zeros(1,M); 
for j = 1:M
    S = [S0 zeros(1,n)];
    for i = 1:n
        S(i+1) = S(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*randn); 
    end
    C(j) = exp(-r*T)*max(mean(S(2:n+1))-K,0);
    P(j) = exp(-r*T)*(max((prod(S(2:n+1))^(1/n)-K),0)); 
end
PTrue = GMAsianTrue(S0, K, r, T, dt, sigma);
Cb = C - bstarhat*(P-PTrue); 
MCAsianArithPrice = mean(Cb); 
MCstd = std(Cb)/sqrt(M); 
format short g
MCConfInt = MCAsianArithPrice + norminv(0.5+alpha/2)*MCstd*[-1 1]; 
end



function [correlcoef, bstarhat] = AsianMCCtrlVarEstimateb2(S0, K, sigma, r, T, dt, Mb)

n = T/dt; 
C = zeros(1,Mb);  
P = zeros(1,Mb); 
for j = 1:Mb
    S = [S0 zeros(1,n)];
    for i = 1:n
        S(i+1) = S(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*randn); 
    end
    C(j) = exp(-r*T)*max(mean(S(2:n+1))-K,0);
    P(j) = exp(-r*T)*(max((prod(S(2:n+1))^(1/n)-K),0)); 
end

regress_coefs = regress(C', [ones(Mb,1) P']); 
bstarhat = sum((C-mean(C)).*(P-mean(P)))/sum((P-mean(P)).^2); 
correlcoef = sum((C-mean(C)).*(P-mean(P)))/sqrt(sum((C-mean(C)).^2)*sum((P-mean(P)).^2));
end


function GnTrue = GMAsianTrue(S0, K, r, T, dt, sigma)
 n = T/dt;
 sigmaBar = sigma*sqrt((2*n+1)/(3*n));
 Tbar = (n+1)*dt/2;
 d = (log(S0/K)+(r-sigma^2/2+sigmaBar^2)*Tbar)/(sigmaBar * sqrt(Tbar));
 GnTrue = S0*exp((r-(sigma^2)/2+(sigmaBar^2)/2)*Tbar-r*T)*normcdf(d) - K*exp(-r*T)*normcdf(d-sigmaBar*sqrt(Tbar));

end
%----------------------------------------------------------------------------------------------------------------------------------
