clear all
close all
clc


S0 = 100; 
K = [90 100 110]; 
sigma = 0.3; 
r = 0.04;
T = 1; 
dt = [T/4 T/12 T/50]; 
M = 1E7; 
alpha = 0.95;


% Preallocate matrices for direct arithmetic Asian option simulation
PriceArithDirect = zeros(3,3);
MCStdArithDirect = zeros(3,3);
ConfIntDirectLower = zeros(3,3);
ConfIntDirectUpper = zeros(3,3);

% Direct arithmetic Asian option price simulation loop
for i = 1:3
    for j = 1:3
        tic;
        [MCstd, MCAsianArithPrice, MCConfInt] = AsianArithMC(S0, K(i), sigma, r, T, dt(j), M, alpha);
        simTime = toc

        PriceArithDirect(i,j) = MCAsianArithPrice;
        MCStdArithDirect(i,j) = MCstd;
        ConfIntDirectLower(i,j) = MCConfInt(1);
        ConfIntDirectUpper(i,j) = MCConfInt(2);
    end
end

% Creating K and dt vectors
[K_mesh, dt_mesh] = meshgrid(K, dt);
K_vals = K_mesh(:);
dt_vals = dt_mesh(:);

% Generate results table
DirectTable = table(K_vals, dt_vals, ...
    reshape(PriceArithDirect', [], 1), reshape(MCStdArithDirect', [], 1), ...
    reshape(ConfIntDirectLower', [], 1), reshape(ConfIntDirectUpper', [], 1), ...
    'VariableNames', {'K', 'dt', 'Price_Direct', 'MCStd_Direct', 'ConfInt_Lower', 'ConfInt_Upper'});

disp('Direct Arithmetic Asian Option Price Table:')
disp(DirectTable)


% ----- Function for Direct Arithmetic Asian Option Pricing -----
function [MCstd, MCAsianArithPrice, MCConfInt] = AsianArithMC(S0, K, sigma, r, T, dt, M, alpha)

n = T/dt;
C = zeros(1,M);

for j = 1:M
    S = zeros(1,n+1);
    S(1) = S0;
    for i = 1:n
        S(i+1) = S(i)*exp((r - sigma^2/2)*dt + sigma*sqrt(dt)*randn);
    end
    C(j) = exp(-r*T)*max(mean(S(2:end))-K,0);
end

MCAsianArithPrice = mean(C);
MCstd = std(C)/sqrt(M);
MCConfInt = MCAsianArithPrice + norminv(0.5 + alpha/2)*MCstd*[-1 1];

end
