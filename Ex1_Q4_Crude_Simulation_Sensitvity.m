% Question 4 Part 2 - Direct simulation without control variate
clear all
close all
clc

S0 = 100; K = [90]; sigma = 0.3; r = 0.04; T = 1; dt = [T/4]; 
M = 1E7; alpha = 0.95;

DeltaArithDirect = [];
MCStdArithDirect = [];
MCConfintArithDirect = [];

col = 0;
for i = 1
    for j = 1
        col = (j - 1) * 2;
        tic;
        [MCstd, MCAsianArithDelta, MCConfInt] = AsianArithMC(S0, K(i), sigma, r, T, dt(j), M, alpha);
        DeltaArithDirect(i,j) = MCAsianArithDelta;
        MCStdArithDirect(i,j) = MCstd;
        MCConfintArithDirect(i,col+1:col+2) = MCConfInt;
        toc
    end
end

% Create vectors for K and dt pairing
[K_mesh, dt_mesh] = meshgrid(K, dt);
K_vals = K_mesh(:);  % Flattened vector of K values
dt_vals = dt_mesh(:);  % Flattened vector of dt values

% Correct table for Direct simulation
DirectTable = table(K_vals, dt_vals, ...
    reshape(DeltaArithDirect', [], 1), reshape(MCStdArithDirect', [], 1), ...
    reshape(MCConfintArithDirect(:,1:2:end)', [], 1), reshape(MCConfintArithDirect(:,2:2:end)', [], 1), ...
    'VariableNames', {'K', 'dt', 'Delta_Direct', 'MCStd_Direct', 'ConfInt_Lower', 'ConfInt_Upper'});

% Displaying the corrected table
disp('Table for Direct Arithmetic Simulation:')
disp(DirectTable)

%------------------------------------------------------------------------------------------------------------------------------
function [MCstd, MCAsianArithDelta, MCConfInt] = AsianArithMC(S0, K, sigma, r, T, dt, M, alpha)
b = r - sigma^2/2;
n = T/dt; % number of monitoring dates
LRDeltaAn = zeros(1,M);

for j = 1:M
    S = [S0 zeros(1,n)];
    for i = 1:n
        S(i+1) = S(i)*exp((r - sigma^2/2)*dt + sigma*sqrt(dt)*randn);
    end
    zetaS1 = (log(S(2)/S(1))-b*dt)/(sigma*sqrt(dt));
    LRDeltaAn(j) = exp(-r*T)*max(mean(S(2:n+1))-K,0)*(zetaS1/(S0*sigma*sqrt(dt)));
end

MCAsianArithDelta = mean(LRDeltaAn);
MCstd = std(LRDeltaAn)/sqrt(M);
format long g
MCConfInt = MCAsianArithDelta + norminv(0.5+alpha/2)*MCstd*[-1 1];
end