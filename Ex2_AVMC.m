% Monte Carlo simulation with antithetic variates to estimate the price of the discretely monitored UOC barrier option 
clear all % clears memory
close all
clc

S0 = 110; K = 100;
UP_Barrier = [160 170]; % barrier level
sigma = 0.1; r = 0.05;
T = 2;
n_vec = 2.^(2:10); % option is path-dependent: n gives the no. of monitoring dates
M = 1E6; % no. of simulations

num_n = length(n_vec);
MCPrices = zeros(num_n, length(UP_Barrier));
MCstd = zeros(num_n, length(UP_Barrier));
MCConfInt = zeros(num_n, length(UP_Barrier),2);
AVPrices = zeros(num_n, length(UP_Barrier));
AVstd = zeros(num_n, length(UP_Barrier));
AVConfInt = zeros(num_n, length(UP_Barrier),2);

tic
for idx_n = 1:num_n
    n = n_vec(idx_n);
    dt = T/n;

    for idx_U = 1:length(UP_Barrier)
        U = UP_Barrier(idx_U);
       
        UO = zeros(1,M);
        UOAntVar = zeros(1,M);
        
        for j = 1:M
            S = [S0 zeros(1,n)]; % pre-allocate array to dimension n+1 by filling with zeros
            SAntVar = [S0 zeros(1,n)];
            Z = randn(1,n);

            for i = 1:n
            S(i+1) = S(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*Z(i)); % exact iterative simulation of the geoBM
            SAntVar(i+1) = SAntVar(i)*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*(-Z(i))); % ANTITHETIC stock price trajectory 
            end
           
            barrier_breached = max(S)>=U;
            terminal_payoff= max(S(end)-K,0);
            
            if ~barrier_breached %this flip the logic order,if the barrier is not breached, then calculate the payoff
               UO(j)=exp(-r*T)*terminal_payoff;
            else % if the barrier is breached, then the payoff will be zero
               UO(j)=0;
            end
            
            AV_barrier_breached = max(SAntVar)>=U;
            AV_terminal_payoff= max(SAntVar(end)-K,0);
            
            if ~AV_barrier_breached %this flip the logic order,if the barrier is not breached, then calculate the payoff
               UOAntVar(j)=exp(-r*T)*AV_terminal_payoff;
            else % if the barrier is breached, then the payoff will be zero
               UOAntVar(j)=0;
            end
        end
        
        r3 = corrcoef(UO,UOAntVar);
        AVPrice = mean((UO+UOAntVar)/2); %Option price (AV)
        MCPrice = mean(UO); %Option price (MC)
        AVStdVal = std((UO+UOAntVar)/2)/sqrt(M); % standard errors (AV)
        MCStdVal = std(UO)/sqrt(M); % standard errors (MC)

        alpha = 0.95; % confidence interval parameter
        AVCI = AVPrice + norminv(0.5+alpha/2)*AVStdVal*[-1 1]; % 100*alpha% confidence interval (AV)
        MCCI = MCPrice + norminv(0.5+alpha/2)*MCStdVal*[-1 1]; % 100*alpha% confidence interval (AV)

        AVPrices(idx_n, idx_U)= AVPrice;
        AVstd(idx_n, idx_U) = AVStdVal;
        AVConfInt(idx_n, idx_U,:) = AVCI;

        MCPrices(idx_n, idx_U)= MCPrice;
        MCstd(idx_n, idx_U) = MCStdVal;
        MCConfInt(idx_n, idx_U,:) = MCCI;
    end
end
exectime1 = toc
format short g

% This is to calculate the price of the continuously monitored up-and-out barrier call option
closed_160 = up_and_out_call(S0, K, 160, sigma, r, T)
closed_170 = up_and_out_call(S0, K, 170, sigma, r, T)

function closed_price = up_and_out_call(S0, K, U, sigma, r, T) % The given closed-form formula 
    

    N = @(x) normcdf(x);
    

    d_plus = @(x) (log(x) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T));
    d_minus = @(x) (log(x) + (r - 0.5 * sigma^2) * T) / (sigma * sqrt(T));
    

    term1 = S0 * (N(d_plus(S0 / K)) - N(d_plus(S0 / U)));
    term2 = -K * exp(-r * T) * (N(d_minus(S0 / K)) - N(d_minus(S0 / U)));
    
    term3 = -U * (S0 / U)^(-2 * r / sigma^2) * (N(d_plus(U^2 / (K * S0))) - N(d_plus(U / S0)));
    term4 = K * exp(-r * T) * (S0 / U)^(1 - 2 * r / sigma^2) * (N(d_minus(U^2 / (K * S0))) - N(d_minus(U / S0)));
    

    closed_price = term1 + term2 + term3 + term4;
end

% Below is to display the simulation result
for idx_U = 1:length(UP_Barrier)
    n_column = n_vec(:);
    MCprice_column = MCPrices(:, idx_U);
    MCstd_column = MCstd(:, idx_U);
    MCCI_lower = MCConfInt(:, idx_U, 1);
    MCCI_upper = MCConfInt(:, idx_U, 2);
    
    AVprice_column = AVPrices(:, idx_U);
    AVstd_column = AVstd(:, idx_U);
    AVCI_lower = AVConfInt(:, idx_U, 1);
    AVCI_upper = AVConfInt(:, idx_U, 2);

    % Create a table for this barrier level
    ResultsTable = table(n_column, MCprice_column, MCstd_column, MCCI_lower, MCCI_upper, AVprice_column, AVstd_column, AVCI_lower, AVCI_upper, ...
        'VariableNames', {'n', 'MCPrice', 'MCStdError', 'MCCI_Lower', 'MCCI_Upper', 'AVPrice', 'AVStdError', 'AVCI_Lower', 'AVCI_Upper'});

    disp(['Results for Barrier U = ', num2str(UP_Barrier(idx_U))]);
    disp(ResultsTable);
end


figure; % for Convergence plot
hold on;

% Plot prices for U=160 and U=170
h1 = semilogx(n_vec, AVPrices(:, 1), 'o-', 'LineWidth', 1.5, 'DisplayName', 'Discrete U=160');
h2 = semilogx(n_vec, AVPrices(:, 2), 's-', 'LineWidth', 1.5, 'DisplayName', 'Discrete U=170');

% Plot continuous-monitoring prices as horizontal lines
h3 = yline(closed_160, '--r', 'Continuous U=160', 'LineWidth', 1.5);
h4 = yline(closed_170, '--b', 'Continuous U=170', 'LineWidth', 1.5);

xlabel('Number of Monitoring Dates (n)');
ylabel('Option Price');
title('Convergence of Discretely Monitored UOC Prices');
legend([h1 h2 h3 h4], {'Discrete U=160', 'Discrete U=170', 'Continuous U=160', 'Continuous U=170'}, 'Location','best');
grid on;
hold off

% Plot the payoff correlation
figure;

scatter(UO, UOAntVar, 10, 'b', 'o'); % 10 is the marker size
xlabel('MC Payoff');
ylabel('Antithetic Payoff');
title('Antithetic option payoff samples');
grid on;