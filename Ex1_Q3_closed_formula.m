%Question 3

clear all
close all
clc

format short

%Parameter initialisation
k = [90; 100; 110]; %Vector k
n = [4; 12; 50]; %Vector n
S0 = 100;
r = 0.04;
sigm = 0.3;
T = 1;

rcol = [size(n,1), size(k,1)]; %Size of full matrix using n rows by k cols

%Creating initial matrices for the final values of all equations
finalExLBn = zeros(rcol(1),rcol(2));
finalFSenLBn = zeros(rcol(1),rcol(2));
finalExGn = zeros(rcol(1),rcol(2));
finalFSenGn = zeros(rcol(1),rcol(2));
%Part a & b
%Nested loop to run all the functions and find the table values for each k and n value
for y = 1:1:rcol(1) %rows organised in n values
    for x = 1:1:rcol(2) %columns organised in k values
        finalExLBn(y,x) = expLBn(S0, r, T, n(y), sigm, k(x)); %Matrix of expected LB
        finalFSenLBn(y,x) = fSensLBn(S0, r, T, n(y), sigm, k(x)); %Matrix of first sensitivity LB
        finalExGn(y,x) = expGn(S0, r, T, n(y), sigm, k(x)); %Matrix of expected geometric Asian call
        finalFSenGn(y,x) = fSensGn(S0, r, T, n(y), sigm, k(x)); %Matrix of first sensitivity Asian call
    end
end

%Print values for the final values of each equation
%Creating the row and column headers for the tables
rowNames = arrayfun(@(x) sprintf('n = %d', x), n, 'UniformOutput', false);
colNames = arrayfun(@(x) sprintf('K = %d', x), k, 'UniformOutput', false);

%Displaying all the various matrices for each equation9
disp('Expected Lower Bound Prices:')
lbTable = array2table(finalExLBn, "RowNames", rowNames, 'VariableNames', colNames)

disp('Lower Bound Deltas:')
lbDeltaTable = array2table(finalFSenLBn, "RowNames", rowNames, 'VariableNames', colNames)

disp('Expected Geometric Asian Call Option Price:')
gnTable = array2table(finalExGn, "RowNames", rowNames, 'VariableNames', colNames)

disp('Geometric Asian Call Option Delta:')
gnDeltaTable = array2table(finalFSenGn, "RowNames", rowNames, 'VariableNames', colNames)


%All of the following are functions to the equations of LB, geo Asian
%option and their first sensitivities
%Function to calculate the expected lower bound
function ans = expLBn(S0, r, T, n, sigm, k)
delta = T/n;

sigBar = sigm*sqrt((2*n+1)/(3*n));
tBar = (n+1)*delta/2;

b = (log(S0/k)+(r-sigm^2/2)*tBar)/(sigBar*sqrt(tBar));

%Initialisation of the summation within the differential equation
total = 0;
for m = 1:1:n %Loop to calculate the sums within the final equation
    %m represents the current step within the sum
    muK = (r - sigm^2/2)*m*delta;
    sigK = sigm*sqrt(m*delta);
    aK = sigm*sqrt(delta)*(m*(n+1-(m+1)/2))/(sqrt(n*(n+1)*(2*n+1)/6));

    total = total + exp(muK+sigK^2/2)*normcdf(b+aK); %Total represents the sum of e^(mu+sigm^2/2)*N(b+a) within the equation
end

ans = (S0*exp(-r*T)/n)*total-k*exp(-r*T)*normcdf(b);
end


%Function to calculate the first sensitivity of the lower bound
function ans = fSensLBn(S0, r, T, n, sigm, k)
delta = T/n; %Delta being 

sigBar = sigm*sqrt((2*n+1)/(3*n)); %Sigma bar
tBar = (n+1)*delta/2; %T bar

b = (log(S0/k)+(r-sigm^2/2)*tBar)/(sigBar*sqrt(tBar)); %b

%Initialisation of the summation within the differential equation
total = 0;
for m = 1:1:n %Loop to calculate the sums within the final equation
    %m represents the current step within the sum
    muK = (r - sigm^2/2)*m*delta;
    sigK = sigm*sqrt(m*delta);
    aK = sigm*sqrt(delta)*(m*(n+1-(m+1)/2))/(sqrt(n*(n+1)*(2*n+1)/6));

    total = total + exp(muK+sigK^2/2)*(normcdf(b+aK)+exp(-(b+aK)^2/2)*(1/(sigBar*sqrt(tBar*2*pi))));
end

ans = exp(-r*T)*((1/n)*total-k*exp(-(b)^2/2)*(1/(S0*sigBar*sqrt(tBar*2*pi)))); %Final equation of delta of LB
end


%Function to calculate the expected geometric Asian call
function ans = expGn(S0, r, T, n, sigm, k)
delta = T/n;

sigBar = sigm*sqrt((2*n+1)/(3*n));
tBar = (n+1)*delta/2;

d = (log(S0/k)+(r-sigm^2/2+sigBar^2)*tBar)/(sigBar*sqrt(tBar));

ans = S0*exp((r-sigm^2/2+sigBar^2/2)*tBar-r*T)*normcdf(d)-k*exp(-r*T)*normcdf(d-sigBar*sqrt(tBar));
end


%Function to calculate the first sensitivity of the geometric Asian call
function ans = fSensGn(S0, r, T, n, sigm, k)
delta = T/n;

sigBar = sigm*sqrt((2*n+1)/(3*n));
tBar = (n+1)*delta/2;

d = (log(S0/k)+(r-sigm^2/2+sigBar^2)*tBar)/(sigBar*sqrt(tBar));

ans = exp(-r*T)*(exp((r-sigm^2/2+sigBar^2/2)*tBar)*(normcdf(d)+exp(-(d)^2/2)*(1/(sigBar*sqrt(tBar*2*pi))))-k*exp(-(d-sigBar*sqrt(tBar))^2/2)*(1/(S0*sigBar*sqrt(tBar*2*pi))));
end