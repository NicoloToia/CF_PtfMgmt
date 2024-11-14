clc
clear all
close all
warning off

%% PART A 

% Load data from Excel files
prices = readtable('prices.xlsx');
capitalizations = readtable('capitalizations.xlsx');

% Extract dates and data
dates = prices{:,1}; % First column contains dates
prices_data = prices{:,2:end}; % Data starts from the second column

% Convert dates to MATLAB date format if needed
dates = datetime(dates);

% Filter prices for 2023
start_date = datetime(2023,1,1);
end_date = datetime(2023,12,31);
prices_2023 = prices_data(dates >= start_date & dates <= end_date, :);

% Names
assetNames = capitalizations.Properties.VariableNames(2:end);

% Calculate daily returns for each index in 2023
returns_2023 = diff(log(prices_2023));

% Calculate mean returns and covariance matrix
mean_returns = mean(returns_2023)';
cov_matrix = cov(returns_2023);

% Define number of assets
num_assets = length(mean_returns);

% Define the risk-free rate (annual, assuming 4% risk-free rate)
risk_free_rate = 0.04 / 365; % Convert to daily (365 is correct)

%% Portfolio object
% p = Portfolio('AssetList', assetNames);
% p = setDefaultConstraints(p);
% P = setAssetMoments(p, mean_returns, cov_matrix);
% pwgt = estimateFrontier(P, 100);
% [pf_risk, pf_Retn] = estimatePortMoments(P, pwgt);
% plot(pf_risk, pf_Retn)
% pf_Retn(pf_risk == min(pf_risk))
% min(pf_risk)
% max(pf_risk)
% 
% (mean_returns' * pwgt(:,1));
% sum(pwgt(:,1));
% pwgt(:,1);
% 
% sr = (pf_Retn - risk_free_rate) ./ pf_risk;
% max(sr)

% Initialize a Portfolio object with the list of asset names
p = Portfolio('AssetList', assetNames);

% Set default constraints for the portfolio: weights sum to 1
% and there are no short positions (non-negative weights)
p = setDefaultConstraints(p);

% Set the expected returns and covariance matrix for the assets
% to define the portfolio's risk and return characteristics
P = setAssetMoments(p, mean_returns, cov_matrix);

% Estimate the efficient frontier with 100 portfolios
pwgt = estimateFrontier(P, 100);

% Calculate the risk (standard deviation) and expected returns for each portfolio on the frontier
[pf_risk, pf_Retn] = estimatePortMoments(P, pwgt);

% Plot the efficient frontier, which shows the risk-return trade-off of optimized portfolios
figure()
plot(pf_risk, pf_Retn)
xlabel('Volatility')
ylabel('Return')
legend('Portfolio froniter with object Portfolio')

% Display the return of the portfolio with the minimum risk on the frontier
ptf_a_return_viaobj = pf_Retn(pf_risk == min(pf_risk));

% Display the minimum risk (volatility) of the portfolios on the frontier
ptf_a_vol_viaobj = min(pf_risk);

% Display the maximum risk (standard deviation) of the portfolios on the frontier
% max(pf_risk)

% Calculate the return of the first portfolio on the frontier
% frontier_returns = (mean_returns' * pwgt(:,:));

% Check that the weights of the first portfolio sum to 1, confirming a fully invested portfolio
if sum(pwgt(:,1))>1-1e-8 && sum(pwgt(:,1))<1+1e-8 && min(pwgt(:,1))>=0
    disp('Standard constraints satisfied by Portfolio A')
else
    disp('Standard constraints NOT satisfied by Portfolio A')
end

% Display the weights of the first portfolio on the frontier
weights_a_viaobj = pwgt(:,1);

% Calculate the Sharpe Ratio for each portfolio on the frontier
sr = (pf_Retn - risk_free_rate) ./ pf_risk;

% Display the maximum Sharpe Ratio, representing the portfolio with the best risk-adjusted return
ptf_b_sharpe_viaobj = max(sr);

% Display results
disp('Portfolio A (Minimum Variance Portfolio) USING THE OBJECT PORTFOLIO:');
disp(['Expected Return: ', num2str(ptf_a_return_viaobj)]);
disp(['Volatility: ', num2str(ptf_a_vol_viaobj)]);
disp(['Sharpe Ratio: ', num2str(sr(pf_risk == min(pf_risk)))]);
disp('Weights:');
disp(weights_a_viaobj');

disp('Portfolio B (Maximum Sharpe Ratio Portfolio):');
disp(['Expected Return: ', num2str(pf_Retn(sr==ptf_b_sharpe_viaobj))]);
disp(['Volatility: ', num2str(pf_risk(sr==ptf_b_sharpe_viaobj))]);
disp(['Sharpe Ratio: ', num2str(ptf_b_sharpe_viaobj)]);
disp('Weights:');
disp(pwgt(:,find((sr==ptf_b_sharpe_viaobj)))');

%% 1) 

% Set up optimization problem
lb = zeros(num_assets, 1);  % Lower bounds on weights
ub = ones(num_assets, 1);   % Upper bounds on weights
Aeq = ones(1, num_assets);  % Equality constraint for sum of weights
beq = 1;

% 1. Portfolio A: Minimum Variance Portfolio
options = optimoptions('fmincon', ...     
        'Algorithm', 'sqp', ... % Specify the algorithm     
        'StepTolerance', 1e-6, ...       % Smaller than default StepTolerance     
        'Display', 'off');               % Show iteration information
min_variance_func = @(w) w' * cov_matrix * w;
initial_guess = ones(num_assets, 1) / num_assets;


[weights_a, min_var] = fmincon(min_variance_func, initial_guess, [], [], Aeq, beq, lb, ub, [], options);
portfolio_a_return = mean_returns' * weights_a;
portfolio_a_std = sqrt(min_var);
portfolio_a_sharpe = (portfolio_a_return - risk_free_rate) / portfolio_a_std;

% 2. Portfolio B: Maximum Sharpe Ratio Portfolio
sharpe_ratio_func = @(w) -((mean_returns' * w - risk_free_rate) / sqrt(w' * cov_matrix * w));

[weights_b, ~] = fmincon(sharpe_ratio_func, initial_guess, [], [], Aeq, beq, lb, ub, [], options);
portfolio_b_return = mean_returns' * weights_b;
portfolio_b_std = sqrt(weights_b' * cov_matrix * weights_b);
portfolio_b_sharpe = (portfolio_b_return - risk_free_rate) / portfolio_b_std;

% Display results
disp('Portfolio A (Minimum Variance Portfolio):');
disp(['Expected Return: ', num2str(portfolio_a_return)]);
disp(['Volatility: ', num2str(portfolio_a_std)]);
disp(['Sharpe Ratio: ', num2str(portfolio_a_sharpe)]);
disp('Weights:');
disp(weights_a');

disp('Portfolio B (Maximum Sharpe Ratio Portfolio):');
disp(['Expected Return: ', num2str(portfolio_b_return)]);
disp(['Volatility: ', num2str(portfolio_b_std)]);
disp(['Sharpe Ratio: ', num2str(portfolio_b_sharpe)]);
disp('Weights:');
disp(weights_b');

%% 2) 



% Define constraints based on problem statement
lb = zeros(num_assets, 1);  % Lower bounds on weights
ub = ones(num_assets, 1);   % Upper bounds on weights
Aeq = ones(1, num_assets);  % Equality constraint for sum of weights
beq = 1;

% Set Consumer Staples and Low Volatility weights to 0
% Assuming Consumer Staples is asset #8 and Low Volatility is asset #16
lb(8) = 0; ub(8) = 0; % Consumer Staples
lb(16) = 0; ub(16) = 0; % Low Volatility

% Additional constraints for sectors:
% Sensible sectors: e.g., Health Care, Utilities, Consumer Staples
% Cyclical sectors: e.g., Consumer Discretionary, Industrials, Energy

SectorConstraints = zeros(11,16);
for i =1:11
    SectorConstraints(i,i) = 1;
end
SectorConstraints = zeros(1,16);
SectorConstraints(1:11) = 1;
Aineq = [
    % Sensible sectors total weight >= 10%
    -[1 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0];
    
    % Cyclical sectors total weight <= 30%
    [0 1 0 1 0 1 0 0 0 1 1 0 0 0 0 0];
    % Consumer staple at 0
    [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
    % Low volatility at 0
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    % Maximum exposure in any sector ≤ 80%
    SectorConstraints;
];

bineq = [
    -0.1;    % Sensible sectors >= 10%
    0.3;     % Cyclical sectors <= 30%
    0;
    0;
    %0.8 * ones(11, 1)  % Max exposure for each sector
    0.8
];
%% 2. su oggetto
P.AInequality = Aineq;
P.bInequality = bineq;
P.AEquality = Aeq;
P.bEquality = beq;
pwgt = estimateFrontier(P, 100);
[pf_risk, pf_Retn] = estimatePortMoments(P, pwgt);
plot(pf_risk, pf_Retn)
pf_Retn(pf_risk == min(pf_risk))
disp('Min Var')
min(pf_risk)
(mean_returns' * pwgt(:,1));
sum(pwgt(:,1));
disp('Min Variance Port')
pwgt(:,pf_risk == min(pf_risk))

sr = (pf_Retn - risk_free_rate) ./ pf_risk;
disp('Max Sr')
max(sr)
disp('Port SR')
pwgt(:, sr == max(sr))
%%




% Optimization options
options = optimoptions('fmincon', ...     
        'Algorithm', 'sqp', ...         % Specify the algorithm     
        'StepTolerance', 1e-12, ...      % Smaller than default StepTolerance     
        'Display', 'off');              % Show iteration information

% 1. Portfolio C: Minimum Variance Portfolio with constraints
min_variance_func = @(w) w' * cov_matrix * w;
initial_guess = 0.5 * ones(num_assets, 1) / num_assets;

[weights_c, min_var] = fmincon(min_variance_func, initial_guess, Aineq, bineq, Aeq, beq, lb, ub, [], options);
portfolio_c_return = mean_returns' * weights_c;
portfolio_c_std = sqrt(min_var);
portfolio_c_sharpe = (portfolio_c_return - risk_free_rate) / portfolio_c_std;

% 2. Portfolio D: Maximum Sharpe Ratio Portfolio with constraints
sharpe_ratio_func = @(w) -((mean_returns' * w - risk_free_rate) / sqrt(w' * cov_matrix * w));

[weights_d, ~] = fmincon(sharpe_ratio_func, initial_guess, Aineq, bineq, Aeq, beq, lb, ub, [], options);
portfolio_d_return = mean_returns' * weights_d;
portfolio_d_std = sqrt(weights_d' * cov_matrix * weights_d);
portfolio_d_sharpe = (portfolio_d_return - risk_free_rate) / portfolio_d_std;

% Display results
disp('Portfolio C (Minimum Variance Portfolio):');
disp(['Expected Return: ', num2str(portfolio_c_return)]);
disp(['Volatility: ', num2str(portfolio_c_std)]);
disp(['Sharpe Ratio: ', num2str(portfolio_c_sharpe)]);
disp('Weights:');
disp(weights_c');

disp('Portfolio D (Maximum Sharpe Ratio Portfolio):');
disp(['Expected Return: ', num2str(portfolio_d_return)]);
disp(['Volatility: ', num2str(portfolio_d_std)]);
disp(['Sharpe Ratio: ', num2str(portfolio_d_sharpe)]);
disp('Weights:');
disp(weights_d');
