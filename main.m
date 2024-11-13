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

% Calculate daily returns for each index in 2023
returns_2023 = diff(log(prices_2023));

% Calculate mean returns and covariance matrix
mean_returns = mean(returns_2023)';
cov_matrix = cov(returns_2023);

% Define number of assets
num_assets = length(mean_returns);

% Define the risk-free rate (annual, assuming 4% risk-free rate)
risk_free_rate = 0.04 / 252; % Convert to daily

%% 1) 

% Set up optimization problem
lb = zeros(num_assets, 1);  % Lower bounds on weights
ub = ones(num_assets, 1);   % Upper bounds on weights
Aeq = ones(1, num_assets);  % Equality constraint for sum of weights
beq = 1;

% 1. Portfolio A: Minimum Variance Portfolio
options = optimoptions('fmincon', 'Display', 'off');
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

Aineq = [
    % Sensible sectors total weight >= 10%
    -[0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0];
    
    % Cyclical sectors total weight <= 30%
    [0 0 0 0 1 0 1 0 0 0 0 0 1 0 0 0];
    
    % Maximum exposure in any sector ≤ 80%
    eye(num_assets)
];

bineq = [
    -0.1;    % Sensible sectors >= 10%
    0.3;     % Cyclical sectors <= 30%
    0.8 * ones(num_assets, 1)  % Max exposure for each sector
];

% Optimization options
options = optimoptions('fmincon', 'Display', 'off');

% 1. Portfolio C: Minimum Variance Portfolio with constraints
min_variance_func = @(w) w' * cov_matrix * w;
initial_guess = ones(num_assets, 1) / num_assets;

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
