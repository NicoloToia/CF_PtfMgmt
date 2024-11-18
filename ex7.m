clc
clear all
close all
warning off
% Load data (Assumes you have 'capitalizations.xlsx' and 'prices.xlsx')
capitalizations = readtable('capitalizations.xlsx', 'Sheet', 'Sheet1');
prices = readtable('prices.xlsx', 'Sheet', 'Sheet1');

% Extract prices and calculate daily returns
dates = prices{:, 1};
prices_data = prices{:, 2:end};
returns = diff(log(prices_data));

% Compute expected returns (mean of daily returns)
expected_returns = mean(returns)';

% Compute the variance-covariance matrix
var_cov_matrix = cov(returns);

% Extract capitalizations and calculate initial weights
capitalization_values = capitalizations{1, 2:end}';
initial_weights = capitalization_values / sum(capitalization_values);

% Risk-free rate (annualized to daily)
risk_free_rate = 0.04 / 365;

% Confidence level for VaR
confidence_level = 0.95;
z = norminv(1 - confidence_level);

% Define the objective function (negative VaR-modified Sharpe Ratio)
objective_function = @(weights) -((weights' * expected_returns - risk_free_rate) / ...
                                   (-z * sqrt(weights' * var_cov_matrix * weights)));

% Constraints: weights sum to 1 and are between 0 and 1
Aeq = ones(1, length(initial_weights));
beq = 1;
lb = zeros(length(initial_weights), 1);
ub = ones(length(initial_weights), 1);

% Optimization using fmincon
options = optimoptions('fmincon', 'Display', 'iter');
[optimized_weights, optimized_sharpe_ratio] = fmincon(objective_function, ...
                                                      initial_weights, [], [], ...
                                                      Aeq, beq, lb, ub, [], options);

% After running the optimization script and obtaining the results

% Assuming 'optimized_weights', 'expected_returns', and 'var_cov_matrix' are available
% Recalculate the portfolio metrics for visualization
portfolio_return = optimized_weights' * expected_returns;
portfolio_variance = optimized_weights' * var_cov_matrix * optimized_weights;
portfolio_volatility = sqrt(portfolio_variance);
sharpe_ratio = (portfolio_return - risk_free_rate) / ...
               (-z * portfolio_volatility);

% Asset names (assuming they are the column headers in the prices file)
names = prices.Properties.VariableNames(2:end);

% Display results in the specified format
disp('===========================================================================')
disp(' VaR modified sharpe ratio Portfolio (Portfolio Q)')
disp('===========================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(optimized_weights)
    fprintf('%-25s %.4f\n', names{i}, optimized_weights(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', portfolio_return);
fprintf('%-25s %.4f\n', 'Volatility', portfolio_volatility);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', sharpe_ratio);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(optimized_weights));
disp('  ')

