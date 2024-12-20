clc
clear all
close all
warning on
% Load data (Assumes you have 'capitalizations.xlsx' and 'prices.xlsx')
capitalizations = readtable('capitalizations.xlsx', 'Sheet', 'Sheet1');
caps = capitalizations{1,2:end};
prices = readtable('prices.xlsx', 'Sheet', 'Sheet1');
% Extract prices and calculate daily returns
dates = prices{:, 1};
prices_data = prices{:, 2:end};
start_date = datetime(2023,1,1);
end_date = datetime(2023,12,31);
prices_2023 = prices_data(dates >= start_date & dates <= end_date, :);
prices_2024 = prices_data(dates > end_date, :);
returns = tick2ret(prices_2023,'Method','continuous');
% Compute expected returns (mean of daily returns)
mean_ret = mean(returns);
% Compute the variance-covariance matrix
cov_matrix = cov(returns);
% Extract capitalizations and calculate initial weights
capitalization_values = capitalizations{1, 2:end}';
initial_weights = capitalization_values / sum(capitalization_values);
% Risk-free rate (annualized to daily)
risk_free_rate = 0;
% Confidence level for VaR
confidence_level = .95;
VaR_N = @(x) - (mean_ret * x - ...
    sqrt(x' * cov_matrix * x) * norminv(confidence_level));
% Define the objective function (negative VaR-modified Sharpe Ratio)
objective_function_N = @(x) -(mean_ret * x - risk_free_rate) / VaR_N(x);
% Constraints: weights sum to 1 and are between 0 and 1
Aeq = ones(1, length(initial_weights));
beq = 1;
lb = zeros(length(initial_weights), 1);
ub = ones(length(initial_weights), 1);
% Optimization using fmincon
options = optimoptions('fmincon', 'Display', 'None');
w_N = fmincon(objective_function_N, ...
              initial_weights, [], [], ...
              Aeq, beq, lb, ub, [], options);
% Asset names (assuming they are the column headers in the prices file)
names = prices.Properties.VariableNames(2:end);
ptfs = [w_N ones(16,1)/16 caps'/sum(caps)];
% Display results in the specified format
disp('===========================================================================')
disp(' VaR modified sharpe ratio Portfolio (Portfolio Q)')
disp('===========================================================================')
weightsTable = array2table(round(ptfs,4), "RowNames",names,...
    "VariableNames",{'95_N','EW','Caps'})
[~, table] = getEquityandMetrices(weightsTable, prices_2023, "2023")
[~, table] = getEquityandMetrices(weightsTable, prices_2024, "2024")

