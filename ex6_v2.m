clc;
clear;
close all;

% Load data from Excel files
prices = readtable('prices.xlsx');
capitalizations = readtable('capitalizations.xlsx');
names = capitalizations.Properties.VariableNames(2:end);

% Extract dates and price data
dates = datetime(prices{:, 1}); % Convert to MATLAB datetime
prices_data = prices{:, 2:end};

% Filter prices for the year 2023
start_date = datetime(2023, 1, 1);
end_date = datetime(2023, 12, 31);
prices_2023 = prices_data(dates >= start_date & dates <= end_date, :);

% Define risk-free rate
risk_free_rate = 0.04 / 365; % Daily risk-free rate

% Calculate daily returns
returns_2023 = diff(log(prices_2023));

%% Standardize returns for PCA
returns_std = zscore(returns_2023);

% Perform PCA
[coeff, score, latent, ~, explained] = pca(returns_std);

% Select components explaining >90% variance
cumulative_variance = cumsum(explained);
num_components = find(cumulative_variance > 90, 1);

% Display PCA results
fprintf('Selected Components: %d\n', num_components);
fprintf('Cumulative Variance Explained: %.2f%%\n', cumulative_variance(num_components));

% Reduce data to selected components
reduced_data = score(:, 1:num_components);
mean_returns = mean(reduced_data, 1)';
cov_matrix = cov(reduced_data);


% Define target volatility
target_volatility = 0.1;

% Initial weights (equal weights)
initial_weights = ones(num_components, 1) / num_components;

% Bounds: weights between 0 and 1 (long-only portfolio)
lb = zeros(num_components, 1);
ub = ones(num_components, 1);

% Equality constraint: weights sum to 1
Aeq = ones(1, num_components);
beq = 1;

% Non-linear constraint: portfolio volatility <= target volatility
nonlincon = @(w) deal([], sqrt(w' * cov_matrix * w) - target_volatility);

% Objective function: minimize negative expected return
objective = @(w) -w' * mean_returns;

% Solve optimization problem using fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
[weights_opt, fval, exitflag, output] = fmincon(objective, initial_weights, [], [], Aeq, beq, lb, ub, nonlincon, options);

% Results
expected_return_opt = weights_opt' * mean_returns;
volatility_opt = sqrt(weights_opt' * cov_matrix * weights_opt);
sharpe_ratio_opt = (expected_return_opt - risk_free_rate) / volatility_opt;

% Display Results
disp('===========================================================================')
disp('    Maximum expected return with PCA  (Portfolio P)   ')
disp('===========================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:num_components
    fprintf('%-25s %.4f\n', names{i}, weights_opt(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', expected_return_opt);
fprintf('%-25s %.4f\n', 'Volatility', volatility_opt);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', sharpe_ratio_opt);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(weights_opt));
disp('  ')
