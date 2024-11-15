clc;
clear;
close all;

% Load data from Excel files
prices = readtable('prices.xlsx');
capitalizations = readtable('capitalizations.xlsx');

% Extract dates and price data
dates = datetime(prices{:, 1}); % Convert to MATLAB datetime
prices_data = prices{:, 2:end};

% Filter prices for the year 2023
start_date = datetime(2023, 1, 1);
end_date = datetime(2023, 12, 31);
prices_2023 = prices_data(dates >= start_date & dates <= end_date, :);


% Calculate daily returns
returns_2023 = diff(log(prices_2023));

% Standardize returns for PCA
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
reduced_mean_returns = mean(reduced_data, 1)';
reduced_cov_matrix = cov(reduced_data);


% Define target volatility
target_volatility = 0.1;

% Check achievable volatility range
min_vol = sqrt(min(eig(reduced_cov_matrix)));
max_vol = sqrt(max(eig(reduced_cov_matrix)));
fprintf('Min Volatility: %.4f, Max Volatility: %.4f\n', min_vol, max_vol);

% Adjust target volatility if necessary
if target_volatility < min_vol
    fprintf('Adjusting target volatility to %.4f (minimum achievable).\n', min_vol);
    target_volatility = min_vol;
elseif target_volatility > max_vol
    fprintf('Adjusting target volatility to %.4f (maximum achievable).\n', max_vol);
    target_volatility = max_vol;
end

% Create and configure portfolio object
p = Portfolio();
p = setAssetMoments(p, reduced_mean_returns, reduced_cov_matrix);
p = setDefaultConstraints(p);


% Optimize portfolio
[weights_pf_P, pf_risk_P, Exp_return_pf_P] = estimateFrontierByRisk(p, target_volatility);

Exp_return_pf_P = Exp_return_pf_P' * weights_pf_P;

pf_risk_P = pf_risk_P'*weights_pf_P;

% Compute Sharpe Ratio
risk_free_rate = 0.04 / 365; % Daily risk-free rate
sr_pf_P = (Exp_return_pf_P - risk_free_rate) / pf_risk_P;


% Display header
disp('===========================================================================');
disp('Optimized Portfolio with Constraints');
disp('===========================================================================');

% Create a table for portfolio weights
asset_table = table((1:length(weights_pf_P))', weights_pf_P, ...
    'VariableNames', {'Principal_Component', 'Weight'});

% Display the table
disp(asset_table);

% Display performance metrics
disp('-------------------------------------------');
fprintf('%-25s %.4f\n', 'Expected Return', Exp_return_pf_P);
fprintf('%-25s %.4f\n', 'Volatility', pf_risk_P);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', sr_pf_P);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(weights_pf_P));
disp('===========================================================================');