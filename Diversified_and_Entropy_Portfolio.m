function [weights_m, portfolio_m_return, portfolio_m_std, portfolio_m_SR, ...
          weights_n, portfolio_n_return, portfolio_n_std, portfolio_n_SR] = ...
          Diversified_and_Entropy_Portfolio(mean_returns, cov_matrix, capitalizations, names, risk_free_rate, mkt, P2,num_assets)
% Compute the Maximum Diversified Portfolio (Portfolio M) and the
% Maximum Entropy (in asset volatility) Portfolio (Portfolio N), under
% the following constraints (to be considered all at once):
% - Standard constraints,
% - The total exposure on cyclicals has to be greater than 20%,
% - Assuming that you have a benchmark portfolio (capitalization
%   weighted portfolio), the sum of the difference (in absolute value)
%   of the weights in the benchmark portfolio and the optimal weights
%   has to be greater than 20%

cyclicalIdx = ismember(P2.AssetList, mkt.sector.cyclical);

% Weights of the benchmark portfolio (capitalization weighted portfolio)
cap_wghtd_ptf = capitalizations{1, names} / sum(capitalizations{1, names});

% Set up optimization problem
% Iniial guess
initial_guess = ones(num_assets, 1) / num_assets;
% Set default constraints for the portfolio: weights sum to 1
% and there are no short positions (non-negative weights)
lb = zeros(num_assets, 1);  % Lower bounds on weights
ub = ones(num_assets, 1);   % Upper bounds on weights
Aeq = ones(1, num_assets);  % Equality constraint for sum of weights
beq = 1;

% Additional constraints for sectors:
% Cyclical sectors total weight >= 20%
A = -cyclicalIdx;
b = -0.2;

% Additional constraint wrt the benchmark portfolio
nonlinconstr =  @(weights) customAbsDiffConstraint(weights, cap_wghtd_ptf); % Non linear constraint function
% Options
options = optimoptions('fmincon', ...     
        'Algorithm', 'sqp', ... % Specify the algorithm     
        'StepTolerance', 1e-6, ...       % Smaller than default StepTolerance     
        'Display', 'off');               % Show iteration information

% Portfolio M: Maximum Diversified Portfolio
diversification_ratio = @(w) -log(w' * sqrt(diag(cov_matrix)) / sqrt(w' * cov_matrix * w));

[weights_m, ptf_m_var] = fmincon(diversification_ratio, initial_guess, A, b, Aeq, beq, lb, ub, nonlinconstr, options);
portfolio_m_return = mean_returns' * weights_m;

% portfolio_m_std = sqrt(ptf_m_var);
portfolio_m_std = sqrt(weights_m' * cov_matrix * weights_m); %NEW

portfolio_m_SR = (portfolio_m_return - risk_free_rate) / portfolio_m_std;

% Portfolio N:  Maximum Entropy (in asset volatility) Portfolio
% entropy = @(w) sum(w.^2' * diag(cov_matrix)/ sum(w.^2' * diag(cov_matrix)) * log(w.^2' * diag(cov_matrix)/ sum(w.^2' * diag(cov_matrix))));
entropy = @(w) sum(w.^2 .* diag(cov_matrix)/ sum(w.^2 .* diag(cov_matrix)) .* log(w.^2 .* diag(cov_matrix)/ sum(w.^2 .* diag(cov_matrix))));%NEW Matte

[weights_n, ptf_n_var] = fmincon(entropy, initial_guess, A, b, Aeq, beq, lb, ub, nonlinconstr, options);
portfolio_n_return = mean_returns' * weights_n;

% portfolio_n_std = sqrt(ptf_n_var);
portfolio_n_std = sqrt(weights_n' * cov_matrix * weights_n); %NEW

portfolio_n_SR = (portfolio_n_return - risk_free_rate) / portfolio_n_std;

% Display Portfolio M - Maximum Diversified Portfolio
disp('==============================================================================================')
disp('Maximum Diversified Portfolio (Portfolio M)')
disp('==============================================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(weights_m)
    fprintf('%-25s %.4f\n', names{i}, weights_m(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', portfolio_m_return);
fprintf('%-25s %.4f\n', 'Volatility', portfolio_m_std);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', portfolio_m_SR);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(weights_m));
disp('              ')

% Display Portfolio N - Maximum Entropy Portfolio
disp('==============================================================================================')
disp('Maximum Entropy Portfolio (Portfolio N)')
disp('==============================================================================================')
disp('Asset Name                Weight')    
disp('-------------------------------------------')
for i = 1:length(weights_n)
    fprintf('%-25s %.4f\n', names{i}, weights_n(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', portfolio_n_return);
fprintf('%-25s %.4f\n', 'Volatility', portfolio_n_std);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', portfolio_n_SR);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(weights_n));
disp('              ')
