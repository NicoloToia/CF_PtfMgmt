clc
clear all
close all
warning off

%% PART A 

% Load data from Excel files
prices = readtable('prices.xlsx');
capitalizations = readtable('capitalizations.xlsx');

% Names
names = capitalizations.Properties.VariableNames(2:end);

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

% cyclical = ["ConsumerDiscretionary", "Financials", "Materials", "RealEstate", "Industrials"];
% defensive = ["ConsumerStaples", "Utilities", "HealthCare"];
% sensible = ["Energy", "InformationTechnology", "CommunicationServices"];

% factor = ["Momentum","Value","Growth","Quality","LowVolatility"];

% groups = {cyclical, defensive, sensible, factor};

mean_returns = mean(returns_2023)';
cov_matrix = cov(returns_2023);

% Define number of assets
num_assets = length(mean_returns);

% Define the risk-free rate (annual, assuming 4% risk-free rate)
risk_free_rate = 0.04 / 365; % Convert to daily (365 is correct)
% DA AGGIUNGERE DISPLAY DI EDO
%% 1. Efficient Frontier

% Create a portfolio object
P1 = Portfolio('AssetList', names);
P1 = setDefaultConstraints(P1); % all weights sum to 1, no shorting, and 100% investment in risky assets
P1 = setAssetMoments(P1, mean_returns, cov_matrix); % set mean returns and covariance matrix
% estimate efficient frontier, via object method and 100 P1oints
pwgt1 = estimateFrontier(P1, 100);
% estimate portfolio moments
[pf_risk_P1, pf_Retn_P1] = estimatePortMoments(P1, pwgt1);
% plot efficient frontier
% plot(pf_risk_P1, pf_Retn_P1)

% find minimum risk portfolio (Portfolio A - Minimum Variance Portfolio)
minRisk_P1 = min(pf_risk_P1);
minRiskWgt_P1 = pwgt1(:, pf_risk_P1 == minRisk_P1);
minRiskRet_P1 = mean_returns' * minRiskWgt_P1;
minRiskSR_P1 = (minRiskRet_P1 - risk_free_rate) / minRisk_P1;

% find maximum Sharpe ratio portfolio (Portfolio B - Maximum Sharpe Ratio Portfolio)
maxSharpeWgt_P1 = estimateMaxSharpeRatio(P1);
maxSharpeRet_P1 = mean_returns' * maxSharpeWgt_P1;
maxSharpeRisk_P1 = sqrt(maxSharpeWgt_P1' * cov_matrix * maxSharpeWgt_P1);
maxSharpeSR_P1 = (maxSharpeRet_P1 - risk_free_rate) / maxSharpeRisk_P1;

% Display Portfolio A - Minimum Variance Portfolio
disp('===========================================================================')
disp('    Minimum Risk Portfolio (Portfolio A)   ')
disp('===========================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(minRiskWgt_P1)
    fprintf('%-25s %.4f\n', names{i}, minRiskWgt_P1(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', minRiskRet_P1);
fprintf('%-25s %.4f\n', 'Volatility', minRisk_P1);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', minRiskSR_P1);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(minRiskWgt_P1));
disp('  ')

% Display Portfolio B - Maximum Sharpe Ratio Portfolio
disp('===========================================================================')
disp(' Maximum Sharpe Ratio Portfolio (Portfolio B)')
disp('===========================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(maxSharpeWgt_P1)
    fprintf('%-25s %.4f\n', names{i}, maxSharpeWgt_P1(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', maxSharpeRet_P1);
fprintf('%-25s %.4f\n', 'Volatility', maxSharpeRisk_P1);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', maxSharpeSR_P1);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(maxSharpeWgt_P1));
disp('  ')

%% 2. Efficient Frontier with additional constraints

% Create a market structure
mkt = struct();

% Define sectors
mkt.sector.cyclical = ["ConsumerDiscretionary", "Financials", "Materials", "RealEstate", "Industrials"];
mkt.sector.defensive = ["ConsumerStaples", "Utilities", "HealthCare"];
mkt.sector.sensible = ["Energy", "InformationTechnology", "CommunicationServices"];
% Define factors
mkt.factor = ["Momentum", "Value", "Growth", "Quality", "LowVolatility"];

sectors = [mkt.sector.sensible, mkt.sector.cyclical, mkt.sector.defensive];
factor = mkt.factor;

% Create a portfolio object
P2 = Portfolio('AssetList', names);
P2 = setAssetMoments(P2, mean_returns, cov_matrix);
% Define the additional constraints
% 1. Standard constraints, all weights sum to 1, no shorting, and 100% investment in risky assets
% 2. Total exposure to sensible sectors (Energy, Information Technology, Communication Services) 
%    should be greater than 10% and total exposure to cyclical sectors should be lower than 30%
% 3. Consumer Staples and Low Volatility should not be included in the portfolio (weights set to 0)
% 4. Maximum exposure to sectors (Sensible, Cyclical, Defensive) should be less than 80%

% Find the indeces
sensibleIdx = ismember(P2.AssetList, mkt.sector.sensible);
cyclicalIdx = ismember(P2.AssetList, mkt.sector.cyclical);
excludeIdx_CS = ismember(P2.AssetList, "ConsumerStaples");
excludeIdx_LV = ismember(P2.AssetList, "LowVolatility");
sectorIdx = ismember(P2.AssetList, sectors);
%%
Aineq = [
    -sensibleIdx;
    cyclicalIdx;
    excludeIdx_CS;
    excludeIdx_LV;
    sectorIdx;
];

bineq = [-0.1; 0.3; 0; 0; 0.8];
% set the constraints
P2 = setDefaultConstraints(P2);
P2.AInequality = Aineq;
P2.bInequality = bineq;

% Estimate efficient frontier with constraints
pwgt2 = estimateFrontier(P2, 100);
% Estimate portfolio moments
[pf_risk_P2, pf_Retn_P2] = estimatePortMoments(P2, pwgt2);
% Plot efficient frontier with constraints
% plot(pf_risk_P2, pf_Retn_P2)

% Find minimum risk portfolio (Portfolio C - Minimum Variance Portfolio with constraints)
minRisk_P2 = min(pf_risk_P2);
minRiskWgt_P2 = pwgt2(:, pf_risk_P2 == minRisk_P2);
minRiskRet_P2 = mean_returns' * minRiskWgt_P2;
minRiskSR_P2 = (minRiskRet_P2 - risk_free_rate) / minRisk_P2;

% Find maximum Sharpe ratio portfolio (Portfolio D - Maximum Sharpe Ratio Portfolio with constraints)
maxSharpeWgt_P2 = estimateMaxSharpeRatio(P2);
maxSharpeRet_P2 = mean_returns' * maxSharpeWgt_P2;
maxSharpeRisk_P2 = sqrt(maxSharpeWgt_P2' * cov_matrix * maxSharpeWgt_P2);
maxSharpeSR_P2 = (maxSharpeRet_P2 - risk_free_rate) / maxSharpeRisk_P2;

% Display Portfolio C - Minimum Variance Portfolio with constraints
disp('===========================================================================')
disp('Minimum Risk Portfolio with Constraints (Portfolio C)')
disp('===========================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(minRiskWgt_P2)
    fprintf('%-25s %.4f\n', names{i}, minRiskWgt_P2(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', minRiskRet_P2);
fprintf('%-25s %.4f\n', 'Volatility', minRisk_P2);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', minRiskSR_P2);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(minRiskWgt_P2));
disp('  ')

% Display Portfolio D - Maximum Sharpe Ratio Portfolio with constraints
disp('===========================================================================')
disp('Maximum Sharpe Ratio Portfolio with Constraints (Portfolio D)')
disp('===========================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(maxSharpeWgt_P2)
    fprintf('%-25s %.4f\n', names{i}, maxSharpeWgt_P2(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', maxSharpeRet_P2);
fprintf('%-25s %.4f\n', 'Volatility', maxSharpeRisk_P2);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', maxSharpeSR_P2);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(maxSharpeWgt_P2));
disp('  ')

%% 3. Efficient Frontier and Resampling Method
clc

% number of resampling simulations
N = 100;
% # of points along the efficient frontier
num_frontier_points = 100;

% Initialize variables
Ret_P1_sim = zeros(num_frontier_points, N);
Risk_P1_sim = zeros(num_frontier_points, N);
Weights_P1_sim = zeros(num_assets, num_frontier_points, N);

Ret_P2_sim = zeros(num_frontier_points, N);
Risk_P2_sim = zeros(num_frontier_points, N);
Weights_P2_sim = zeros(num_assets, num_frontier_points, N);

for n = 1:N
    
    resampledReturns = mvnrnd(mean_returns, cov_matrix, 252); % 252 mi suggerisce chat per resaplare i rendimenti giornalieri i.e.  252 giorni
    New_mean_returns = mean(resampledReturns)';               % altrimenti mi dice che non ho un resampling completo ma non sono convinto
    NewCov = cov(resampledReturns);

    P1_sim = setAssetMoments(P1, New_mean_returns, NewCov);
    P2_sim = setAssetMoments(P2, New_mean_returns, NewCov);

    w1_sim = estimateFrontier(P1_sim, num_frontier_points);
    w2_sim = estimateFrontier(P2_sim, num_frontier_points);

    [pf_risk_P1_sim, pf_Retn_P1_sim] = estimatePortMoments(P1_sim, w1_sim);
    [pf_risk_P2_sim, pf_Retn_P2_sim] = estimatePortMoments(P2_sim, w2_sim);

    Ret_P1_sim(:,n) = pf_Retn_P1_sim;
    Risk_P1_sim(:, n) = pf_risk_P1_sim;
    Weights_P1_sim(:,:,n) = w1_sim;

    Ret_P2_sim(:,n) = pf_Retn_P2_sim;
    Risk_P2_sim(:, n) = pf_risk_P2_sim;
    Weights_P2_sim(:,:,n) = w2_sim;
end

pwgt_1 = mean(Weights_P1_sim, 3);
pf_risk_1 = mean(Risk_P1_sim, 2);
pf_Retn_1 = mean(Ret_P1_sim, 2);

pwgt_2 = mean(Weights_P2_sim, 3);
pf_risk_2 = mean(Risk_P2_sim, 2);
pf_Retn_2 = mean(Ret_P2_sim, 2);

% Find minimum risk portfolio (Portfolio E - Minimum Variance Portfolio with resampling)
[minRisk_P1_Rsim, idx_minRisk_P1] = min(pf_risk_1);
minRiskWgt_P1_Rsim = pwgt_1(:, idx_minRisk_P1);
minRiskRet_P1_Rsim = pf_Retn_1(idx_minRisk_P1);
minRiskSR_P1_Rsim = (minRiskRet_P1_Rsim - risk_free_rate) / minRisk_P1_Rsim;

% Find minimum risk portfolio (Portfolio F - Minimum Variance Portfolio with resampling and constraints)
[minRisk_P2_Rsim, idx_minRisk_P2] = min(pf_risk_2);
minRiskWgt_P2_Rsim = pwgt_2(:, idx_minRisk_P2);
minRiskRet_P2_Rsim = pf_Retn_2(idx_minRisk_P2);
minRiskSR_P2_Rsim = (minRiskRet_P2_Rsim - risk_free_rate) / minRisk_P2_Rsim;

% Find maximum Sharpe ratio portfolio (Portfolio G - Maximum Sharpe Ratio Portfolio with resampling)
sharpeRatio_P1_Rsim = (pf_Retn_1 - risk_free_rate) ./ pf_risk_1;
[maxSharpeSR_P1_Rsim, idx_maxSharpe_P1] = max(sharpeRatio_P1_Rsim);
maxSharpeWgt_P1_Rsim = pwgt_1(:, idx_maxSharpe_P1);
maxSharpeRet_P1_Rsim = pf_Retn_1(idx_maxSharpe_P1);
maxSharpeRisk_P1_Rsim = pf_risk_1(idx_maxSharpe_P1);

% Find maximum Sharpe ratio portfolio (Portfolio H - Maximum Sharpe Ratio Portfolio with resampling and constraints)
sharpeRatio_P2_Rsim = (pf_Retn_2 - risk_free_rate) ./ pf_risk_2;
[maxSharpeSR_P2_Rsim, idx_maxSharpe_P2] = max(sharpeRatio_P2_Rsim);
maxSharpeWgt_P2_Rsim = pwgt_2(:, idx_maxSharpe_P2);
maxSharpeRet_P2_Rsim = pf_Retn_2(idx_maxSharpe_P2);
maxSharpeRisk_P2_Rsim = pf_risk_2(idx_maxSharpe_P2);

% Display Portfolio E - Minimum Variance Portfolio with resampling
disp('==============================================================================================')
disp('Minimum Risk Portfolio with Resampling (Portfolio E)')
disp('==============================================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(minRiskWgt_P1_Rsim)
    fprintf('%-25s %.4f\n', names{i}, minRiskWgt_P1_Rsim(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', minRiskRet_P1_Rsim);
fprintf('%-25s %.4f\n', 'Volatility', minRisk_P1_Rsim);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', minRiskSR_P1_Rsim);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(minRiskWgt_P1_Rsim));
disp('              ')

% Display Portfolio F - Minimum Variance Portfolio with resampling and constraints
disp('==============================================================================================')
disp('Minimum Risk Portfolio with Resampling and Constraints (Portfolio F)')
disp('==============================================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(minRiskWgt_P2_Rsim)
    fprintf('%-25s %.4f\n', names{i}, minRiskWgt_P2_Rsim(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', minRiskRet_P2_Rsim);
fprintf('%-25s %.4f\n', 'Volatility', minRisk_P2_Rsim);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', minRiskSR_P2_Rsim);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(minRiskWgt_P2_Rsim));
disp('              ')

% Display Portfolio G - Maximum Sharpe Ratio Portfolio with resampling
disp('==============================================================================================')
disp('Maximum Sharpe Ratio Portfolio with Resampling (Portfolio G)')
disp('==============================================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(maxSharpeWgt_P1_Rsim)
    fprintf('%-25s %.4f\n', names{i}, maxSharpeWgt_P1_Rsim(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', maxSharpeRet_P1_Rsim);
fprintf('%-25s %.4f\n', 'Volatility', maxSharpeRisk_P1_Rsim);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', maxSharpeSR_P1_Rsim);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(maxSharpeWgt_P1_Rsim));
disp('                      ')

% Display Portfolio H - Maximum Sharpe Ratio Portfolio with resampling and constraints
disp('==============================================================================================')
disp('Maximum Sharpe Ratio Portfolio with Resampling and Constraints (Portfolio H)')
disp('==============================================================================================')
disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(maxSharpeWgt_P2_Rsim)
    fprintf('%-25s %.4f\n', names{i}, maxSharpeWgt_P2_Rsim(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', maxSharpeRet_P2_Rsim);
fprintf('%-25s %.4f\n', 'Volatility', maxSharpeRisk_P2_Rsim);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', maxSharpeSR_P2_Rsim);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(maxSharpeWgt_P2_Rsim));
disp('          ')

%% 5.  Compute the Maximum Diversified Portfolio and the Maximum Entropy (in asset volatility) Portfolio
% Compute the Maximum Diversified Portfolio (Portfolio M) and the
% Maximum Entropy (in asset volatility) Portfolio (Portfolio N), under
% the following constraints (to be considered all at once):
% - Standard constraints,
% - The total exposure on cyclicals has to be greater than 20%,
% - Assuming that you have a benchmark portfolio (capitalization
%   weighted portfolio), the sum of the difference (in absolute value)
%   of the weights in the benchmark portfolio and the optimal weights
%   has to be greater than 20%

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
entropy = @(w) sum(w.^2' * diag(cov_matrix)/ sum(w.^2' * diag(cov_matrix)) * log(w.^2' * diag(cov_matrix)/ sum(w.^2' * diag(cov_matrix))));

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


%% Output carino dei portafogli insieme 

% Combine the weights for each portfolio
weights_A = array2table(minRiskWgt_P1, 'RowNames', names, 'VariableNames', {'Portfolio_A'});
weights_B = array2table(maxSharpeWgt_P1, 'RowNames', names, 'VariableNames', {'Portfolio_B'});
weights_C = array2table(minRiskWgt_P2, 'RowNames', names, 'VariableNames', {'Portfolio_C'});
weights_D = array2table(maxSharpeWgt_P2, 'RowNames', names, 'VariableNames', {'Portfolio_D'});
weights_E = array2table(minRiskWgt_P1_Rsim, 'RowNames', names, 'VariableNames', {'Portfolio_E'});
weights_F = array2table(minRiskWgt_P2_Rsim, 'RowNames', names, 'VariableNames', {'Portfolio_F'});
weights_G = array2table(maxSharpeWgt_P1_Rsim, 'RowNames', names, 'VariableNames', {'Portfolio_G'});
weights_H = array2table(maxSharpeWgt_P2_Rsim, 'RowNames', names, 'VariableNames', {'Portfolio_H'});
weights_M = array2table(weights_m, 'RowNames', names, 'VariableNames', {'Portfolio_M'});
weights_N = array2table(weights_n, 'RowNames', names, 'VariableNames', {'Portfolio_N'});

weightsTable = [weights_A, weights_B, weights_C, weights_D, weights_E, weights_F, weights_G, weights_H, weights_M, weights_N];

metricsTable = table(...
    [minRiskRet_P1; minRisk_P1; minRiskSR_P1; sum(minRiskWgt_P1)], ...
    [maxSharpeRet_P1; maxSharpeRisk_P1; maxSharpeSR_P1; sum(maxSharpeWgt_P1)], ...
    [minRiskRet_P2; minRisk_P2; minRiskSR_P2; sum(minRiskWgt_P2)], ...
    [maxSharpeRet_P2; maxSharpeRisk_P2; maxSharpeSR_P2; sum(maxSharpeWgt_P2)], ...
    [minRiskRet_P1_Rsim; minRisk_P1_Rsim; minRiskSR_P1_Rsim; sum(minRiskWgt_P1_Rsim)], ...
    [minRiskRet_P2_Rsim; minRisk_P2_Rsim; minRiskSR_P2_Rsim; sum(minRiskWgt_P2_Rsim)], ...
    [maxSharpeRet_P1_Rsim; maxSharpeRisk_P1_Rsim; maxSharpeSR_P1_Rsim; sum(maxSharpeWgt_P1_Rsim)], ...
    [maxSharpeRet_P2_Rsim; maxSharpeRisk_P2_Rsim; maxSharpeSR_P2_Rsim; sum(maxSharpeWgt_P2_Rsim)], ...
    [portfolio_m_return; portfolio_m_std; portfolio_m_SR; sum(weights_m)], ...
    [portfolio_n_return; portfolio_n_std; portfolio_n_SR; sum(weights_n)], ...
    'VariableNames', {'Portfolio_A', 'Portfolio_B', 'Portfolio_C', 'Portfolio_D', ...
                      'Portfolio_E', 'Portfolio_F', 'Portfolio_G', 'Portfolio_H', ...
                        'Portfolio_M', 'Portfolio_N'}, ...
    'RowNames', {'Expected Return', 'Volatility', 'Sharpe Ratio', 'Sum of Weights'} ...
);

% Display the weights table
disp('==============================================================================================')
disp('                                    Portfolio Weights Table                                   ')
disp('==============================================================================================')
disp(weightsTable)

% Display the metrics table
disp('==============================================================================================')
disp('                                    Portfolio Metrics Table                                   ')
disp('==============================================================================================')
disp(metricsTable)


