% Portfolio Managemet project - Main script
% Authors: 
% Matteo Torba
% Edoardo Pariani
% Nicolò Toia
% Andrea Tarditi
% Giacomo Manfredi
%
% Date: 2024-11-13

clc, clear, close all;
% warning off

rng(42); % For reproducibility

% start timer
tic

% Set the flag for plotting the data analysis, 1: plot, 0: do not plot
flag_plot = 1;

% Add paths
addpath('Data')
addpath('Portfolio_Methods')
addpath('Print_and_Display')
addpath('Analysis_tools')

%% PART A 

% Load data from Excel files
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames'); % Suppress the warning on variables' names
prices = readtable('prices.xlsx');
capitalizations = readtable('capitalizations.xlsx');
caps = capitalizations{1,2:end};
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames'); % Re-enable the warning on variables' names

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
dates_2023 = dates(dates >= start_date & dates <= end_date, :);
% Filter prices for 2024
prices_2024 = prices_data(dates > end_date, :);
dates_2024 = dates(dates > end_date, :);

% Calculate daily log-returns for each index in 2023
returns_2023 = diff(log(prices_2023));

mean_returns = mean(returns_2023)';
cov_matrix = cov(returns_2023);

% Define number of assets
num_assets = length(names);

% Set the risk-free rate = 0
risk_free_rate = 0; 

if flag_plot == 1
    plotData(returns_2023,prices_2023, names);
end

%% Market structure 2023

% Create a market structure
mkt = struct();

% Define sectors
mkt.sector.cyclical = ["ConsumerDiscretionary", "Financials", "Materials", "RealEstate", "Industrials"];
mkt.sector.defensive = ["ConsumerStaples", "Utilities", "HealthCare"];
mkt.sector.sensible = ["Energy", "InformationTechnology", "CommunicationServices"];
% Define factors
mkt.factor = ["Momentum", "Value", "Growth", "Quality", "LowVolatility"];

% Define the desired order of columns to cluster the data in sectors and factors
desiredOrder = [mkt.sector.cyclical, mkt.sector.defensive, mkt.sector.sensible, mkt.factor];

% Define a Portfolio object for the All Market
Ptf_AMkt = Portfolio('AssetList', names);
sectors = [mkt.sector.sensible, mkt.sector.cyclical, mkt.sector.defensive];

% Define indices for sectors and factors used in the constraints
sensibleIdx = ismember(Ptf_AMkt.AssetList, mkt.sector.sensible);
cyclicalIdx = ismember(Ptf_AMkt.AssetList, mkt.sector.cyclical);
excludeIdx_CS = ismember(Ptf_AMkt.AssetList, "ConsumerStaples");
excludeIdx_LV = ismember(Ptf_AMkt.AssetList, "LowVolatility");
sectorIdx = ismember(Ptf_AMkt.AssetList, sectors);

%% EW Portfolio
% Create an Equally Weighted Portfolio
Portfolio_EW = struct();
Portfolio_EW.Weights = ones(num_assets,1)/num_assets;
Portfolio_EW.Return = mean_returns' * Portfolio_EW.Weights;
Portfolio_EW.Volatility = ...
    sqrt(Portfolio_EW.Weights' * cov_matrix * Portfolio_EW.Weights);
Portfolio_EW.Sharpe_Ratio = (Portfolio_EW.Return - risk_free_rate)...
    / Portfolio_EW.Volatility;
%% Caps Portfolio
% Create a Portfolio with weights proportional to the capitalizations
Portfolio_Caps = struct();
Portfolio_Caps.Weights = caps' ./ sum(caps);
Portfolio_Caps.Return = mean_returns' * Portfolio_Caps.Weights;
Portfolio_Caps.Volatility = ...
    sqrt(Portfolio_Caps.Weights' * cov_matrix * Portfolio_Caps.Weights);
Portfolio_Caps.Sharpe_Ratio = (Portfolio_Caps.Return - risk_free_rate)...
    / Portfolio_Caps.Volatility;

%% 1. Efficient Frontier

%   Compute the efficient frontier under the standard constraints.
%   Compute the Minimum Variance Portfolio, named Portfolio A, and the Maximum Sharpe 
%   Ratio Portfolio, named Portfolio B, of the frontier.

% Create a struct for the constraints
const_std = struct();
% set constraints
const_std.Aineq = [];
const_std.bineq = [];

% Create a portfolio object
P1 = Portfolio('AssetList', names);
P1 = setDefaultConstraints(P1); % all weights sum to 1, no shorting, and 100% investment in risky assets
P1 = setAssetMoments(P1, mean_returns, cov_matrix); % set mean returns and covariance matrix
P1.Name = 'Standard Portfolio';
P1.RiskFreeRate = risk_free_rate;
% set constraints of the portfolio
P1.AInequality = const_std.Aineq;
P1.bInequality = const_std.bineq;

% estimate efficient frontier, via object method and 100 Ptf points
pwgt1 = estimateFrontier(P1, 100);
% estimate portfolio moments
[pf_risk_Ptf_1, pf_Retn_Ptf_1] = estimatePortMoments(P1, pwgt1);

% Portfolio A: Minimum Variance Portfolio
Portfolio_A = minRiskPortfolio(P1, pwgt1, pf_risk_Ptf_1, 'Minimum Risk Portfolio (A)');

% Portfolio B: Maximum Sharpe Ratio Portfolio
Portfolio_B = maxSharpPortfolio(P1, 'Max sharpe ratio Portfolio (B)');

% Plot the efficient frontier
if flag_plot
    % Plot the efficient frontier
    portfolios_vector = [Portfolio_A; Portfolio_B];
    titleEntry = 'Efficient Frontier Under Standard Constraints';
    legendEntries = {'Efficient Frontier', 'Minimum Variance Portfolio (A)', 'Max sharpe ratio Portfolio (B)'};
    plot_frontier(pf_risk_Ptf_1, pf_Retn_Ptf_1, portfolios_vector, titleEntry, legendEntries)
    
    % Plot the weights for each asset
    plotFrontierWeights(pf_risk_Ptf_1, pwgt1, names, desiredOrder)
end

%% 2. Efficient Frontier with additional constraints

%  Compute the efficient frontier under the following constraints (to be
%  considered all at once):
%  • Standard constraints,
%  • The total exposure to sensible sectors has to be greater than 10%
%  and the total exposure on cyclical sectors has to be lower than
%  30%,
%  • The weights of the sector Consumer Staples and the factor
%  Low Volatility have to be set equal to 0,
%  • The maximum exposure on sectors has to be lower than 80%.
%  Compute the Minimum Variance Portfolio, named Portfolio C, and the
%  Maximum Sharpe Ratio Portfolio, named Portfolio D, of the frontier.

% Create a struct for the constraints
const_spec = struct();
% set constraints
const_spec.Aineq = [
    -sensibleIdx;
    cyclicalIdx;
    excludeIdx_CS;
    excludeIdx_LV;
    sectorIdx;
    ];
const_spec.bineq = [-0.1; 0.3; 0; 0; 0.8];

% Create a portfolio object
P2 = Portfolio('AssetList', names);
P2 = setDefaultConstraints(P2); % all weights sum to 1, no shorting, and 100% investment in risky assets
P2 = setAssetMoments(P2, mean_returns, cov_matrix); % set mean returns and covariance matrix
P2.Name = 'Portfolio with constraints';
P2.RiskFreeRate = risk_free_rate;
% set constraints of the portfolio
P2.AInequality = const_spec.Aineq;
P2.bInequality = const_spec.bineq;

% estimate efficient frontier, via object method and 100 Ptf points
pwgt2 = estimateFrontier(P2, 100);
% estimate portfolio moments
[pf_risk_Ptf_2, pf_Retn_Ptf_2] = estimatePortMoments(P2, pwgt2);

% Portfolio C: Minimum Variance Portfolio with specific constraints
Portfolio_C = minRiskPortfolio(P2, pwgt2, pf_risk_Ptf_2, 'Minimum risk Portfolio with constraints (C)');

% Portfolio D: Maximum Sharpe Ratio Portfolio with specific constraints
Portfolio_D = maxSharpPortfolio(P2, 'Max sharpe ratio Portfolio with constraints (D)');

% Plot the efficient frontier
if flag_plot
    % Plot the efficient frontier
    risk_frontiers = [pf_risk_Ptf_1 pf_risk_Ptf_2];
    return_frontiers = [pf_Retn_Ptf_1 pf_Retn_Ptf_2];
    portfolios_vector = [Portfolio_A; Portfolio_B; Portfolio_C; Portfolio_D];
    titleEntry = 'Efficient Frontier Under Additonal Constraints';
    legendEntries = {'Standard Efficient Frontier', 'Efficient Frontier with constraints', 'Minimum Variance Portfolio (A)', 'Max sharpe ratio Portfolio (B)',...
        'Minimum risk Portfolio with constraints (C)', 'Max sharpe ratio Portfolio with constraints (D)'};
    plot_frontier(risk_frontiers, return_frontiers, portfolios_vector, titleEntry, legendEntries)

    % Plot the weights for each asset
    plotFrontierWeights(pf_risk_Ptf_2, pwgt2, names, desiredOrder)
end

%% 3. Efficient Frontier and Resampling Method

% Compute the frontiers in step 1 and 2 using the resampling method in
% order to obtain 2 robust frontier. For each frontier save the Minimum
% Variance Portfolios, named Portfolios E and F, and the Maximum
% Sharpe Ratio Portfolios, named Portfolios G and H, of the frontiers.

% Set the flag for Minimum Variance Portfolio or Maximum Sharpe Ratio Portfolio
% 0: Minimum Variance Portfolio
% 1: Maximum Sharpe Ratio Portfolio

% Portfolio E: Minimum Variance Portfolio with resampling
flag = 0;
Portfolio_E = resampling_method(P1, 'Minimum risk Portfolio with resampling (E)', flag);

% Portfolio F: Minimum Variance Portfolio with resampling and constraints
Portfolio_F = resampling_method(P2, 'Minimum risk Portfolio with resampling and constraints (F)', flag);

% Portfolio G: Maximum Sharpe Ratio Portfolio with resampling
flag = 1;
Portfolio_G = resampling_method(P1, 'Max sharpe ratio Portfolio with resampling (G)', flag);

% Portfolio H: Maximum Sharpe Ratio Portfolio with resampling and constraints
Portfolio_H = resampling_method(P2, 'Max sharpe ratio Portfolio with resampling and constraints (H)', flag);

% Plot the efficient frontier
if flag_plot
    risk_frontiers = [pf_risk_Ptf_1 pf_risk_Ptf_2];
    return_frontiers = [pf_Retn_Ptf_1 pf_Retn_Ptf_2];
    portfolios_vector = [Portfolio_A; Portfolio_B; Portfolio_C; Portfolio_D; Portfolio_E; Portfolio_F; Portfolio_G; Portfolio_H];
    titleEntry = 'Robust Frontier Portoflios against Original Data Frontiers';
    legendEntries = {'Standard Efficient Frontier', 'Efficient Frontier with constraints', 'Minimum Variance Portfolio (A)', 'Max sharpe ratio Portfolio (B)',...
        'Minimum risk Portfolio with constraints (C)', 'Max sharpe ratio Portfolio with constraints (D)', 'Minimum risk Portfolio with resampling (E)',...
        'Minimum risk Portfolio with resampling and constraints (F)', 'Max sharpe ratio Portfolio with resampling (G)', 'Max sharpe ratio Portfolio with resampling and constraints (H)'};
    plot_frontier(risk_frontiers, return_frontiers, portfolios_vector, titleEntry, legendEntries)
end

%% 4. Black-Litterman Model

% Compute the portfolio frontier, under standard constraints, using the
% Black-Litterman model with the following views (to be considered all
% at once):
% • View on Technology vs. Financials: Given the growing im
% portance of technology, you think that the Technology sector will
% outperform the Financial sector of the 2% (annual).
% • View on the Momentum vs. Low Volatility factor: You
% might assume that Momentum will outperform Low Volatility in
% a bull market scenario (annual overperformance of 1%)

% Set the flag for Minimum Variance Portfolio or Maximum Sharpe Ratio Portfolio
% 0: Minimum Variance Portfolio
% 1: Maximum Sharpe Ratio Portfolio

% View 1: Information Technology - Financials = 2%
% View 2: Momentum - Volatility = 1%
% Number of views
v = 2;
% Describe views with a struct
% View 1
view1 = struct();
view1.overperformer = "InformationTechnology";
view1.underperformer = "Financials";
view1.delta = 2/100;
% View 2
view2 = struct();
view2.overperformer = "Momentum";
view2.underperformer = "LowVolatility";
view2.delta = 1/100;

% Market data
returns_lnr = tick2ret(prices_2023); % Linear returns
% Calculate mean returns and covariance matrix
mean_returns_lnr = mean(returns_lnr)';
cov_matrix_lnr = cov(returns_lnr);
% Set the risk aversion
lambda = 1.2;

% Build a portfolio object
P3 = Portfolio('AssetList', names);
P3 = setDefaultConstraints(P3); % all weights sum to 1, no shorting, and 100% investment in risky assets
P3 = setAssetMoments(P3, mean_returns_lnr, cov_matrix_lnr); % set mean returns and covariance matrix
P3.Name = 'Black-Litterman Portfolio';
P3.RiskFreeRate = risk_free_rate;
P3.NumAssets = num_assets;

% Portfolio I: Minimum Variance Portfolio with Black-Litterman
flag = 0;
[Portfolio_I, ~] = BlackLitterman(P3, returns_lnr, caps, lambda, [view1, view2], 'Minimum Variance Portfolio with Black Litterman (I)', flag);

% Portfolio L: Maximum Sharpe Ratio Portfolio with Black-Litterman
flag = 1;
[Portfolio_L, P3] = BlackLitterman(P3, returns_lnr, caps, lambda, [view1, view2], 'Maximum Sharpe Ratio Portfolio with Black Litterman (L)', flag);

%% 5.  Compute Maximum Diversified Portfolio and the Maximum Entropy (in asset volatility) Portfolio

% Compute the Maximum Diversified Portfolio (Portfolio M) and the
% Maximum Entropy (in asset volatility) Portfolio (Portfolio N), under
% the following constraints (to be considered all at once):
% • Standard constraints,
% • The total exposure on cyclicals has to be greater than 20%,
% • Assuming that you have a benchmark portfolio (capitalization
% weighted portfolio), the sum of the difference (in absolute value)
% of the weights in the benchmark portfolio and the optimal weights
% has to be greater than 20%

% Set default constraints for the portfolio: weights sum to 1
% and there are no short positions (non-negative weights)
const.lb = zeros(num_assets, 1);  % Lower bounds on weights
const.ub = ones(num_assets, 1);   % Upper bounds on weights
const.Aeq = ones(1, num_assets);  % Equality constraint for sum of weights
const.beq = 1;

% Additional constraints for sectors:
% Cyclical sectors total weight >= 20%
const.A = -cyclicalIdx;
const.b = -0.2;

% set the non-linear constraint using the weights of the benchmark portfolio (capitalization weighted portfolio)
const.nonlinconstr =  @(weights) customAbsDiffConstraint(weights, Portfolio_Caps.Weights); % Non linear constraint function

% Build a portfolio object
P4 = Portfolio('AssetList', names);
P4 = setDefaultConstraints(P4); % all weights sum to 1, no shorting, and 100% investment in risky assets
P4 = setAssetMoments(P4, mean_returns, cov_matrix); % set mean returns and covariance matrix
P4.Name = 'Portfolio with constraints';
P4.RiskFreeRate = risk_free_rate;
P4.NumAssets = num_assets;

% Portfolio M: Maximum Diversified Portfolio
Portfolio_M = Max_Diversified_Portfolio(P4, const, 'Max diversified Portfolio (M)');

% Portfolio N: Maximum Entropy Portfolio
Portfolio_N = Max_Entropy_Portfolio(P4, const, 'Max Entropy Portfolio (N)');

% Plot the efficient frontier
% if flag_plot
%     risk_frontiers = [pf_risk_Ptf_1 pf_risk_Ptf_2];
%     return_frontiers = [pf_Retn_Ptf_1 pf_Retn_Ptf_2];
%     portfolios_vector = [Portfolio_A; Portfolio_B; Portfolio_C; Portfolio_D; Portfolio_M; Portfolio_N];
%     titleEntry = 'Efficient Frontier Under Additonal Constraints';
%     legendEntries = {'Standard Efficient Frontier', 'Efficient Frontier with constraints', 'Minimum Variance Portfolio (A)', 'Max sharpe ratio Portfolio (B)',...
%         'Minimum risk Portfolio with constraints (C)', 'Max sharpe ratio Portfolio with constraints (D)', 'Max diversified Portfolio (M)', 'Max Entropy Portfolio (N)'};
%     plot_frontier(risk_frontiers, return_frontiers, portfolios_vector, titleEntry, legendEntries)
% end

%% 6. Portfolio with PCA

% Compute the portfolio (Portfolio P), using the Principal Component
%  Analysis, that maximizes its expected return under the following con
% straints (to be considered all at once):
%  • Standard constraints,
%  • Thevolatility of the portfolio has to be equal or less than a target
%  volatility of σtgt = 0.7
%  You have to use the minimum number of factors that explains more
%  than the 90% of the cumulative variance.

% create a portfolio object
P5 = Portfolio('AssetList', names);
P5 = setDefaultConstraints(P5); % all weights sum to 1, no shorting, and 100% investment in risky assets
P5.Name = 'Portfolio with PCA';
P5.RiskFreeRate = risk_free_rate;
P5.NumAssets = num_assets;

Portfolio_P = PCA_Ptf(P5, returns_2023, 'Portfolio with PCA (P)');

%% 7. Portfolio with VaR modified Sharpe Ratio

% Compute the Portfolio that maximizes, under standard constraints,
%  the VaR-modified Sharpe Ratio (i.e. the risk in the formula of Sharpe
%  Ratio is the VaR), named Portfolio Q, using the Variance-Covariance
%  method.

confidence_level = 0.95; % Confidence level for VaR

% build a portfolio object
P6 = Portfolio('AssetList', names);
P6 = setDefaultConstraints(P6); % all weights sum to 1, no shorting, and 100% investment in risky assets
P6 = setAssetMoments(P6, mean_returns, cov_matrix); % set mean returns and covariance matrix
P6.Name = 'Portfolio with VaR-modified Sharpe Ratio';
P6.RiskFreeRate = risk_free_rate;
P6.NumAssets = num_assets;

Portfolio_Q = VaR_modified_SR(P6, confidence_level, caps, 'VaR modified sharpe ratio Portfolio (Portfolio Q)');

%% Display the results of the portfolios: weights and performance metrics 

% Combine the weights for each portfolio
weights_A = array2table(Portfolio_A.Weights, 'RowNames', names, 'VariableNames', {'Portfolio A'});
weights_B = array2table(Portfolio_B.Weights, 'RowNames', names, 'VariableNames', {'Portfolio B'});
weights_C = array2table(Portfolio_C.Weights, 'RowNames', names, 'VariableNames', {'Portfolio C'});
weights_D = array2table(Portfolio_D.Weights, 'RowNames', names, 'VariableNames', {'Portfolio D'});
weights_E = array2table(Portfolio_E.Weights, 'RowNames', names, 'VariableNames', {'Portfolio E'});
weights_F = array2table(Portfolio_F.Weights, 'RowNames', names, 'VariableNames', {'Portfolio F'});
weights_G = array2table(Portfolio_G.Weights, 'RowNames', names, 'VariableNames', {'Portfolio G'});
weights_H = array2table(Portfolio_H.Weights, 'RowNames', names, 'VariableNames', {'Portfolio H'});
weights_I = array2table(Portfolio_I.Weights, 'RowNames', names, 'VariableNames', {'Portfolio I'});
weights_L = array2table(Portfolio_L.Weights, 'RowNames', names, 'VariableNames', {'Portfolio L'});
weights_M = array2table(Portfolio_M.Weights, 'RowNames', names, 'VariableNames', {'Portfolio M'});
weights_N = array2table(Portfolio_N.Weights, 'RowNames', names, 'VariableNames', {'Portfolio N'});
weights_P = array2table(Portfolio_P.Weights, 'RowNames', names, 'VariableNames', {'Portfolio P'});
weights_Q = array2table(Portfolio_Q.Weights, 'RowNames', names, 'VariableNames', {'Portfolio Q'});
weights_EW = array2table(Portfolio_EW.Weights, 'RowNames', names, 'VariableNames', {'Portfolio EW'});
weights_Caps = array2table(Portfolio_Caps.Weights, 'RowNames', names, 'VariableNames', {'Portfolio Caps'});

weightsTable = [weights_A, weights_B, weights_C, weights_D,...
                weights_E, weights_F, weights_G, weights_H,...
                weights_I, weights_L, weights_M, weights_N,...
                weights_P, weights_Q, weights_EW, weights_Caps];

metricsTable = table(...
    [Portfolio_A.Return; Portfolio_A.Volatility; Portfolio_A.Sharpe_Ratio; sum(Portfolio_A.Weights)], ...
    [Portfolio_B.Return; Portfolio_B.Volatility; Portfolio_B.Sharpe_Ratio; sum(Portfolio_B.Weights)], ...
    [Portfolio_C.Return; Portfolio_C.Volatility; Portfolio_C.Sharpe_Ratio; sum(Portfolio_C.Weights)], ...
    [Portfolio_D.Return; Portfolio_D.Volatility; Portfolio_D.Sharpe_Ratio; sum(Portfolio_D.Weights)], ...
    [Portfolio_E.Return; Portfolio_E.Volatility; Portfolio_E.Sharpe_Ratio; sum(Portfolio_E.Weights)], ...
    [Portfolio_F.Return; Portfolio_F.Volatility; Portfolio_F.Sharpe_Ratio; sum(Portfolio_F.Weights)], ...
    [Portfolio_G.Return; Portfolio_G.Volatility; Portfolio_G.Sharpe_Ratio; sum(Portfolio_G.Weights)], ...
    [Portfolio_H.Return; Portfolio_H.Volatility; Portfolio_H.Sharpe_Ratio; sum(Portfolio_H.Weights)], ...
    [Portfolio_I.Return; Portfolio_I.Volatility; Portfolio_I.Sharpe_Ratio; sum(Portfolio_I.Weights)], ...
    [Portfolio_L.Return; Portfolio_L.Volatility; Portfolio_L.Sharpe_Ratio; sum(Portfolio_L.Weights)], ...
    [Portfolio_M.Return; Portfolio_M.Volatility; Portfolio_M.Sharpe_Ratio; sum(Portfolio_M.Weights)], ...
    [Portfolio_N.Return; Portfolio_N.Volatility; Portfolio_N.Sharpe_Ratio; sum(Portfolio_N.Weights)], ...
    [Portfolio_P.Return; Portfolio_P.Volatility; Portfolio_P.Sharpe_Ratio; sum(Portfolio_P.Weights)], ...
    [Portfolio_Q.Return; Portfolio_Q.Volatility; Portfolio_Q.Sharpe_Ratio; sum(Portfolio_Q.Weights)], ...
    [Portfolio_EW.Return; Portfolio_EW.Volatility; Portfolio_EW.Sharpe_Ratio; sum(Portfolio_EW.Weights)], ...
    [Portfolio_Caps.Return; Portfolio_Caps.Volatility; Portfolio_Caps.Sharpe_Ratio; sum(Portfolio_Caps.Weights)], ...
    'VariableNames', {'Portfolio A', 'Portfolio B', 'Portfolio C', 'Portfolio D', ...
                      'Portfolio E', 'Portfolio F', 'Portfolio G', 'Portfolio H', ...
                      'Portfolio I', 'Portfolio L', 'Portfolio M', 'Portfolio N', ...
                      'Portfolio P', 'Portfolio Q', 'Portfolio EW', 'Portfolio Caps'}, ...
    'RowNames', {'Expected Return', 'Volatility', 'Sharpe Ratio', 'Sum of Weights'} ...
);

% Define the threshold for the weights, e.g. if the weight is smaller than the threshold, it is considered 0
threshold = 1e-4;

% Set the weights smaller than the threshold to 0
weightsTable{:, :} = arrayfun(@(x) x * (abs(x) >= threshold), weightsTable{:, :});

% Display the weights table
disp('==============================================================================================')
disp('                                    Portfolio Weights Table                                   ')
disp('==============================================================================================')
disp(weightsTable)
disp_weights(weightsTable, desiredOrder, mkt)

% Display the metrics table
disp('==============================================================================================')
disp('                                    Portfolio Metrics Table                                   ')
disp('==============================================================================================')
disp(metricsTable)

% Display the metrics table 2023
disp('==============================================================================================')
disp('                                    Portfolio Performance Metrics 2023                                 ')
disp('==============================================================================================')
[~, performancesMetrics2023] = getEquityandMetrices(weightsTable, prices_2023, "2023", dates_2023(2:end));
disp(performancesMetrics2023)
plotMetrics(performancesMetrics2023)

% Display the relRC table 2023
disp('==============================================================================================')
disp('                                    Relative Risk Contributions 2023                                 ')
disp('==============================================================================================')
tableRC_2023 = getRiskContributionsTables(weightsTable, prices_2023);
disp(tableRC_2023)

% Display the metrics table 2024
disp('==============================================================================================')
disp('                                    Portfolio Performance Metrics 2024                                 ')
disp('==============================================================================================')
[eq, performancesMetrics2024] = getEquityandMetrices(weightsTable, prices_2024, "2024", dates_2024(2:end));
disp(performancesMetrics2024)
plotMetrics(performancesMetrics2024)

% Display the relRC table 2024
disp('==============================================================================================')
disp('                                    Relative Risk Contributions 2024                                 ')
disp('==============================================================================================')
tableRC_2024 = getRiskContributionsTables(weightsTable, prices_2024);
disp(tableRC_2024)

% Display the time taken
disp('==============================================================================================')
disp(['Time taken: ', num2str(toc), ' seconds'])
disp('==============================================================================================')