% Portfolio Managemet project - Main script
% Authors: 
% Matteo Torba
% Edoardo Pariani
% Nicolò Toia
% Andrea Tarditi
% Giacomo Manfredi
%
% Date: 2024-11-13

clc
% clear all
close all
warning off

rng(42); % For reproducibility

% start timer
tic

flag_plot = 0;
%% PART A 

% Load data from Excel files
prices = readtable('prices.xlsx');
capitalizations = readtable('capitalizations.xlsx');
caps = capitalizations{1,2:end};

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

mean_returns = mean(returns_2023)';
cov_matrix = cov(returns_2023);

% Define number of assets
num_assets = length(names);

% Set the risk-free rate = 0
risk_free_rate = 0; 

if flag_plot == 1
    plotData(returns_2023,prices_2023, names);
end

%% market structure

% Create a market structure
mkt = struct();

% Define sectors
mkt.sector.cyclical = ["ConsumerDiscretionary", "Financials", "Materials", "RealEstate", "Industrials"];
mkt.sector.defensive = ["ConsumerStaples", "Utilities", "HealthCare"];
mkt.sector.sensible = ["Energy", "InformationTechnology", "CommunicationServices"];
% Define factors
mkt.factor = ["Momentum", "Value", "Growth", "Quality", "LowVolatility"];

Ptf_AMkt = Portfolio('AssetList', names);
sectors = [mkt.sector.sensible, mkt.sector.cyclical, mkt.sector.defensive];

sensibleIdx = ismember(Ptf_AMkt.AssetList, mkt.sector.sensible);
cyclicalIdx = ismember(Ptf_AMkt.AssetList, mkt.sector.cyclical);
excludeIdx_CS = ismember(Ptf_AMkt.AssetList, "ConsumerStaples");
excludeIdx_LV = ismember(Ptf_AMkt.AssetList, "LowVolatility");
sectorIdx = ismember(Ptf_AMkt.AssetList, sectors);

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

% estimate efficient frontier, via object method and 100 Ptfoints
pwgt1 = estimateFrontier(P1, 100);
% estimate portfolio moments
[pf_risk_Ptf_1, pf_Retn_Ptf_1] = estimatePortMoments(P1, pwgt1);

% Portfolio A: Minimum Variance Portfolio
Portfolio_A = minRiskPortfolio(P1, pwgt1, pf_risk_Ptf_1, 'Minimum Risk Portfolio (A)');

% Portfolio B: Maximum Sharpe Ratio Portfolio
Portfolio_B = maxSharpPortfolio(P1, 'Max sharpe ratio Portfolio (B)');

% plot efficient frontier
% plot(pf_risk_Ptf_1, pf_Retn_Ptf_1)

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

% estimate efficient frontier, via object method and 100 Ptfoints
pwgt2 = estimateFrontier(P2, 100);
% estimate portfolio moments
[pf_risk_Ptf_2, pf_Retn_Ptf_2] = estimatePortMoments(P2, pwgt2);

% Portfolio C: Minimum Variance Portfolio with specific constraints
Portfolio_C = minRiskPortfolio(P2, pwgt2, pf_risk_Ptf_2, 'Minimum risk Portfolio with constrints (C)');

% Portfolio D: Maximum Sharpe Ratio Portfolio with specific constraints
Portfolio_D = maxSharpPortfolio(P2, 'Max sharpe ratio Portfolio with constrints (D)');

% plot efficient frontier
% plot(pf_risk_Ptf_2, pf_Retn_Ptf_2)

%% 3. Efficient Frontier and Resampling Method

% Compute the frontiers in step 1 and 2 using the resampling method in
%  order to obtain 2 robust frontier. For each frontier save the Minimum
%  Variance Portfolios, named Portfolios E and F, and the Maximum
%  Sharpe Ratio Portfolios, named Portfolios G and H, of the frontiers.

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

% [ptf_I, ptf_L] = BlackLitterman_2(P3, returns_lnr, caps, lambda, [view1, view2], 0);

% Portfolio I: Minimum Variance Portfolio with Black-Litterman
flag = 0;
[Portfolio_I, ~] = BlackLitterman_2(P3, returns_lnr, caps, lambda, [view1, view2], 'Minimum Variance Portfolio with Black Litterman (I)', flag);

% Portfolio L: Maximum Sharpe Ratio Portfolio with Black-Litterman
flag = 1;
[Portfolio_L, P3] = BlackLitterman_2(P3, returns_lnr, caps, lambda, [view1, view2], 'Maximum Sharpe Ratio Portfolio with Black Litterman (L)', flag);

%% 5.  Compute the Maximum Diversified Portfolio and the Maximum Entropy (in asset volatility) Portfolio

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

% Weights of the benchmark portfolio (capitalization weighted portfolio)
%  cap_wghtd_ptf = capitalizations{1, names} / sum(capitalizations{1, names});

cap_wghtd_ptf = caps / sum(caps);

 const.nonlinconstr =  @(weights) customAbsDiffConstraint(weights, cap_wghtd_ptf); % Non linear constraint function

% Portfolio M: Maximum Diversified Portfolio
[weights_m, portfolio_m_return, portfolio_m_std, portfolio_m_SR] = ...
          Max_Diversified_Portfolio(mean_returns, cov_matrix, capitalizations, names, risk_free_rate, mkt, P2,num_assets,const);

% Portfolio N: Maximum Entropy Portfolio
[weights_n, portfolio_n_return, portfolio_n_std, portfolio_n_SR] = ...
    Max_Entropy_Portfolio(mean_returns, cov_matrix, capitalizations, names, risk_free_rate, mkt, P2,num_assets,const);

%% 6. 

% Compute the portfolio (Portfolio P), using the Principal Component
%  Analysis, that maximizes its expected return under the following con
% straints (to be considered all at once):
%  • Standard constraints,
%  • Thevolatility of the portfolio has to be equal or less than a target
%  volatility of σtgt = 0.1
%  You have to use the minimum number of factors that explains more
%  than the 90% of the cumulative variance.

%% 7.

% Compute the Portfolio that maximizes, under standard constraints,
%  the VaR-modified Sharpe Ratio (i.e. the risk in the formula of Sharpe
%  Ratio is the VaR), named Portfolio Q, using the Variance-Covariance
%  method.


%% Output carino dei portafogli insieme 

% % Combine the weights for each portfolio
% weights_A = array2table(minRiskWgt_P1, 'RowNames', names, 'VariableNames', {'Portfolio A'});
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
% weights_B = array2table(maxSharpeWgt_P1, 'RowNames', names, 'VariableNames', {'Portfolio B'});
% weights_C = array2table(minRiskWgt_P2, 'RowNames', names, 'VariableNames', {'Portfolio C'});
% weights_D = array2table(maxSharpeWgt_P2, 'RowNames', names, 'VariableNames', {'Portfolio D'});
% weights_E = array2table(minRiskWgt_P1_Rsim, 'RowNames', names, 'VariableNames', {'Portfolio E'});
% weights_F = array2table(minRiskWgt_P2_Rsim, 'RowNames', names, 'VariableNames', {'Portfolio F'});
% weights_G = array2table(maxSharpeWgt_P1_Rsim, 'RowNames', names, 'VariableNames', {'Portfolio G'});
% weights_H = array2table(maxSharpeWgt_P2_Rsim, 'RowNames', names, 'VariableNames', {'Portfolio H'});
% weights_I = array2table(ptf_I.w, 'RowNames', names, 'VariableNames', {'Portfolio I'});
% weights_L = array2table(ptf_L.w, 'RowNames', names, 'VariableNames', {'Portfolio L'});


% weights_M = array2table(weights_m, 'RowNames', names, 'VariableNames', {'Portfolio M'});
% weights_N = array2table(weights_n, 'RowNames', names, 'VariableNames', {'Portfolio N'});

% weightsTable = [weights_A, weights_B, weights_C, weights_D,...
%                 weights_E, weights_F, weights_G, weights_H,...
%                 weights_I, weights_L, weights_M, weights_N];
weightsTable = [weights_A, weights_B, weights_C, weights_D,...
                weights_E, weights_F, weights_G, weights_H,...
                weights_I, weights_L];

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
    'VariableNames', {'Portfolio A', 'Portfolio B', 'Portfolio C', 'Portfolio D', ...
                      'Portfolio E', 'Portfolio F', 'Portfolio G', 'Portfolio H', ...
                      'Portfolio I', 'Portfolio L'}, ...
    'RowNames', {'Expected Return', 'Volatility', 'Sharpe Ratio', 'Sum of Weights'} ...
);

% metricsTable = table(...
%     [minRiskRet_P1; minRisk_P1; minRiskSR_P1; sum(minRiskWgt_P1)], ...
%     [maxSharpeRet_P1; maxSharpeRisk_P1; maxSharpeSR_P1; sum(maxSharpeWgt_P1)], ...
%     [minRiskRet_P2; minRisk_P2; minRiskSR_P2; sum(minRiskWgt_P2)], ...
%     [maxSharpeRet_P2; maxSharpeRisk_P2; maxSharpeSR_P2; sum(maxSharpeWgt_P2)], ...
%     [minRiskRet_P1_Rsim; minRisk_P1_Rsim; minRiskSR_P1_Rsim; sum(minRiskWgt_P1_Rsim)], ...
%     [minRiskRet_P2_Rsim; minRisk_P2_Rsim; minRiskSR_P2_Rsim; sum(minRiskWgt_P2_Rsim)], ...
%     [maxSharpeRet_P1_Rsim; maxSharpeRisk_P1_Rsim; maxSharpeSR_P1_Rsim; sum(maxSharpeWgt_P1_Rsim)], ...
%     [maxSharpeRet_P2_Rsim; maxSharpeRisk_P2_Rsim; maxSharpeSR_P2_Rsim; sum(maxSharpeWgt_P2_Rsim)], ...
%     [ptf_I.ret; ptf_I.std; ptf_I.sr; sum(ptf_I.w)], ...
%     [ptf_L.ret; ptf_L.std; ptf_L.sr; sum(ptf_L.w)], ...
%     [portfolio_m_return; portfolio_m_std; portfolio_m_SR; sum(weights_m)], ...
%     [portfolio_n_return; portfolio_n_std; portfolio_n_SR; sum(weights_n)], ...
%     'VariableNames', {'Portfolio A', 'Portfolio B', 'Portfolio C', 'Portfolio D', ...
%                       'Portfolio E', 'Portfolio F', 'Portfolio G', 'Portfolio H', ...
%                       'Portfolio I', 'Portfolio L', 'Portfolio M', 'Portfolio N'}, ...
%     'RowNames', {'Expected Return', 'Volatility', 'Sharpe Ratio', 'Sum of Weights'} ...
% );

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

disp_pie_weights(weightsTable)

% Display the metrics table 2023
disp('==============================================================================================')
disp('                                    Portfolio Performance Metrics 2023                                 ')
disp('==============================================================================================')
[eq, performancesMetrics] = getEquityandMetrices(weightsTable, prices_2023, "2023");
disp(performancesMetrics)

% % Display the metrics table 2024
% disp('==============================================================================================')
% disp('                                    Portfolio Performance Metrics 2024                                 ')
% disp('==============================================================================================')
% [eq, performancesMetrics] = getEquityandMetrices(weightsTable, prices_2024, "2024");
% disp(performancesMetrics)


% Display the time taken
disp('==============================================================================================')
disp(['Time taken: ', num2str(toc), ' seconds'])
disp('==============================================================================================')