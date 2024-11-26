% Portfolio Managemet project - Main script
% Authors: 
% Matteo Torba
% Edoardo Pariani
% Nicolo Toia
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
risk_free_rate = 0; % Convert to daily (365 is correct)

if flag_plot == 1
    plotData(returns_2023,prices_2023, names);
end

const = struct();

%% 1. Efficient Frontier

%   Compute the efficient frontier under the standard constraints.
%   Compute the Minimum Variance Portfolio, named Portfolio A, and the Maximum Sharpe 
%   Ratio Portfolio, named Portfolio B, of the frontier.

flag_constraints = 0;

% set constraints
const.Aineq = [];
const.bineq = [];

% Portfolio A: Minimum Variance Portfolio
[~,minRisk_P1, minRiskWgt_P1, minRiskRet_P1, minRiskSR_P1] = ...
minRiskPortfolio(names, mean_returns, cov_matrix, risk_free_rate,const,flag_constraints);

% Portfolio B: Maximum Sharpe Ratio Portfolio
[P1,maxSharpeRisk_P1, maxSharpeWgt_P1, maxSharpeRet_P1, maxSharpeSR_P1] = ...
 maxSharpPortfolio(names, mean_returns, cov_matrix, risk_free_rate,const,flag_constraints);

%% market structure

% Create a market structure
mkt = struct();

% Define sectors
mkt.sector.cyclical = ["ConsumerDiscretionary", "Financials", "Materials", "RealEstate", "Industrials"];
mkt.sector.defensive = ["ConsumerStaples", "Utilities", "HealthCare"];
mkt.sector.sensible = ["Energy", "InformationTechnology", "CommunicationServices"];
% Define factors
mkt.factor = ["Momentum", "Value", "Growth", "Quality", "LowVolatility"];

P = Portfolio('AssetList', names);
sectors = [mkt.sector.sensible, mkt.sector.cyclical, mkt.sector.defensive];

sensibleIdx = ismember(P.AssetList, mkt.sector.sensible);
cyclicalIdx = ismember(P.AssetList, mkt.sector.cyclical);
excludeIdx_CS = ismember(P.AssetList, "ConsumerStaples");
excludeIdx_LV = ismember(P.AssetList, "LowVolatility");
sectorIdx = ismember(P.AssetList, sectors);


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

flag_constarints = 1;

% set constraints
const.Aineq = [
    -sensibleIdx;
    cyclicalIdx;
    excludeIdx_CS;
    excludeIdx_LV;
    sectorIdx;
];

const.bineq = [-0.1; 0.3; 0; 0; 0.8];

% Portfolio C and D: Minimum Variance Portfolio and Maximum Sharpe Ratio Portfolio with constraints
[~,minRisk_P2, minRiskWgt_P2, minRiskRet_P2, minRiskSR_P2] = ...
    minRiskPortfolio(names, mean_returns, cov_matrix, risk_free_rate,const,flag_constraints);

[P2,maxSharpeRisk_P2, maxSharpeWgt_P2, maxSharpeRet_P2, maxSharpeSR_P2] = ...
    maxSharpPortfolio(names, mean_returns, cov_matrix, risk_free_rate,const,flag_constraints);


%% 3. Efficient Frontier and Resampling Method

% Portfolio E and F: Minimum Variance Portfolio with resampling
[minRisk_P1_Rsim, minRiskWgt_P1_Rsim, minRiskRet_P1_Rsim, minRiskSR_P1_Rsim, ...
    minRisk_P2_Rsim, minRiskWgt_P2_Rsim, minRiskRet_P2_Rsim, minRiskSR_P2_Rsim, ...
        maxSharpeSR_P1_Rsim, maxSharpeWgt_P1_Rsim, maxSharpeRet_P1_Rsim, maxSharpeRisk_P1_Rsim, ...
            maxSharpeSR_P2_Rsim, maxSharpeWgt_P2_Rsim, maxSharpeRet_P2_Rsim, maxSharpeRisk_P2_Rsim] = ...
                resampling_method(mean_returns, cov_matrix, P1, P2, risk_free_rate,num_assets,names);

%% 4. Black-Litterman Model
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
% Portfolio I: Minimum Variance Portfolio
% Portfolio L: Max Sharpe Ratio
[ptf_I, ptf_L] = ...
        BlackLitterman( prices_2023,...
                        caps,...
                        names,...
                        [view1, view2],...
                        0);

%% 5.  Compute the Maximum Diversified Portfolio and the Maximum Entropy (in asset volatility) Portfolio

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
 cap_wghtd_ptf = capitalizations{1, names} / sum(capitalizations{1, names});

 const.nonlinconstr =  @(weights) customAbsDiffConstraint(weights, cap_wghtd_ptf); % Non linear constraint function



% Portfolio M: Maximum Diversified Portfolio
[weights_m, portfolio_m_return, portfolio_m_std, portfolio_m_SR] = ...
          Max_Diversified_Portfolio(mean_returns, cov_matrix, capitalizations, names, risk_free_rate, mkt, P2,num_assets,const);

% Portfolio N: Maximum Entropy Portfolio
[weights_n, portfolio_n_return, portfolio_n_std, portfolio_n_SR] = ...
    Max_Entropy_Portfolio(mean_returns, cov_matrix, capitalizations, names, risk_free_rate, mkt, P2,num_assets,const);


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
weights_I = array2table(ptf_I.w, 'RowNames', names, 'VariableNames', {'Portfolio_I'});
weights_L = array2table(ptf_L.w, 'RowNames', names, 'VariableNames', {'Portfolio_L'});
weights_M = array2table(weights_m, 'RowNames', names, 'VariableNames', {'Portfolio_M'});
weights_N = array2table(weights_n, 'RowNames', names, 'VariableNames', {'Portfolio_N'});

weightsTable = [weights_A, weights_B, weights_C, weights_D,...
                weights_E, weights_F, weights_G, weights_H,...
                weights_I, weights_L, weights_M, weights_N];

metricsTable = table(...
    [minRiskRet_P1; minRisk_P1; minRiskSR_P1; sum(minRiskWgt_P1)], ...
    [maxSharpeRet_P1; maxSharpeRisk_P1; maxSharpeSR_P1; sum(maxSharpeWgt_P1)], ...
    [minRiskRet_P2; minRisk_P2; minRiskSR_P2; sum(minRiskWgt_P2)], ...
    [maxSharpeRet_P2; maxSharpeRisk_P2; maxSharpeSR_P2; sum(maxSharpeWgt_P2)], ...
    [minRiskRet_P1_Rsim; minRisk_P1_Rsim; minRiskSR_P1_Rsim; sum(minRiskWgt_P1_Rsim)], ...
    [minRiskRet_P2_Rsim; minRisk_P2_Rsim; minRiskSR_P2_Rsim; sum(minRiskWgt_P2_Rsim)], ...
    [maxSharpeRet_P1_Rsim; maxSharpeRisk_P1_Rsim; maxSharpeSR_P1_Rsim; sum(maxSharpeWgt_P1_Rsim)], ...
    [maxSharpeRet_P2_Rsim; maxSharpeRisk_P2_Rsim; maxSharpeSR_P2_Rsim; sum(maxSharpeWgt_P2_Rsim)], ...
    [ptf_I.ret; ptf_I.std; ptf_I.sr; sum(ptf_I.w)], ...
    [ptf_L.ret; ptf_L.std; ptf_L.sr; sum(ptf_L.w)], ...
    [portfolio_m_return; portfolio_m_std; portfolio_m_SR; sum(weights_m)], ...
    [portfolio_n_return; portfolio_n_std; portfolio_n_SR; sum(weights_n)], ...
    'VariableNames', {'Portfolio_A', 'Portfolio_B', 'Portfolio_C', 'Portfolio_D', ...
                      'Portfolio_E', 'Portfolio_F', 'Portfolio_G', 'Portfolio_H', ...
                      'Portfolio_I', 'Portfolio_L', 'Portfolio_M', 'Portfolio_N'}, ...
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

% Display the time taken
disp('==============================================================================================')
disp(['Time taken: ', num2str(toc), ' seconds'])
disp('==============================================================================================')