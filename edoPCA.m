clear;
close all; 
clc;

%% PART A 

% Load data from Excel files
prices = readtable('prices.xlsx');
capitalizations = readtable('capitalizations.xlsx');
caps = capitalizations{1,2:end};

% Extract dates and data
dates = prices{:,1}; % First column contains dates
prices_data = prices{:,2:end}; % Data starts from the second column

% Convert dates to MATLAB date format if needed
dates = datetime(dates);

% Filter prices for 2023
start_date = datetime(2023,1,1);
end_date = datetime(2023,12,31);
prices_2023 = prices_data(dates >= start_date & dates <= end_date, :);
start_test = datetime(2024, 1, 1);
prices_2024 = prices_data(dates >= start_test, :);

% Names
assetNames = capitalizations.Properties.VariableNames(2:end);

% Calculate daily returns for each index in 2023
returns_2023 = tick2ret(prices_2023);

% Calculate mean returns and covariance matrix
mean_returns = mean(returns_2023)';
cov_matrix = cov(returns_2023);

% Define number of assets
num_assets = length(mean_returns);

cyclical = ["ConsumerDiscretionary", "Financials", "Materials", "RealEstate", "Industrials"];
defensive = ["ConsumerStaples", "Utilities", "HealthCare"];
sensible = ["Energy", "InformationTechnology", "CommunicationServices"];
factor = ["Momentum","Value","Growth","Quality","LowVolatility"];
% Groups vector
groups = {cyclical, defensive, sensible, factor};


%% 6.
% Compute the portfolio (Portfolio P), using the Principal Component 
% Analysis, that maximizes its expected return under the following 
% constraints (to be considered all at once):
% • Standard constraints,
% • The volatility of the portfolio has to be equal or less than a target
% volatility of σtgt = 0.1
% You have to use the minimum number of factors that explains more
% than the 90% of the cumulative variance.

ret = tick2ret(prices_2023, 'Method','continuous');
mean_ret = mean(ret);
std_ret = (ret - mean_ret) ./ std(ret);

k = 16;
[~,~,latent,~,~,~] =...
    pca(std_ret, 'NumComponents',k);
% coeff -> loadings
% score -> scores matrix
% latent
% explained -> Percentage of total variance explained column vector
% Find cumulative explained variance
TotVar = sum(latent);
explainedVar = latent(1:k)/TotVar;
CumExplainedVar = cumsum(explainedVar);
find(CumExplainedVar >= .9, 1)
% Choosing k = 6
%%
h = figure();
subplot(1,2,1)
bar(explainedVar)
title("% Explained Variance");
subplot(1,2,2)
bar(CumExplainedVar)
title("% Cumulative Variance");
%%
k = 6;
[factorLoading,factorRetn,latent,r,explained,mu] =...
    pca(std_ret, 'NumComponents',k);
covarFactor = cov(factorRetn);
%% PTF OPTIMIZATION
% Reconstruct 
reconReturn = factorRetn * factorLoading' + mean_ret;
unexplainedRetn = ret - reconReturn; % IDIOSYNCRATHIC PART (epsilon)
% 
unexplainedCov = diag(cov(unexplainedRetn));
D = diag(unexplainedCov);
covarAsset = factorLoading*covarFactor*factorLoading' + D;

targetRisk = 0.1;  % Standard deviation of portfolio return
tRisk = targetRisk*targetRisk;  % Variance of portfolio return

optimProb = optimproblem('Description',...
    'Portfolio with factor covariance matrix','ObjectiveSense','max');
wgtAsset = optimvar('asset_weight', num_assets, 1,...
    'Type', 'continuous', 'LowerBound', 0, 'UpperBound', 1);
wgtFactor = optimvar('factor_weight', k, 1, 'Type', 'continuous');

optimProb.Objective = sum(mean_ret*wgtAsset);

optimProb.Constraints.asset_factor_weight = ...
    factorLoading'*wgtAsset - wgtFactor == 0;
optimProb.Constraints.risk = ...
    wgtFactor'*covarFactor*wgtFactor + ...
    wgtAsset'*D*wgtAsset <= tRisk;
optimProb.Constraints.budget = sum(wgtAsset) == 1;

x0.asset_weight = ones(num_assets, 1)/num_assets;
x0.factor_weight = zeros(k, 1);
opt = optimoptions("fmincon", "Algorithm","sqp", "Display", "off", ...
    'ConstraintTolerance', 1.0e-8, 'OptimalityTolerance', 1.0e-8, 'StepTolerance', 1.0e-8);
x = solve(optimProb,x0, "Options",opt);
portfolio_P = x.asset_weight;
factorWgt1 = x.factor_weight;
%%
portfolio_EW = ones(num_assets,1)/num_assets;
equity_P = getEquityandMetrices([portfolio_P portfolio_EW], prices_2023, ...
    ["P", "EW"], "2023");
figure;
pie(portfolio_P (portfolio_P >= 0.001), assetNames(portfolio_P >= 0.001))
T2 = table(assetNames', round(portfolio_P,4))

%% Using PTF Object
covarAsset = factorLoading*covarFactor*factorLoading'+D;
port = Portfolio("AssetMean", mean_ret, 'AssetCovar', round(covarAsset,13),...
    'LowerBound', 0, 'UpperBound', 1, ...
    'Budget', 1);
portfolio_P = estimateFrontierByRisk(port, targetRisk);
T2 = table(assetNames', round(portfolio_P,4))
portfolio_EW = ones(num_assets,1)/num_assets;
equity_P = getEquityandMetrices([portfolio_P portfolio_EW], prices_2023, ...
    ["P", "EW"], "2023");
figure
pie(portfolio_P (portfolio_P >= 0.01), assetNames(portfolio_P >= 0.01))


%% fmincon trial

func = @(x) -((mean_ret*x)- ...
             ((factorLoading'*x)'*covarFactor*(factorLoading'*x) + ...
             x'*D*x));
func = @(x) -((mean_ret*x));
lb = zeros(1,num_assets);
ub = ones(1,num_assets);
Aeq = ones(1,num_assets);
beq = 1;
x0 = ones(num_assets, 1)/num_assets;
nonlcon =  @(x) nonlconPCA(x, ...
                           factorLoading, ...
                           D, ...
                           covarFactor,...
                           0.1); % Non linear constraint function

options = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e8,...
    'Algorithm','sqp','MaxIterations',1e5);
[portfolio_P, fval] = fmincon(func, x0, [],[],Aeq, beq, lb, ub,nonlcon,...
    options);
T2 = table(assetNames', round(portfolio_P,4))
equity_P = getEquityandMetrices([portfolio_P portfolio_EW], prices_2023, ...
    ["P", "EW"], "2023");
figure;
pie(portfolio_P (portfolio_P >= 0.01), assetNames(portfolio_P >= 0.01))


