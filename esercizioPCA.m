clear;
close all; 
clc;
warning on;
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
returns_2023 = tick2ret(prices_2023, 'Method','continuous');

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
n_samples = size(ret, 1);
mean_ret = mean(ret);
sd_ret = std(ret);
std_ret = (ret - repmat(mean_ret,n_samples,1)) ./ repmat(sd_ret,n_samples,1);

k = 16;
[~,~,latent,~,~,~] =...
    pca(std_ret, 'NumComponents',k);
% coeff -> loadings
% score -> scores matrix
% latent
% explained -> Percentage of total variance explained column vector
% Find cumulative explained variance
TotVar = sum(latent);
explainedVar = latent./TotVar;
CumExplainedVar = cumsum(explainedVar);
find(CumExplainedVar >= .9, 1)
% Choosing k = 6
%%
h = figure();
cmap = parula(length(explainedVar));
subplot(1,2,1)
b = bar(explainedVar,'FaceColor','white','EdgeColor','k');
for i = 1:length(explainedVar)
    b.FaceColor = 'flat'; % Enable individual bar coloring
    b.CData(i, :) = cmap(i, :); % Assign color to each bar
end
title("% Explained Variance");
subplot(1,2,2)
b = bar(CumExplainedVar,'FaceColor','white','EdgeColor','k');
for i = 1:length(explainedVar)
    b.FaceColor = 'flat'; % Enable individual bar coloring
    b.CData(i, :) = cmap(i, :); % Assign color to each bar
end
yline(0.90, '--k', 'LineWidth',1.5);
title("% Cumulative Variance");
%%
k = 6;
[factorLoading,factorRetn,latent,r,explained,mu] =...
    pca(std_ret, 'NumComponents',k);
covarFactor = cov(factorRetn);
%% PTF OPTIMIZATION
% Reconstruct 
samples = size(ret,1);
reconReturn = (factorRetn * factorLoading') .* repmat(sd_ret,samples,1) + ...
    repmat(mean_ret,samples,1);
unexplainedRetn = ret - reconReturn; % IDIOSYNCRATHIC PART (epsilon)
% 
unexplainedCov = diag(cov(unexplainedRetn));
D = diag(unexplainedCov);
covarAsset = factorLoading*covarFactor*factorLoading' + D;

%% Using PTF Object
port = Portfolio("AssetMean", mean_ret, 'AssetCovar', round(covarAsset,13),...
    'LowerBound', 0, 'UpperBound', 1, ...
    'Budget', 1);
portfolio_P = estimateFrontierByRisk(port, 0.7);
T2 = table(assetNames', round(portfolio_P,4))
figure;
pie(portfolio_P (portfolio_P >= 0.01), assetNames(portfolio_P >= 0.01))

%%
% fmincon

func = @(x) -(mean_ret*x);
Aeq = ones(1, num_assets);
beq = 1;
lb = zeros(num_assets, 1);
ub = ones(num_assets, 1);
x0 = rand(num_assets, 1); x0 = x0/sum(x0);


options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
    'MaxFunctionEvaluations', 1e6, 'MaxIterations',1e6,...
    'ConstraintTolerance', 1e-2, 'Display', 'none');

[w, fval] = fmincon(func, x0, [], [], Aeq, beq, lb, ub,...
    @(x) nonlconPCA(x, 0.7, factorLoading, covarFactor, D),...
    options);

ptfs = [portfolio_P ones(16,1)/16 caps'/sum(caps)];
% Display results in the specified format
disp('===========================================================================')
disp(' PCA Portfolio (Portfolio P)')
disp('===========================================================================')
weightsTable = array2table(round(ptfs,4), "RowNames",assetNames,...
    "VariableNames",{'PCA','EW','Caps'})
[~, table] = getEquityandMetrices(weightsTable, prices_2023, "2023")
[~, table] = getEquityandMetrices(weightsTable, prices_2024, "2024")

