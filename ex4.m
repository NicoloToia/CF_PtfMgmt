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

% Names
assetNames = capitalizations.Properties.VariableNames(2:end);

% Calculate daily returns for each index in 2023
returns_2023 = diff(log(prices_2023));

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

%% Build the views
% Number of views
v = 2;
tau = 1/length(returns_2023);
P = zeros(v, num_assets);
q = zeros(v, 1);
Omega = zeros(v);
% View on Technology vs. Financials: Given the growing importance of 
% technology, you think that the Technology sector will outperform the 
% Financial sector of the 2% (annual).
P(1, find(assetNames == "InformationTechnology")) = 1;
P(1, find(assetNames == "Financials")) = -1;
q(1) = 0.02;
% View on the Momentum vs. Low Volatility factor: You might assume that 
% Momentum will outperform Low Volatility in a bull market scenario 
% (annual overperformance of 1%)
P(2, find(assetNames == "Momentum")) = 1;
P(2, find(assetNames == "LowVolatility")) = -1;
q(2) = 0.01;
% Fill matrix Omega
Omega(1,1) = tau.*P(1,:)*cov_matrix*P(1,:)';
Omega(2,2) = tau.*P(2,:)*cov_matrix*P(2,:)';
% Change to daily returns
bizyear2bizday = 1/252;
q = q*bizyear2bizday;
Omega = Omega*bizyear2bizday;
% Plot views distribution
X_views = mvnrnd(q, Omega, 20000);
hist3(X_views,'CDataMode','auto','FaceColor','interp','Nbins',[50 50])
%% Capitalizations Weighted PTF
figure;
bar(assetNames,caps);
weightsCaps = caps(1:16)/sum(caps(1:16)); weightsCaps = weightsCaps';
lambda = 1.2;
mu_market = lambda.*cov_matrix*weightsCaps;
cov_market = tau.*cov_matrix;
% plot prior distribution
figure;
X = mvnrnd(mu_market, cov_market, 20000);
histogram(X)
%% Black Litterman Formulas
muBL = inv(inv(cov_market)+P'*inv(Omega)*P)*...
    (P'*inv(Omega)*q + inv(cov_market)*mu_market); 
covBL = inv(P'*inv(Omega)*P + inv(cov_market));
table(assetNames', mu_market*252, muBL*252,...
    'VariableNames', ["AssetNames", "Prior Belief on Exp Ret", "BL ExpREt"])
%% Create Portfolio for Black & Litterman
portBL = Portfolio('NumAssets', num_assets, 'Name', 'MV with BL');
portBL = setDefaultConstraints(portBL);
portBL = setAssetMoments(portBL, muBL, cov_matrix+covBL);
% Estimate Frontier
pwBL = estimateFrontier(portBL, 10000);
[risksBL, retBL] = estimatePortMoments(portBL, pwBL);
portfolio_I = pwBL(:,1);
return_I = retBL(1);
risk_I = risksBL(1);
% Max Sharpe Ratio 
portfolio_L = estimateMaxSharpeRatio(portBL);
[risk_L, return_L] = estimatePortMoments(portBL, portfolio_L);
figure;
plot(risksBL, retBL, 'LineWidth', 2);
hold on;
scatter(risk_I, return_I, 'marker', 'p', 'SizeData', 400,...
    'MarkerFaceColor','green', 'MarkerEdgeColor','k');
scatter(risk_L, return_L, 'marker', 'o', 'SizeData', 200,...
    'MarkerFaceColor','red', 'MarkerEdgeColor','k');
legend('Frontier', 'MVP', 'Max Sharpe Ratio', 'location','southeast');
%% Pie plot of weights
figure;
ax1 = subplot(1,2,1);
idx = portfolio_I > 0.0001;
pie(ax1, portfolio_I(idx), assetNames(idx));
title(ax1, 'Minimum Variance Portfolio', 'Position', [-0.05, 1.6, 0]);

ax2 = subplot(1,2,2);
idx = portfolio_L > 0.0001;
pie(ax2, portfolio_L(idx), assetNames(idx));
title(ax2, "Max Sharpe Ratio", 'Position', [-0.05, 1.6, 0]);

%% Stacked weigths Frontier
figure;
ba = bar(pwBL','stacked', 'FaceColor','flat');
ba(1).CData = [0.3 0.3 0.7];
ba(2).CData = [1 1 1]*0.8;
legend(assetNames)