%% 3)
clc
clear all
close all
warning off


% Load data from Excel files
prices = readtable('prices.xlsx');
capitalizations = readtable('capitalizations.xlsx');

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

% Define the risk-free rate (annual, assuming 4% risk-free rate)
risk_free_rate = 0.04 / 365; % Convert to daily (365 is correct)

% Initialize a Portfolio object with the list of asset names
p_1 = Portfolio('AssetList', assetNames);
p_2 = Portfolio('AssetList', assetNames);


% Set default constraints for the portfolio: weights sum to 1
% and there are no short positions (non-negative weights)
P_1 = setDefaultConstraints(p_1);
P_2 = setDefaultConstraints(p_2);



% Define constraints based on problem statement
lb = zeros(num_assets, 1);  % Lower bounds on weights
ub = ones(num_assets, 1);   % Upper bounds on weights
Aeq = ones(1, num_assets);  % Equality constraint for sum of weights
beq = 1;

% Set Consumer Staples and Low Volatility weights to 0
% Assuming Consumer Staples is asset #8 and Low Volatility is asset #16
lb(8) = 0; ub(8) = 0; % Consumer Staples
lb(16) = 0; ub(16) = 0; % Low Volatility

% Additional constraints for sectors:
% Sensible sectors: e.g., Health Care, Utilities, Consumer Staples
% Cyclical sectors: e.g., Consumer Discretionary, Industrials, Energy

SectorConstraints = zeros(11,16);
for i =1:11
    SectorConstraints(i,i) = 1;
end
SectorConstraints = zeros(1,16);
SectorConstraints(1:11) = 1;
Aineq = [
    % Sensible sectors total weight >= 10%
    -[1 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0];
    
    % Cyclical sectors total weight <= 30%
    [0 1 0 1 0 1 0 0 0 1 1 0 0 0 0 0];
    % Consumer staple at 0
    [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
    % Low volatility at 0
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    % Maximum exposure in any sector ≤ 80%
    SectorConstraints;
];


bineq = [
    -0.1;    % Sensible sectors >= 10%
    0.3;     % Cyclical sectors <= 30%
    0;
    0;
    %0.8 * ones(11, 1)  % Max exposure for each sector
    0.8
];

P_2.AInequality = Aineq;
P_2.bInequality = bineq;
P_2.AEquality = Aeq;
P_2.bEquality = beq;

for n = 1:50

    R = mvnrnd(mean_returns, cov_matrix);
    New_mean_returns = R;
    NewCov = iwishrnd(cov_matrix, num_assets);


    Psim_1 = setAssetMoments(P_1, New_mean_returns, NewCov);
    Psim_2 = setAssetMoments(P_2, New_mean_returns, NewCov);


    w_sim_1 = estimateFrontier(Psim_1, 100);
    w_sim_2 = estimateFrontier(Psim_2, 100);


    [pf_riskSim_1, pf_RetnSim_1] = estimatePortMoments(Psim_1, w_sim_1);
    [pf_riskSim_2, pf_RetnSim_2] = estimatePortMoments(Psim_2, w_sim_2);


    RetPtfSim_1(:,n) = pf_RetnSim_1;
    RiskPtfSim_1(:, n) = pf_riskSim_1;
    Weights_1(n,:,:) = w_sim_1;

    RetPtfSim_2(:,n) = pf_RetnSim_2;
    RiskPtfSim_2(:, n) = pf_riskSim_2;
    Weights_2(n,:,:) = w_sim_2;

end


pwgt_1 = mean(Weights_1,3);
pf_risk_1 = mean(RiskPtfSim_1,1);
pf_Retn_1 = mean(RetPtfSim_1,1);

pwgt_2 = mean(Weights_2,3);
pf_risk_2 = mean(RiskPtfSim_2,1);
pf_Retn_2 = mean(RetPtfSim_2,1);



% Display the return of the portfolio with the minimum risk on the frontier
ptf_E_return = pf_Retn_1(pf_risk_1 == min(pf_risk_1));
ptf_G_return = pf_Retn_2(pf_risk_2 == min(pf_risk_2));


% Display the minimum risk (volatility) of the portfolios on the frontier
ptf_E_vol = min(pf_risk_1);
ptf_G_vol = min(pf_risk_2);


% Display the weights of the first portfolio on the frontier
weights_E = pwgt_1(:,1);
weights_G = pwgt_2(:,1);


% Calculate the Sharpe Ratio for each portfolio on the frontier
sr_1 = (pf_Retn_1 - risk_free_rate) ./ pf_risk_1;
sr_2 = (pf_Retn_2 - risk_free_rate) ./ pf_risk_2;


% Display the maximum Sharpe Ratio, representing the portfolio with the best risk-adjusted return
ptf_F_sharpe = max(sr_1);
ptf_H_sharpe = max(sr_2);


% Display results
disp('Portfolio E (Minimum Variance Robust Portfolio) USING THE OBJECT PORTFOLIO:');
disp(['Expected Return: ', num2str(ptf_E_return)]);
disp(['Volatility: ', num2str(ptf_E_vol)]);
disp(['Sharpe Ratio: ', num2str(sr_1(pf_risk_1 == min(pf_risk_1)))]);
disp('Weights:');
disp(weights_E');

disp('Portfolio F (Maximum Sharpe Ratio Robust Portfolio):');
disp(['Expected Return: ', num2str(pf_Retn_1(sr_1==ptf_F_sharpe))]);
disp(['Volatility: ', num2str(pf_risk_1(sr_1==ptf_F_sharpe))]);
disp(['Sharpe Ratio: ', num2str(ptf_F_sharpe)]);
disp('Weights:');
disp(pwgt_1(:,find((sr_1==ptf_F_sharpe)))');

% Display results
disp('Portfolio G (Minimum Variance Robust Portfolio) USING THE OBJECT PORTFOLIO:');
disp(['Expected Return: ', num2str(ptf_G_return)]);
disp(['Volatility: ', num2str(ptf_G_vol)]);
disp(['Sharpe Ratio: ', num2str(sr_2(pf_risk_2 == min(pf_risk_2)))]);
disp('Weights:');
disp(weights_G');

disp('Portfolio H (Maximum Sharpe Ratio Robust Portfolio):');
disp(['Expected Return: ', num2str(pf_Retn_2(sr_2==ptf_H_sharpe))]);
disp(['Volatility: ', num2str(pf_risk_2(sr_2==ptf_H_sharpe))]);
disp(['Sharpe Ratio: ', num2str(ptf_H_sharpe)]);
disp('Weights:');
disp(pwgt_2(:,find((sr_1==ptf_H_sharpe)))');



