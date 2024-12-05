clc
clear all
close all
warning off
% Load data (Assumes you have 'capitalizations.xlsx' and 'prices.xlsx')
capitalizations = readtable('capitalizations.xlsx', 'Sheet', 'Sheet1');
prices = readtable('prices.xlsx', 'Sheet', 'Sheet1');
% Extract prices and calculate daily returns
dates = prices{:, 1};
prices_data = prices{:, 2:end};
start_date = datetime(2023,1,1);
end_date = datetime(2023,12,31);
prices_2023 = prices_data(dates >= start_date & dates <= end_date, :);
returns = tick2ret(prices_2023);
% Compute expected returns (mean of daily returns)
mean_ret = mean(returns);
% Compute the variance-covariance matrix
cov_matrix = cov(returns);
% Extract capitalizations and calculate initial weights
capitalization_values = capitalizations{1, 2:end}';
initial_weights = capitalization_values / sum(capitalization_values);
% Risk-free rate (annualized to daily)
risk_free_rate = 0;
% Confidence level for VaR
confidence_levels = [0.95 0.99];
ptfs = [];

for i = 1:2
    confidence_levels(i)
    VaR_H = @(x) - quantile(returns * x, 1 - confidence_levels(i));
    VaR_N = @(x)  - mean_ret * x + sqrt(x' * cov_matrix * x) * norminv(confidence_levels(i));
    % Define the objective function (negative VaR-modified Sharpe Ratio)
    objective_function_H = @(x) -(mean_ret * x - risk_free_rate) / VaR_H(x);
    objective_function_N = @(x) -(mean_ret * x - risk_free_rate) / VaR_N(x);
    % Constraints: weights sum to 1 and are between 0 and 1
    Aeq = ones(1, length(initial_weights));
    beq = 1;
    lb = zeros(length(initial_weights), 1);
    ub = ones(length(initial_weights), 1);
    % Optimization using fmincon
    options = optimoptions('fmincon', 'Display', 'None');
    w_H = fmincon(objective_function_H, ...
                  initial_weights, [], [], ...
                  Aeq, beq, lb, ub, [], options);
    w_N = fmincon(objective_function_N, ...
                  initial_weights, [], [], ...
                  Aeq, beq, lb, ub, [], options);
    ptfs = [ptfs w_H]; ptfs = [ptfs w_N];
end
% Asset names (assuming they are the column headers in the prices file)
names = prices.Properties.VariableNames(2:end);
ptfs = [ptfs ones(16,1)/16];
rets = returns * ptfs;
figure;
hold on;
histogram(rets(:,1), 20)
histogram(rets(:,2), 20)
histogram(rets(:,3), 20)
histogram(rets(:,4), 20)
histogram(rets(:,5), 20)
legend('95_H','95_N','99_H','99_N','EW')
% Display results in the specified format
disp('===========================================================================')
disp(' VaR modified sharpe ratio Portfolio (Portfolio Q)')
disp('===========================================================================')
weightsTable = array2table(ptfs, "RowNames",names,...
    "VariableNames",{'95_H','95_N','99_H','99_N','EW'})
[~, table] = getEquityandMetrices(weightsTable, prices_data, "2023");
figure;
subplot(2,2,1);
ws = ptfs(:,1);
pie(ws(ws>=1/1000),names(ws>=1/1000))
subplot(2,2,2);
ws = ptfs(:,2);
pie(ws(ws>=1/100),names(ws>=1/100))
subplot(2,2,3);
ws = ptfs(:,3);
pie(ws(ws>=1/100),names(ws>=1/100))
subplot(2,2,4);
ws = ptfs(:,4);
pie(ws(ws>=1/100),names(ws>=1/100))
table