function [equities,metricesTable] = getEquityandMetrices(Ws, prices, Title, dates)
% This function plots the equity of the portfolios and calculates the
% performance metrics for each portfolio.
%
% INPUTS:
% Ws:       Weights of the portfolios -- Table
% prices:   Prices of the assets
% ptfNames: Names of the portfolios
% Title:    Title of the plot
% dates:    Dates of the assets
%
% OUTPUTS:
% equities:         Equity of the portfolios
% metricesTable:    Table of the performance metrics

% Define a cell array of hexadecimal color codes
hexColors = {
    '#ff0000'; '#e81e63'; '#9c27b0'; 
    '#673ab7'; '#3f51b5'; '#2196f3'; 
    '#03a9f4'; '#00bcd4'; '#009688'; 
    '#4caf50'; '#8bc34a'; '#cddc39'; 
    '#ffeb3b'; '#ffc107'; '#ff9800'; 
    '#ff5722'
};

% Convert hexadecimal codes to RGB triplets
colors = hexToRGB(hexColors);
% Extract portfolio names
ptfNames = Ws.Properties.VariableNames;
% Convert the portfolio weights table into a numeric array
Ws = table2array(Ws);
%  Determine the number of portfolios
numPtfs = size(Ws,2);

% Calculate simple returns from asset prices
ret = prices(2:end, :) ./ prices(1:end-1,:);
% Calculate log returns using continuous compounding
logret = tick2ret(prices, 'Method','continuous');
% Compute the covariance matrix of log returns
cov_matrix = cov(logret);

% Initialize the equity array to store equity values for each portfolio
equities = zeros(size(ret,1), numPtfs);
% Initialize the metrics matrix to store performance metrics for each portfolio
metricsMatrix = zeros(7, numPtfs);

% Create a new figure for plotting equity curves
equity_fig = figure;
set(equity_fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on;

% Loop through each portfolio to compute and plot its equity curve
for col = 1:numPtfs
    % Compute cumulative equity using portfolio returns and weights
    equity = cumprod(ret * Ws(:,col));
    equity = 100 .* equity / equity(1);
    equities(:,col) = equity;

    % Calculate performance metrics for the current portfolio
    [annRet, annVol, Sharpe, MaxDD, Calmar] = getPerformanceMetrics(equity);
    % Diversification Ration
    DR = getDiversificationRatio( Ws(:,col), logret);
    % Entropy
    Entropy = getEntropy( Ws(:,col), cov_matrix);
    
    % Store all metrics in the metrics matrix
    metricsMatrix(:, col) = [annRet; annVol; Sharpe;...
        MaxDD; Calmar; DR; Entropy];

    % Plot the equity curve for the current portfolio
    plot(dates, equity, 'Color', colors(col, :), 'LineWidth', 1.5)
end

% Add legend, labels and title
legend(ptfNames, 'Location','northwest');
xlabel("Dates", 'FontSize', 15);
ylabel("Equity", 'FontSize', 15);
title(Title, 'FontSize', 25);

% Define row names for the performance metrics table
rowNames = {'Annual Return', 'Annual Volatility', 'Sharpe Ratio',...
    'Max Drawdown', 'Calmar Ratio', 'DivRatio', 'Entropy'};

% Create a table to store performance metrics with appropriate row and column names
metricesTable = array2table(metricsMatrix,...
    'RowNames', rowNames,...
    'VariableNames', ptfNames);
end