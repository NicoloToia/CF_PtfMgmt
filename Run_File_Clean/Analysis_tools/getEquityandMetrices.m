function [equities,metricesTable] = getEquityandMetrices(Ws, prices, Title)

% This function plots the equity of the portfolios and calculates the
% performance metrics for each portfolio.
%
% INPUTS:
% Ws: weights of the portfolios -- Table
% prices: prices of the assets
% ptfNames: names of the portfolios
% Title: title of the plot
%
% OUTPUTS:
% equities: equity of the portfolios
% metricesTable: table of the performance metrics

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

ptfNames = Ws.Properties.VariableNames;
Ws = table2array(Ws);
numPtfs = size(Ws,2);
ret = prices(2:end, :) ./ prices(1:end-1,:);
logret = tick2ret(prices, 'Method','continuous');
equities = zeros(size(ret,1), numPtfs);
metricsMatrix = zeros(7, numPtfs);
equity_fig = figure;
set(equity_fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on;
for col = 1:numPtfs
    equity = cumprod(ret * Ws(:,col));
    equity = 100 .* equity / equity(1);
    equities(:,col) = equity;
    [annRet, annVol, Sharpe, MaxDD, Calmar] = getPerformanceMetrics(equity);
    % Diversification Ration
    DR = getDiversificationRatio( Ws(:,col), logret);
    % Entropy
    Entropy = getEntropy( Ws(:,col));

    metricsMatrix(:, col) = [annRet; annVol; Sharpe;...
        MaxDD; Calmar; DR; Entropy];
    plot(equity, 'Color', colors(col, :), 'LineWidth', 1.5)
end
legend(ptfNames, 'Location','northwest');
title(Title);
% rowNames = {'annRet', 'annVol', 'Sharpe', 'MaxDD', 'Calmar'};
rowNames = {'Annual Return', 'Annual Volatility', 'Sharpe Ratio',...
    'Max Drawdown', 'Calmar Ratio', 'DivRatio', 'Entropy'};
metricesTable = array2table(metricsMatrix,...
    'RowNames', rowNames,...
    'VariableNames', ptfNames);