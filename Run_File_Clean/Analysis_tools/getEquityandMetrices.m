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
    '#0000EE', '#00EE00', '#EE0000', '#EE00EE', '#EEBE00', ...
    '#F78536', '#000000', '#FFCF83', '#CF56A1', '#A66E38', ...
    '#ADD8E6', '#FEBBCC', '#FF0060', '#8B9A46', '#FFD744'};

% Convert hexadecimal codes to RGB triplets
colors = hexToRGB(hexColors);

% Get the names of the portfolios, the weights and the number of portfolios
ptfNames = Ws.Properties.VariableNames;
Ws = table2array(Ws);
numPtfs = size(Ws,2);

% Calculate the returns
ret = prices(2:end, :) ./ prices(1:end-1,:);

% Initialize the equity and metrics matrices
equities = zeros(size(ret,1), numPtfs);
metricsMatrix = zeros(5, numPtfs);

figure;
hold on;
for col = 1:numPtfs
    
    % Plot Portfolios Performance 
    equity = cumprod(ret * Ws(:,col));
    equity = 100 .* equity / equity(1);
    equities(:,col) = equity;
    [annRet, annVol, Sharpe, MaxDD, Calmar] = getPerformanceMetrics(equity);
    metricsMatrix(:, col) = [annRet; annVol; Sharpe; MaxDD; Calmar];
    plot(equity, 'Color', colors(col, :), 'LineWidth', 3 )

end

% Set the plot properties
legend(ptfNames, 'Location', 'northwest');
title(Title);
xlabel('Time (daily)');
ylabel('Portfolio Value (USD)');
hold off;

% Create a table of the performance metrics
rowNames = {'Annual Return', 'Annual Volatility', 'Sharpe Ratio', 'Max Drawdown', 'Calmar Ratio'};
metricesTable = array2table(metricsMatrix,...
    'RowNames', rowNames,...
    'VariableNames', ptfNames);

end