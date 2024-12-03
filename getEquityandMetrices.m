function [equities,metricesTable] = getEquityandMetrices(   Ws,...
                                                            prices,...
                                                            Title)

% GETEQUITYANDMETRICES plots the equity of the portfolios
%
% INPUTS:
% Ws: weights of the portfolios -- Table
% prices: prices of the assets
% ptfNames: names of the portfolios
% Title: title of the plot
% Define a cell array of hexadecimal color codes
hexColors = {
    '#0000EE', '#00EE00', '#EE0000', '#EE00EE', '#EEBE00', ...
    '#F78536', '#000000', '#FFCF83', '#CF56A1', '#A66E38', ...
    '#ADD8E6', '#FEBBCC', '#FF0060', '#8B9A46'};

% Convert hexadecimal codes to RGB triplets
colors = hexToRGB(hexColors);

ptfNames = Ws.Properties.VariableNames;
Ws = table2array(Ws);
numPtfs = size(Ws,2);
ret = prices(2:end, :) ./ prices(1:end-1,:);
equities = zeros(size(ret,1), numPtfs);
metricsMatrix = zeros(5, numPtfs);
figure;
hold on;
for col = 1:numPtfs
    equity = cumprod(ret * Ws(:,col));
    equity = 100 .* equity / equity(1);
    equities(:,col) = equity;
    [annRet, annVol, Sharpe, MaxDD, Calmar] = getPerformanceMetrics(equity);
    metricsMatrix(:, col) = [annRet; annVol; Sharpe; MaxDD; Calmar];
    plot(equity, 'Color', colors(col, :), 'LineWidth', 3 )
end
legend(ptfNames);
title(Title);
rowNames = {'annRet', 'annVol', 'Sharpe', 'MaxDD', 'Calmar'};
metricesTable = array2table(metricsMatrix,...
    'RowNames', rowNames,...
    'VariableNames', ptfNames);

end

