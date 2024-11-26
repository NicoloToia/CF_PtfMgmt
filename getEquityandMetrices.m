function [equities,metricesTable] = getEquityandMetrices(   Ws,...
                                                            prices,...
                                                            ptfNames,...
                                                            Title)

% GETEQUITYANDMETRICES plots the equity of the portfolios
%
% INPUTS:
% Ws: weights of the portfolios
% prices: prices of the assets
% ptfNames: names of the portfolios
% Title: title of the plot

numPtfs = size(Ws,2);
ret = prices(2:end, :) ./ prices(1:end-1,:);
equities = zeros(size(ret,1), numPtfs);
figure;
hold on;
for col = 1:numPtfs
    equity = cumprod(ret * Ws(:,col));
    equity = 100 .* equity / equity(1);
    equities(:,col) = equity;
    plot(equity, 'LineWidth', 3)
end
legend(ptfNames);
metricesTable = 0;
title(Title);
end

