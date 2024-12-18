function plotFrontierWeights(volatility, weights, names, desiredOrder)
% Function to plot how the weights of the efficient frontier evolve for
% different levels of volatility
%
% INPUTS
% volatility:   Vector of volailities of the efficient frontier
% weights:      Matrix of the weights of the portfolios belonging to the
%               efficient frontier
% names:        Asset's names
% desiredOrder: Desired order of the assets

% Define colors for the plot
hexColors = {
    '#AF1740', '#1F4529', '#FFC436', '#526E48', '#DE7C7D', ...
    '#C2FFC7', '#F7E987', '#CC2B52', '#F8DE22', '#9EDF9C', ...
    '#62825D', '#133E87', '#4A628A', '#7AB2D3', '#B9E5E8', ...
    '#DFF2EB'};
colours = hexToRGB(hexColors);
% Order the colours grouping them by sector
[~, order] = ismember(desiredOrder, names);
weights = weights(order, :); % Reorder rows of weights matrix
colours = colours(order, :); % Reorder colors accordingly
names = names(order);        % Reorder asset names

% Stacked Area Plot
figure;
areaHandle = area(volatility, weights', 'EdgeColor', 'none');
% Assign the custom colours to the areas
for i = 1:size(weights,1)
    areaHandle(i).FaceColor = colours(i, :); % Apply the i-th color from the RGB matrix
end

% Adding white lines for separation at the boundary between areas
cumulativeWeights = cumsum(weights, 1); % Cumulative weights to get boundaries
hold on;
for i = 1:size(weights,1)
    % Plot white line at the boundary between areas
    plot(volatility, cumulativeWeights(i, :), 'w-', 'LineWidth', 1);
end
hold off;

% Add labels, title and legend
xlim([min(volatility) max(volatility)]);
ylim([0 1]);
xlabel('Volatility', 'FontSize', 15);
ylabel('Portfolio Weights', 'FontSize', 15);
title('Weights of Efficient Frontier Portfolios', 'FontSize', 25);
legend(names, 'Location', 'bestoutside', 'FontSize', 10, 'NumColumns', 1);
grid on;

end