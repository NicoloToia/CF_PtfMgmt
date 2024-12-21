function [] = plot_frontier(volatilities, returns, portfolios_vector, title_text, legend_text)
% Function to plot the efficient frontiers and some portfolios on the same plot.
%
% INPUTS:
% volatilities:         Matrix in which the i-th column represents the
%                       volatilities' values of the i-th frontier
% returns:              Matrix in which the i-th column represents the
%                       returns' values of the i-th frontier              
% portfolios_vector:    Vector containing the portfolios structs that we
%                       want to plot. In particular, each element of the
%                       vector is a struct with fields Volatility and
%                       Return
% title_text:           String containing the text of the title of the plot
% legend_text:          Cell containing the strings with the legend entries


% Plot the frontier:
figure()
for i=1:size(volatilities, 2)
    plot(volatilities(:,i), returns(:,i), 'LineWidth', 2)
    hold on
end

% Plot the portfolios:
for j=1:length(portfolios_vector)
    plot(portfolios_vector(j).Volatility, portfolios_vector(j).Return, '*', 'LineWidth', 2)
end

% Add axis labels, title and legend
xlabel('Volatility', 'FontSize', 15)
ylabel('Return', 'FontSize', 15)
title(title_text, 'FontSize', 25);
legend(legend_text{:}, 'FontSize', 15, 'Location', 'southeast')
grid on
hold off


end