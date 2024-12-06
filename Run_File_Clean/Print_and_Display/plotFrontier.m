function [] = plotFrontier(volatilities, returns, portfolios_vector, title_text, legend_text)
% Function to plot the efficient frontiers.
%
% INPUTS
% volatilities:
% returns:
% portfolios_vector:

figure()
for i=1:size(volatilities, 2)
    plot(volatilities(:,i), returns(:,i), 'LineWidth', 2)
    hold on
end
for j=1:size(portfolios_vector, 1)
    plot(portfolios_vector(j,1), portfolios_vector(j,2), '*', 'LineWidth', 2)
end
xlabel('Volatility', 'FontSize', 15)
ylabel('Return', 'FontSize', 15)
title(title_text, 'FontSize', 25);
legend(legend_text, 'FontSize', 25, 'Location', 'southeast')
hold off


end