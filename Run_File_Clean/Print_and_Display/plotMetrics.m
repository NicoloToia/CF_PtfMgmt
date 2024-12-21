function plotMetrics(performancesMetrics)
% This function creates a clustered bar chart to compare the performance
% metrics of different portfolios
%
% INPUT:
% performancesMetrics: Performance metrics table

% Filter the data to display only selected metrics
selected_metrics = {'Annual Return', 'Annual Volatility', 'Sharpe Ratio', 'DivRatio', 'Entropy'};
filtered_data = performancesMetrics(selected_metrics, :);

% Extract the data to plot
metrics = filtered_data.Properties.RowNames; 
portfolios = filtered_data.Properties.VariableNames; 
data = table2array(filtered_data); 

% Transpose the data for a better visualization
data_transposed = data';

% Create a clustered bar chart
figure('Color', 'w', 'Position', [100, 100, 1000, 600]); 
b = bar(data_transposed, 'grouped'); 

colors = lines(length(metrics)); 
for k = 1:length(b)
    b(k).FaceColor = colors(k, :); 
    b(k).EdgeColor = 'none'; 
end

% Add labels and title
xlabel('Portfolios', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Value', 'FontSize', 14, 'FontWeight', 'bold');
title('Performance Metrics Benchmarking', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.3, 0.3, 0.3]);

% Force the x-axis to display all portfolios
set(gca, 'XTick', 1:length(portfolios));
set(gca, 'XTickLabel', portfolios, 'FontSize', 12, 'FontWeight', 'bold');
xtickangle(45);

% Add legend with metric names
legend(metrics, 'Location', 'NorthWest', 'FontSize', 12, 'FontWeight', 'bold', 'Box', 'off');

grid on;
set(gca, 'GridLineStyle', ':', 'LineWidth', 0.7); % dashed line for the grid
box on;
ylim([0 max(data_transposed(:)) * 1.1]); % y-axis adjusted limits 
    
end
