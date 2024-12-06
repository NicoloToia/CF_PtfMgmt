function plotFrontierWeights(volatility, weights, names)
% Function to plot how the weights of the efficient frontier evolve for
% different levels of volatility
%
% INPUTS
% volatility:   Vector of volailities of the efficient frontier
% weights:      Matrix of the weights of the portfolios belonging to the
%               efficient frontier

% DA MATCHARE CON I COLORI DI EDO:
customColors = [
    0.50, 0.00, 0.50; % Magenta (InformationTechnology)
    0.40, 0.80, 0.67; % Teal (Financials)
    0.96, 0.36, 0.43; % Pink (HealthCare)
    0.50, 0.00, 0.50; % Magenta (ConsumerDiscretionary)
    0.30, 0.75, 0.93; % Cyan (CommunicationServices)
    0.47, 0.67, 0.19; % Lime (Industrials)
    0.64, 0.08, 0.18; % Maroon (ConsumerStaples)
    0.00, 0.45, 0.74; % Blue (Energy)
    0.98, 0.41, 0.07; % Red-Orange (Utilities)
    0.47, 0.25, 0.80; % Purple (RealEstate)
    0.85, 0.33, 0.10; % Orange (Materials)
    0.93, 0.69, 0.13; % Gold (Momentum)
    1.00, 0.84, 0.00; % Yellow (Value)
    0.20, 0.63, 0.17; % Green (Growth)
    0.25, 0.25, 0.25; % Dark Gray (Quality)
    0.53, 0.81, 0.92; % Light Blue (LowVolatility)
];


figure;
hold on;
% Plot each asset's weights
for i = 1:size(weights,1)
    plot(volatility, weights(i, :), '--', 'LineWidth', 2, 'Color', customColors(i, :));
end

% Add labels, legend, and title
xlabel('Volatility', 'FontSize', 15);
ylabel('Portfolio Weights', 'FontSize', 15);
title('Weights of Efficient Frontier Portfolios', 'FontSize', 25);
legend(names, 'Location', 'northwest', 'FontSize', 10, 'NumColumns', 2);
grid on;
hold off;

end