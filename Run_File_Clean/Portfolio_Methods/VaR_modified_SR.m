function Output_struct = VaR_modified_SR(Ptf, confidence_level, caps, name_ptf)
% This function computes the VaR-modified Sharpe Ratio portfolio
%
% INPUTS:
% Ptf:                  Portfolio object
% confidence_level:     Confidence level for VaR
% caps:                 Capitalizations of the assets
% name_ptf:             Name of the portfolio
%
% OUTPUTS:
% Output_struct:        A struct containing the volatility, weights, return,
%                       and Sharpe ratio of the VaR-modified Sharpe Ratio portfolio


% Extract capitalizations and calculate initial weights
initial_weights = caps / sum(caps);

% VaR-modified Sharpe Ratio optimization
VaR_N = @(x) - (Ptf.AssetMean' * x - ...
    sqrt(x' * Ptf.AssetCovar * x) * norminv(confidence_level));
% Define the objective function (negative VaR-modified Sharpe Ratio)
objective_function_N = @(x) -(Ptf.AssetMean' * x - Ptf.RiskFreeRate) / VaR_N(x);

% Constraints: weights sum to 1 and are between 0 and 1
Aeq = ones(1, length(initial_weights));
beq = 1;
lb = zeros(length(initial_weights), 1);
ub = ones(length(initial_weights), 1);

% Optimization using fmincon
options = optimoptions('fmincon', 'Display', 'None');
weights = fmincon(objective_function_N, initial_weights', [], [], Aeq, beq, lb, ub, [], options);

% Compute the volatility
volatility = sqrt(weights' * Ptf.AssetCovar * weights);

% Compute the return
returns = Ptf.AssetMean' * weights;

% Compute the Sharpe ratio
sharpe_ratio = (returns - Ptf.RiskFreeRate) / volatility;

% Display the minimum risk portfolio
print_portfolio(weights, Ptf.AssetList, returns, volatility, sharpe_ratio, name_ptf);

% Build a struct for the output
Output_struct = struct('Volatility', volatility, 'Weights', weights,...
    'Return', returns, 'Sharpe_Ratio', sharpe_ratio);
Output_struct.Name = name_ptf;
Output_struct.Ptf = Ptf;

end