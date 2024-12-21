function Output_struct = Max_Diversified_Portfolio(Ptf, const, name_ptf)
% This function computes the Maximum Diversified Portfolio under
% the assigned constraints.
% 
% INPUTS
% Ptf:      Portfolio object
% const:    Structure with the constraints
% name_ptf: Name of the portfolio
%
% OUTPUTS
% Output_struct: A struct containing the volatility, weights, return,
%                and Sharpe ratio of the Maximum Diversified Portfolio

% Set up optimization problem
% Initial guess
initial_guess = ones(Ptf.NumAssets, 1) / Ptf.NumAssets;

% Options for optimization
options = optimoptions('fmincon', ...     
        'Algorithm', 'sqp', ...          % Specify the algorithm     
        'StepTolerance', 1e-6, ...       % Smaller than default StepTolerance     
        'Display', 'off');               % Show iteration information

% Diversity Ratio function
diversification_ratio = @(w) -log(w' * sqrt(diag(Ptf.AssetCovar)) / sqrt(w' * Ptf.AssetCovar * w));

% Optimization
[weights, ~] = fmincon(diversification_ratio, initial_guess, const.A, const.b,...
                                 const.Aeq, const.beq, const.lb, const.ub, const.nonlinconstr, options);
                  
% Compute the return, volatility, and Sharpe ratio of the portfolio                                        
portfolio_return = Ptf.AssetMean' * weights;
portfolio_std = sqrt(weights' * Ptf.AssetCovar * weights);
portfolio_SR = (portfolio_return - Ptf.RiskFreeRate) / portfolio_std;

% Build a struct for the output
Output_struct = struct('Volatility', portfolio_std, 'Weights', weights,...
    'Return', portfolio_return, 'Sharpe_Ratio', portfolio_SR);
Output_struct.Name = name_ptf;
Output_struct.Ptf = Ptf;

% Display Portfolio - Maximum Diversified Portfolio
print_portfolio(weights, Ptf.AssetList, portfolio_return, portfolio_std, portfolio_SR, name_ptf)

end