function [c, ceq] = customAbsDiffConstraint(weights, benchmarkWeights)
% Function to impose the non linear constraint that the sum of the
% difference (in absolute value) of the weights in the benchmark
% portfolio and the optimal weights has to be greater than 20%
% 
% INPUTS
% weights:          Weights of the portfolio
% benchmarkWeights: Weights of the benchmark portfolio
%
% OUTPUTS
% c:                Inequality non-linear constraint
% ceq               Equality non-linear constraint

% Calculate the absolute difference in weights
absDiff = abs(weights - benchmarkWeights);

% Sum of absolute differences should be at least 0.2
c = 0.2 - sum(absDiff);  % Inequality constraint (c <= 0)
ceq = [];                % No equality constraint

end