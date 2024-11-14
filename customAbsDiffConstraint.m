function [c, ceq] = customAbsDiffConstraint(weights, benchmarkWeights)
    % Calculate the absolute difference in weights
    absDiff = abs(weights' - benchmarkWeights);
    
    % Sum of absolute differences should be at least 0.2
    c = 0.2 - sum(absDiff);  % Inequality constraint (c <= 0)
    ceq = [];                % No equality constraint
end