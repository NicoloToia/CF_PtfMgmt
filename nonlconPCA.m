function [c, ceq] = nonlconPCA(x, fl, D, sigmaF, tRisk)
    % Calculate sigma of ptf
    sigmaP  = sqrt((fl'*x)'*sigmaF*(fl'*x)+x'*D*x);
    % Constraint on volatility
    c = sigmaP - tRisk;  % Inequality constraint (c <= 0)
    ceq = [];            % No equality constraint
end

