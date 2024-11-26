function [weights_n, portfolio_n_return, portfolio_n_std, portfolio_n_SR] = ...
          Max_Entropy_Portfolio(mean_returns, cov_matrix, capitalizations, names, risk_free_rate, mkt, P2,num_assets,const)

        % Compute the Maximum Diversified Portfolio (Portfolio M) and the
        % Maximum Entropy (in asset volatility) Portfolio (Portfolio N), under
        % the following constraints (to be considered all at once):
        % - Standard constraints,
        % - The total exposure on cyclicals has to be greater than 20%,
        % - Assuming that you have a benchmark portfolio (capitalization
        %   weighted portfolio), the sum of the difference (in absolute value)
        %   of the weights in the benchmark portfolio and the optimal weights
        %   has to be greater than 20%

        % INPUTS
        % mean_returns: vector of mean returns of the assets
        % cov_matrix: covariance matrix of the assets
        % capitalizations: cell array with the capitalizations of the assets
        % names: cell array with the names of the assets
        % risk_free_rate: risk free rate
        % mkt: structure with the market information
        % P2: structure with the information about the sectors
        % num_assets: number of assets
        % const: structure with the constraints



        % Set up optimization problem
        % Initial guess
        initial_guess = ones(num_assets, 1) / num_assets;
        % Options
        options = optimoptions('fmincon', ...     
                'Algorithm', 'sqp', ... % Specify the algorithm     
                'StepTolerance', 1e-6, ...       % Smaller than default StepTolerance     
                'Display', 'off');               % Show iteration information

        % Portfolio N:  Maximum Entropy (in asset volatility) Portfolio
        % entropy = @(w) sum(w.^2' * diag(cov_matrix)/ sum(w.^2' * diag(cov_matrix)) * log(w.^2' * diag(cov_matrix)/ sum(w.^2' * diag(cov_matrix))));
        entropy = @(w) sum(w.^2 .* diag(cov_matrix)/ sum(w.^2 .* diag(cov_matrix)) .* log(w.^2 .* diag(cov_matrix)/ sum(w.^2 .* diag(cov_matrix))));%NEW Matte

        [weights_n, minvalue_n] = fmincon(entropy, initial_guess, const.A, const.b, const.Aeq, const.beq, const.lb, const.ub, const.nonlinconstr, options);
        portfolio_n_return = mean_returns' * weights_n;

        portfolio_n_std = sqrt(weights_n' * cov_matrix * weights_n); %NEW

        portfolio_n_SR = (portfolio_n_return - risk_free_rate) / portfolio_n_std;

        % Display Portfolio N - Maximum Entropy Portfolio

        print_portfolio(weights_n, names, portfolio_n_return, portfolio_n_std, portfolio_n_SR,'N')
end