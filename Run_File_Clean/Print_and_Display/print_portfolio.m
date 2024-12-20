function [] =  print_portfolio(weights,names,returns,risk,SR,flag) 
    
    % Print the portfolio weights, expected return, volatility, and Sharpe ratio
    % for the given portfolio.
    %
    % INPUTS:
    % weights: a vector of portfolio weights
    % names: a cell array of asset names
    % returns: the expected return of the portfolio
    % risk: the volatility of the portfolio
    % SR: the Sharpe ratio of the portfolio
    % flag: a string to identify the portfolio  

    disp('=======================================================')
    fprintf('%s\n', flag)
    disp('=======================================================')

    disp('Asset Name                Weight')
    disp('-------------------------------------------')
    for i = 1:length(weights)
        if round(weights(i),5) == 0
            fprintf('%-25s %s\n', names{i}, '-');
        else
        fprintf('%-25s %.4f\n', names{i}, round(weights(i),14));
        end
    end
    disp('-------------------------------------------')
    fprintf('%-25s %.4f\n', 'Expected Return', returns);
    fprintf('%-25s %.4f\n', 'Volatility', risk);
    fprintf('%-25s %.4f\n', 'Sharpe Ratio', SR);
    fprintf('%-25s %.4f\n', 'Sum of weights', sum(weights));
    disp('  ')

end