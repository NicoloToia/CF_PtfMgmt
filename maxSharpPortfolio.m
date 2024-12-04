function Output_struct = maxSharpPortfolio(Ptf,risk_free_rate, name_ptf)

    % Compute the Maximum Sharpe Ratio Portfolio
    % INPUTS
    % Ptf: Portfolio object
    % risk_free_rate: risk free rate
    % name_ptf: name of the portfolio
    %
    % OUTPUTS
    % Output_struct: a struct containing the volatility, weights, return,
    %                   and Sharpe ratio of the maximum Sharpe ratio portfolio

    % find maximum Sharpe ratio portfolio (Maximum Sharpe Ratio Portfolio)
    maxSharpeWgt_Ptf = estimateMaxSharpeRatio(Ptf);
    maxSharpeRet_Ptf = Ptf.AssetMean' * maxSharpeWgt_Ptf;
    maxSharpeRisk_Ptf = sqrt(maxSharpeWgt_Ptf' * Ptf.AssetCovar * maxSharpeWgt_Ptf);
    maxSharpeSR_Ptf = (maxSharpeRet_Ptf - risk_free_rate) / maxSharpeRisk_Ptf;

    % Display the maximum Sharpe ratio portfolio
    print_portfolio(maxSharpeWgt_Ptf, Ptf.AssetList, maxSharpeRet_Ptf, maxSharpeRisk_Ptf, maxSharpeSR_Ptf, name_ptf)

    % Build a struct for the output
    Output_struct = struct('Volatility', maxSharpeRisk_Ptf, 'Weights', maxSharpeWgt_Ptf,...
        'Return', maxSharpeRet_Ptf, 'Sharpe_Ratio', maxSharpeSR_Ptf );
    Output_struct.Name = name_ptf;
    Output_struct.Ptf = Ptf;
    
end