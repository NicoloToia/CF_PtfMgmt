function Output_structure = maxSharpPortfolio(Ptf,risk_free_rate, name_Ptf)

    %   Compute the efficient frontier under the standard constraints.
    %   Compute the Maximum Sharpe Ratio Portfolio
    %
    % INPUTS
    % Ptf: Portfolio data
    % risk_free_rate: risk free rate
    % name_Ptf: name of the portfolio

    % find maximum Sharpe ratio portfolio (Portfolio B or D - Maximum Sharpe Ratio Portfolio)
    maxSharpeWgt_Ptf = estimateMaxSharpeRatio(Ptf);
    maxSharpeRet_Ptf = Ptf.AssetMean' * maxSharpeWgt_Ptf;
    maxSharpeRisk_Ptf = sqrt(maxSharpeWgt_Ptf' * Ptf.AssetCovar * maxSharpeWgt_Ptf);
    maxSharpeSR_Ptf = (maxSharpeRet_Ptf - risk_free_rate) / maxSharpeRisk_Ptf;

    print_portfolio(maxSharpeWgt_Ptf, Ptf.AssetList, maxSharpeRet_Ptf, maxSharpeRisk_Ptf, maxSharpeSR_Ptf, name_Ptf)

    Output_structure = struct('Volatility', maxSharpeRisk_Ptf, 'Weights', maxSharpeWgt_Ptf, 'Return', maxSharpeRet_Ptf, 'Sharpe_Ratio', maxSharpeSR_Ptf );
    Output_structure.name = name_Ptf;
    Output_structure.Ptf = Ptf;
end