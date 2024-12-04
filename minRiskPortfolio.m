function Output_struct = minRiskPortfolio(Ptf, risk_free_rate, pwgt, pf_risk_Ptf, name_ptf)

    % Compute the Minimum Variance Portfolio of the frontier.
    % INPUTS
    % Ptf: Portfolio object
    % risk_free_rate: risk free rate
    % pwgt: the weights of the portfolios on the frontier
    % pf_risk_Ptf: the volatility of the portfolios on the frontier
    % name_ptf: name of the portfolio
    %
    % OUTPUTS
    % Output_struct: a struct containing the volatility, weights, return, 
    %                and Sharpe ratio of the minimum variance portfolio

    % find minimum risk portfolio (Minimum Variance Portfolio)
    minRisk_Ptf = min(pf_risk_Ptf);
    minRiskWgt_Ptf = pwgt(:, pf_risk_Ptf == minRisk_Ptf);
    minRiskRet_Ptf = Ptf.AssetMean' * minRiskWgt_Ptf;
    minRiskSR_Ptf = (minRiskRet_Ptf - risk_free_rate) / minRisk_Ptf;

    % Display the minimum risk portfolio
    print_portfolio(minRiskWgt_Ptf, Ptf.AssetList, minRiskRet_Ptf, minRisk_Ptf, minRiskSR_Ptf, name_ptf);

    % Build a struct for the output
    Output_struct = struct('Volatility', minRisk_Ptf, 'Weights', minRiskWgt_Ptf,...
        'Return', minRiskRet_Ptf, 'Sharpe_Ratio', minRiskSR_Ptf);
    Output_struct.Name = name_ptf;
    Output_struct.Ptf = Ptf;

end