function [Ptf,minRisk_P1, minRiskWgt_P1, minRiskRet_P1, minRiskSR_P1] = ...
    minRiskPortfolio(Ptf, risk_free_rate , const, name_ptf)

    %   Compute the Minimum Variance Portfolio of the frontier.
    % INPUTS
    % names: cell array of asset names
    % mean_returns: vector of mean returns
    % cov_matrix: covariance matrix
    % risk_free_rate: risk free rate
    % const: constraints struct 


    % estimate efficient frontier, via object method and 100 P1oints
    pwgt1 = estimateFrontier(Ptf, 100);
    % estimate portfolio moments
    [pf_risk_P1, pf_Retn_P1] = estimatePortMoments(Ptf, pwgt1);
    % plot efficient frontier
    % plot(pf_risk_P1, pf_Retn_P1)

    % find minimum risk portfolio (Portfolio A or C - Minimum Variance Portfolio)
    minRisk_P1 = min(pf_risk_P1);
    minRiskWgt_P1 = pwgt1(:, pf_risk_P1 == minRisk_P1);
    minRiskRet_P1 = mean_returns' * minRiskWgt_P1;
    minRiskSR_P1 = (minRiskRet_P1 - risk_free_rate) / minRisk_P1;

    print_portfolio(minRiskWgt_P1, names, minRiskRet_P1, minRisk_P1, minRiskSR_P1, name_ptf);


    % if flag == 0
    %     print_portfolio(minRiskWgt_P1, names, minRiskRet_P1, minRisk_P1, minRiskSR_P1,'Minimum Risk Portfolio (A)')
    % elseif flag == 1
    %     print_portfolio(minRiskWgt_P1, names, minRiskRet_P1, minRisk_P1, minRiskSR_P1,'Minimum risk Portfolio with constrints (C)')
    % end

end