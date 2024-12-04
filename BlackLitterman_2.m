function [Output_struct, Ptf] = BlackLitterman_2(Ptf, returns, caps, lambda, views, name_ptf, flag)

    % Compute the portfolio frontier, under standard constraints, using the
    % Black-Litterman model
    % INPUTS
    % Ptf: Portfolio object
    % returns: matrix of linear returns
    % caps: vector of capitalizations
    % lambda: risk aversion coefficient
    % views: vector of market-views
    % name_ptf: name of the portfolio
    % flag: 0 for minimum risk portfolio, 1 for maximum Sharpe ratio portfolio
    %
    % OUTPUTS
    % Output_struct: a struct containing the volatility, weights, return,
    %                and Sharpe ratio of the minimum variance portfolio
    % Ptf: Portfolio object with the moments set to the Black-Litterman moments

    % Build the views
    % Number of views
    v = length(views);
    tau = 1/length(returns);
    % Initialize Matrixes
    P = zeros(v, Ptf.NumAssets);
    q = zeros(v, 1);
    Omega = zeros(v);

    for i = 1:v
        P(i, Ptf.AssetList == views(i).overperformer) = 1;
        P(i, Ptf.AssetList == views(i).underperformer) = -1;
        q(i) = views(i).delta;
        Omega(i,i) = tau.*P(i,:)*Ptf.AssetCovar*P(i,:)';
    end

    % Change to daily returns
    bizyear2bizday = 1/252;
    q = q*bizyear2bizday;
    Omega = Omega*bizyear2bizday;

    % Capitalizations Weighted PTF - Market Point of View
    weightsCaps = (caps/sum(caps)); 
    weightsCaps = weightsCaps';

    % Market Moments
    mu_market = lambda.*Ptf.AssetCovar*weightsCaps;
    cov_market = tau.*Ptf.AssetCovar;

    % Black Litterman Moments
    % muBL = inv(inv(cov_market)+P'*inv(Omega)*P)*...
    %     (P'*inv(Omega)*q + inv(cov_market)*mu_market); 
    muBL = (cov_market \ eye(size(cov_market)) + P' * (Omega \ P)) \ ...
       (P' * (Omega \ q) + (cov_market \ mu_market));
    % covBL = inv(P'*inv(Omega)*P + inv(cov_market));
    covBL = (P' * (Omega \ P) + cov_market \ eye(size(cov_market))) \ eye(size(P, 2));


    % Set the moments in the Portfolio object
    Ptf = setAssetMoments(Ptf, muBL, Ptf.AssetCovar + covBL);

    % Estimate Frontier
    pwBL = estimateFrontier(Ptf, 100);
    [risksBL, retBL] = estimatePortMoments(Ptf, pwBL);

    if flag == 0

        % Find the minimum risk portfolio (Minimum Variance Portfolio - MVP)
        % idxMVP = find(risksBL == min(risksBL));
        [stdMVP, idxMVP] = min(risksBL);
        wMVP = pwBL(:,idxMVP);
        retMVP = retBL(idxMVP);
        % stdMVP = risksBL(idxMVP);
        srMVP = (retMVP - Ptf.RiskFreeRate)/stdMVP;

        % Display the minimum risk portfolio
        print_portfolio(wMVP, Ptf.AssetList, retMVP, stdMVP, srMVP, name_ptf);

        % Build a struct for the output
        Output_struct = struct('Volatility', stdMVP, 'Weights', wMVP,...
            'Return', retMVP, 'Sharpe_Ratio', srMVP);
        Output_struct.Name = name_ptf;
        Output_struct.Ptf = Ptf;

    elseif flag == 1

        % Find the maximum Sharpe ratio portfolio (Maximum Sharpe Ratio Portfolio - MSR)
        wMSR = estimateMaxSharpeRatio(Ptf);
        [stdMSR, retMSR] = estimatePortMoments(Ptf, wMSR);
        srMSR = (retMSR - Ptf.RiskFreeRate)/stdMSR;

        % Display the maximum Sharpe ratio portfolio
        print_portfolio(wMSR, Ptf.AssetList, retMSR, stdMSR, srMSR, name_ptf);

        % Build a struct for the output
        Output_struct = struct('Volatility', stdMSR, 'Weights', wMSR,...
            'Return', retMSR, 'Sharpe_Ratio', srMSR);
        Output_struct.Name = name_ptf;
        Output_struct.Ptf = Ptf;

    else
        % Error message, not valid flag
        error('Not valid flag, please use 0 for minimum risk portfolio and 1 for maximum Sharpe ratio portfolio')
    end


end