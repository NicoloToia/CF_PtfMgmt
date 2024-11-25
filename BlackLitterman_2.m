function [ptfMVP, ptfMSR] = ...
    BlackLitterman( prices,...
                    caps,...
                    names,...
                    views,...
                    rf)

    % Compute the portfolio frontier, under standard constraints, using the
    % Black-Litterman model with the following views (to be considered all
    % at once):
    % •View on Technology vs. Financials: Given the growing im-
    % portance of technology, you think that the Technology sector will
    % outperform the Financial sector of the 2% (annual).
    % •View on the Momentum vs. Low Volatility factor: You
    % might assume that Momentum will outperform Low Volatility in
    % a bull market scenario (annual overperformance of 1%)
    % Compute the Minimum Variance Portfolio, named Portfolio I, and the
    % Maximum Sharpe Ratio Portfolio, named Portfolio L, of the frontier.

    % INPUTS
    % prices: matrix of prices
    % caps: vector of market capitalizations
    % names: cell array of asset names
    % views: structure array with fields overperformer, underperformer, delta
    % rf: risk-free rate

    returns = tick2ret(prices);
    % Calculate mean returns and covariance matrix
    mean_returns = mean(returns)';
    cov_matrix = cov(returns);
    % Define number of assets
    num_assets = length(names);
    % Build the views
    % Number of views
    v = length(views);
    tau = 1/length(returns);
    % Initialize Matrixes
    P = zeros(v, num_assets);
    q = zeros(v, 1);
    Omega = zeros(v);
    for i = 1:v
        P(i, find(names == views(i).overperformer)) = 1;
        P(i, find(names == views(i).underperformer)) = -1;
        q(i) = views(i).delta;
        Omega(i,i) = tau.*P(i,:)*cov_matrix*P(i,:)';
    end
    % Change to daily returns
    bizyear2bizday = 1/252;
    q = q*bizyear2bizday;
    Omega = Omega*bizyear2bizday;
    % Capitalizations Weighted PTF - Market Point of View
    weightsCaps = caps/sum(caps); weightsCaps = weightsCaps';
    lambda = 1.2;
    mu_market = lambda.*cov_matrix*weightsCaps;
    cov_market = tau.*cov_matrix;
    % Black Litterman Moments
    muBL = inv(inv(cov_market)+P'*inv(Omega)*P)*...
        (P'*inv(Omega)*q + inv(cov_market)*mu_market); 
    covBL = inv(P'*inv(Omega)*P + inv(cov_market));
    % Create Portfolio for Black & Litterman
    ptf = Portfolio('NumAssets', num_assets, 'Name', 'MV with BL');
    ptf = setDefaultConstraints(ptf);
    ptf.RiskFreeRate = rf;
    ptf = setAssetMoments(ptf, muBL, cov_matrix + covBL);
    % Estimate Frontier
    pwBL = estimateFrontier(ptf, 100);
    [risksBL, retBL] = estimatePortMoments(ptf, pwBL);
    % Find MVP
    idxMVP = find(risksBL == min(risksBL));
    wMVP = pwBL(:,idxMVP);
    retMVP = retBL(idxMVP);
    stdMVP = risksBL(idxMVP);
    srMVP = (retMVP - rf)/stdMVP;
    % Max Sharpe Ratio 
    wMSR = estimateMaxSharpeRatio(ptf);
    [stdMSR, retMSR] = estimatePortMoments(ptf, wMSR);
    srMSR = (retMSR - rf)/stdMSR;

    ptfMVP = struct();
    ptfMVP.w = wMVP;
    ptfMVP.ret = retMVP;
    ptfMVP.std = stdMVP;
    ptfMVP.sr = srMVP;

    ptfMSR = struct();
    ptfMSR.w = wMSR;
    ptfMSR.ret = retMSR;
    ptfMSR.std = stdMSR;
    ptfMSR.sr = srMSR;

end

