function DR = getDiversificationRatio(x, Ret)
% This function computes the diversification ratio of a portfolio

    vola = std(Ret);
    V = cov(Ret);
    volaPtf =sqrt(x'*V*x);
    DR = (x'*vola')/volaPtf;
end