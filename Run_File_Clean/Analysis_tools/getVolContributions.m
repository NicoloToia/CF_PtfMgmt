function vols = getVolContributions(x, LogRet)
% This function computes the volatility contributions of a portfolio

    vols = ((x.^2).*(std(LogRet)'.^2))/sum((x.^2).*(std(LogRet)'.^2));
end