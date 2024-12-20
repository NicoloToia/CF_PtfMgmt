function [relRC, RC, mVol] = getRiskContributions(x, Ret)
% This function computes the relative risk contributions of a portfolio

    V = cov(Ret);
    VolaPtf = sqrt(x'*V*x);
    mVol = V*x/VolaPtf;
    RC = mVol.*x;
    relRC = RC/sum(RC);
end