function [relRC, RC, mVol] = getRiskContributions(x, Ret)
    V = cov(Ret);
    VolaPtf = sqrt(x'*V*x);
    mVol = V*x/VolaPtf;
    RC = mVol.*x;
    relRC = RC/sum(RC);
end