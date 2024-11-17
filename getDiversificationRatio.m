function DR = getDiversificationRatio(x, Ret)
    vola = std(Ret);
    V = cov(Ret);
    volaPtf =sqrt(x'*V*x);
    DR = (x'*vola')/volaPtf;
end