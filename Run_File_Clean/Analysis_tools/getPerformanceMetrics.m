function [annRet, annVol, Sharpe, MaxDD, Calmar] = getPerformanceMetrics(x)

    annRet = (x(end) / x(1)) .^ (1 / (length(x) / 252)) - 1;
    
    annVol = std(tick2ret(x)) * sqrt(252);
    
    Sharpe = annRet / annVol; 
    
    dd = zeros(1, length(x));
    
    for i = 1 : length(x)
        dd(i) = (x(i) / max(x(1 : i))) - 1;
    end
    
    MaxDD = min(dd); 
    
    Calmar = annRet / abs(MaxDD); 

end

