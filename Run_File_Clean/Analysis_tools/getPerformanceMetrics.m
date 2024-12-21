function [annRet, annVol, Sharpe, MaxDD, Calmar] = getPerformanceMetrics(x)
% This function computes the annualized metrics
%
% INPUTS:
% x:        Equity Vector
%
% OUTPUTS:
% annRet:   Annualized Return
% annVol:   Annualized Volatility
% Sharpe:   Annualized Sharpe Ratio
% MaxDD:    Maximum Drawdown
% Calmar:   Calmar Ratio

% Annualized return
annRet = (x(end) / x(1)) .^ (1 / (length(x) / 252)) - 1;
   
% Annualized volatility
annVol = std(tick2ret(x)) * sqrt(252);
    
% Sharpe ratio
Sharpe = annRet / annVol; 
    
% Maximum drawdown
dd = zeros(1, length(x));
    
for i = 1 : length(x)
    dd(i) = (x(i) / max(x(1 : i))) - 1;
end
    
MaxDD = min(dd); 
    
% Calmar ratio
Calmar = annRet / abs(MaxDD); 

end

