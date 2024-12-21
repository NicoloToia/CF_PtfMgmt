function [tableRC] = getRiskContributionsTables(w,prices)
% Function to display in a table the relative risk contributions
% for each asset in the portfolio
%
% INPUTS:
% w:        MATLAB Table with weights of different portfolios
% prices:   Historical series of prices
%
% OUTPUT:
% tableRC:  MATLAB Table with relative risk contributions for different
%           portfolios

% Fetch names of ptfs and assets
ptfNames = w.Properties.VariableNames;
assetNames = w.Properties.RowNames;

% Weigths and dimensions
w = table2array(w);
numPtfs = size(w,2);
numAssets = size(w,1);

% Compute the returns
ret = tick2ret(prices);

% Preallocation of memory
relRCs = zeros(numAssets, numPtfs);

% Computation of the relative risk contributions
for i = 1:numPtfs
    [relRC, ~, ~] = getRiskContributions(w(:,i), ret);
    relRCs(:,i) = relRC;
end

% Creation of the output MATLAB Table
tableRC = array2table(relRCs,...
    "RowNames",assetNames,...
    "VariableNames",ptfNames);

end

