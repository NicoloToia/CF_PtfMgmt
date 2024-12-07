function [tableRC] = getRiskContributionsTables(w,prices)
%GETRISKCONTRIBUTIONSTABLES: Displays in a table the rel risk contribution
%for each asset in the portfolio

% Fetch names of ptfs and assets
ptfNames = w.Properties.VariableNames;
assetNames = w.Properties.RowNames;
% Weigths and dimensions
w = table2array(w);
numPtfs = size(w,2);
numAssets = size(w,1);
ret = tick2ret(prices);
relRCs = zeros(numAssets, numPtfs);
for i = 1:numPtfs
[relRC, ~, ~] = getRiskContributions(w(:,i), ret);
relRCs(:,i) = relRC;
end
tableRC = array2table(relRCs,...
    "RowNames",assetNames,...
    "VariableNames",ptfNames);
end

