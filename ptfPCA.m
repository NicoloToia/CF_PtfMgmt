function [ptf] = ptfPCA(ret, tgtExplained, tgtRisk)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Standardize returns
mean_ret = mean(ret);
std_ret = (ret - mean_ret) ./ std(ret);
k = length(mean_ret);
[~,~,latent,~,~,~] = pca(std_ret, 'NumComponents',k);
% Find components by target Expaineed Variance
TotVar = sum(latent);
explainedVar = latent(1:k)/TotVar;
CumExplainedVar = cumsum(explainedVar);
find(CumExplainedVar >= tgtExplained, 1)
k_choice = 6;
[factorLoading,factorRetn] =...
    pca(std_ret, 'NumComponents',k_choice);
% Reconstruct 
covarFactor = cov(factorRetn);
reconReturn = factorRetn * factorLoading' + mean_ret;
unexplainedRetn = ret - reconReturn; % IDIOSYNCRATHIC PART (epsilon)
unexplainedCov = diag(cov(unexplainedRetn));
D = diag(unexplainedCov);
% Using PTF Object -> find portfolio by target risk
covarAsset = factorLoading*covarFactor*factorLoading'+D;
port = Portfolio("AssetMean", mean_ret, 'AssetCovar', round(covarAsset,13),...
    'LowerBound', 0, 'UpperBound', 1, ...
    'Budget', 1);
weights = estimateFrontierByRisk(port, tgtRisk^2);

ptf = struct();
ptf.w = w;
ptf.ret = mean_ret * weights;
ptf.std = weights' * cov(ret) * weights;
ptf.sr = (ptf.ret - rf) / ptf.std;

end

