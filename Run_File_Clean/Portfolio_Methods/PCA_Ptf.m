function Output_struct = PCA_Ptf(Ptf, ret, name_ptf)
% This function performs a PCA analysis on the returns of the assets in the
% portfolio and then constructs a minimum risk portfolio using the factors
% obtained from the PCA analysis.
%
% INPUTS:
% Ptf:      a structure containing the portfolio information
% ret:      a matrix of asset returns
% name_ptf: a string to identify the portfolio
%
% OUTPUT:
% Output_struct: a structure containing the portfolio weights, volatility,
%                return, and Sharpe ratio

% Get the number of samples
samples = size(ret, 1);

mean_ret = mean(ret);
sd_ret = std(ret);
std_ret = (ret - repmat(mean_ret,samples,1)) ./ repmat(sd_ret,samples,1);

% Get the number of assets in the portfolio
k = Ptf.NumAssets;
% Perform PCA analysis
[~,~,latent,~,~,~] = pca(std_ret, 'NumComponents',k);
% coeff -> loadings
% score -> scores matrix
% latent
% explained -> Percentage of total variance explained column vector

% Find cumulative explained variance
TotVar = sum(latent);
explainedVar = latent./TotVar;
CumExplainedVar = cumsum(explainedVar);
k_hat = find(CumExplainedVar >= .9, 1);

% Call the plot function for the pca analysis
figure();
cmap = parula(length(explainedVar)); 

subplot(1,2,1)
b = bar(explainedVar,'FaceColor','white','EdgeColor','k');
for i = 1:length(explainedVar)
    b.FaceColor = 'flat'; % Enable individual bar coloring
    b.CData(i, :) = cmap(i, :); % Assign color to each bar
end
title("% Explained Variance");

subplot(1,2,2)
b = bar(CumExplainedVar,'FaceColor','white','EdgeColor','k');
for i = 1:length(explainedVar)
    b.FaceColor = 'flat'; % Enable individual bar coloring
    b.CData(i, :) = cmap(i, :); % Assign color to each bar
end
title("% Cumulative Variance");
% set the line at 90% cumulative variance
yline(0.90, '--k', 'LineWidth',1.5);

%% PTF OPTIMIZATION

% call the pca function with the new value of k
[factorLoading,factorRetn,~,~,~,~] = pca(std_ret, 'NumComponents',k_hat);

% calculate the covariance matrix of the factors
covarFactor = cov(factorRetn);  

% Reconstruction of the return
reconReturn = (factorRetn * factorLoading') .* repmat(sd_ret,samples,1) + ...
                repmat(mean_ret,samples,1);

% calculate the IDIOSYNCRATHIC PART (epsilon)
unexplainedRetn = ret - reconReturn; 

% calculate the covariance matrix of the idiosyncratic part
unexplainedCov = diag(cov(unexplainedRetn));
D = diag(unexplainedCov);
covarAsset = factorLoading*covarFactor*factorLoading' + D;

% set the portfolio moments
Ptf.AssetMean = mean_ret;
Ptf.AssetCovar = round(covarAsset,10);

weights = estimateFrontierByRisk(Ptf, 0.7);

% compute the volatility of the portfolio
volatility = sqrt(weights' * covarAsset * weights);

% compute the return of the portfolio
returns = mean_ret * weights;

% compute the Sharpe ratio of the portfolio
sharpe_ratio = (returns - Ptf.RiskFreeRate) / volatility;

% Display the minimum risk portfolio
print_portfolio(weights, Ptf.AssetList, returns, volatility, sharpe_ratio, name_ptf);

 % display the value of k that explains 90% of the variance (k = 6)
disp(['The number of factors that explains more than 90% of the cumulative variance is: ', num2str(k_hat)]);

Output_struct = struct('Weights', weights, 'Volatility', volatility,...
                       'Return', returns, 'Sharpe_Ratio', sharpe_ratio);
Output_struct.Name = name_ptf;
Output_struct.Ptf = Ptf;

end