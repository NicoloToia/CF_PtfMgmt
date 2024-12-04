function [] = plotData(returns, prices, names)
% This function plots the returns and prices of the assets in the dataset
% and performs the Shapiro-Wilk and Shapiro-Francia normality tests.
%
% INPUTS:
% returns: returns of the assets
% prices: prices of the assets
% names: names of the assets

cyclical = ["ConsumerDiscretionary", "Financials", "Materials", "RealEstate", "Industrials"];
defensive = ["ConsumerStaples", "Utilities", "HealthCare"];
sensible = ["Energy", "InformationTechnology", "CommunicationServices"];

factor = ["Momentum","Value","Growth","Quality","LowVolatility"];

groups = {cyclical, defensive, sensible, factor};

%% Plot
% Returns
for i = 1:16
    subplot(4,4,i);
    groupIdx = find(cellfun(@(v) ismember(names(i), v), groups));
    colors = ["#77AC30", "#EDB120", "#A2142F", "#4DBEEE"];
    yy = returns(:,i);
    plot(yy, 'Color', colors(groupIdx), 'LineWidth', 1.25);
    title(names(i));
end
% Prices
figure;
for i = 1:16
    subplot(4,4,i);
    groupIdx = find(cellfun(@(v) ismember(names(i), v), groups));
    colors = ["#77AC30", "#EDB120", "#A2142F", "#4DBEEE"];
    yy = prices(:,i);
    plot(yy, 'Color', colors(groupIdx), 'LineWidth', 1.25);
    title(names(i));
end
%% QQ plot
figure;
for i = 1:16
    subplot(4,4,i);
    groupIdx = find(cellfun(@(v) ismember(names(i), v), groups));
    colors = ["#77AC30", "#EDB120", "#A2142F", "#4DBEEE"];
    yy = returns(:,i);
    qq = qqplot(yy);
    set(qq(1),'Marker','o','MarkerFaceColor', colors(groupIdx), ...
        'MarkerEdgeColor', 'none'); % Change points color
    title(names(i));
end
%% Shapiro-Wilk and Shapiro-Francia normality tests
fprintf('\nlegend \nname --> H0: accept normality \nname --> H1: reject normality \n\n')
H = [];
pValue = [];
W = [];
for i = 1 : size(returns,2)
    % swtest performs the Shapiro-Francia test when the series is Leptokurtik (kurtosis > 3), 
    % otherwise it performs the Shapiro-Wilk test.
    [H(i), pValue(i), W(i)] = swtest(returns(:,i)); %%% SI PUO CAMBIARE ALPHA E CAPIRE CON CHE ALPHA SI ACCETTA L'IPOTESI NULLA NORMALE
    %%% ADD THE TTEST TO COMPARE THE TWO TESTS, AND SHOW THA ìT FINANCIAL DATA HAVE HEAVY TAILS
    fprintf('%s --> H%d\n', names{i}, H(i));
end
fprintf('\n')
%% Histograms
figure;
for i = 1:16
    subplot(4,4,i);
    groupIdx = find(cellfun(@(v) ismember(names(i), v), groups));
    colors = ["#77AC30", "#EDB120", "#A2142F", "#4DBEEE"];
    yy = returns(:,i);
    histogram(yy, 'FaceColor',colors(groupIdx));
    title(names(i));
end

%% Scatter
% for i = 1:16
%     for j = 1:16
%     subplot(16,16, (i-1) * 16 + j);
%     xx = returns(:,i);
%     yy = returns(:,j);
%     scatter(xx, yy)
%     end
% end
%% Corr plot
V = corr(returns);
for i =1:16
    V(i,i) = 0;
end
figure;
heatmap(V);
colormap summer;