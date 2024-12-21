function [] = plotData(returns, prices, names)
% This function plots the returns and prices of the assets in the dataset,
% it plots the QQ-plot, the correlation plot and the histograms of the
% returns, and performs the Shapiro-Wilk/Shapiro-Francia and the
% Kolmogorov-Smirnov Normality tests.
%
% INPUTS:
% returns:  Log returns of the assets
% prices:   Prices of the assets
% names:    Names of the assets

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
%% Shapiro-Wilk/Shapiro-Francia and Kolmogorov-Smirnov (K-S Test) normality tests
% Variables declarations
H_SW = zeros(size(returns,2), 1); pValue_SW = zeros(size(returns,2), 1); W_SW = zeros(size(returns,2), 1);
H_KS = zeros(size(returns,2), 1); pValue_KS = zeros(size(returns,2), 1);

disp('====================================================================================================')
fprintf('SW-test and K-S test for normality\n')
disp('====================================================================================================')
fprintf('Legend: \n- name --> H = 0: accept normality \n- name --> H = 1: reject normality \n\n')

for i = 1 : size(returns,2)
    % swtest performs the Shapiro-Francia test when the series is Leptokurtik (kurtosis > 3), 
    % otherwise it performs the Shapiro-Wilk test.
    [H_SW(i), pValue_SW(i), W_SW(i)] = swtest(returns(:,i));

    % Kolmogorov-Smirnov Test (K-S Test) for normality:
    test_cdf = makedist('Normal','mu',mean(returns(:,i)),'sigma',std(returns(:,i)));
    [H_KS(i), pValue_KS(i)] = kstest(returns(:,i), 'CDF', test_cdf);
end

% Create Table
tests_resultTable = table(string(names'), H_SW, pValue_SW, H_KS, pValue_KS, ...
    'VariableNames', {'Asset_Name', 'H_SW', 'pValue_SW', 'H_KS', 'pValue_KS'});

% Display Table
disp(tests_resultTable);
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

%% Corr plot
V = corr(returns);
figure;
heatmap(V);
ax = gca; % Get current axes
ax.XDisplayLabels = names;
ax.YDisplayLabels = names;
colormap summer;
end