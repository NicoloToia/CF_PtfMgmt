function [] = plotData(returns_2023, prices_2023, names)

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
    yy = returns_2023(:,i);
    plot(yy, 'Color', colors(groupIdx), 'LineWidth', 1.25);
    title(names(i));
end
% Prices
figure;
for i = 1:16
    subplot(4,4,i);
    groupIdx = find(cellfun(@(v) ismember(names(i), v), groups));
    colors = ["#77AC30", "#EDB120", "#A2142F", "#4DBEEE"];
    yy = prices_2023(:,i);
    plot(yy, 'Color', colors(groupIdx), 'LineWidth', 1.25);
    title(names(i));
end
%% QQ plot
figure;
for i = 1:16
    subplot(4,4,i);
    groupIdx = find(cellfun(@(v) ismember(names(i), v), groups));
    colors = ["#77AC30", "#EDB120", "#A2142F", "#4DBEEE"];
    yy = returns_2023(:,i);
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
for i = 1 : size(returns_2023,2)
    % swtest performs the Shapiro-Francia test when the series is Leptokurtik (kurtosis > 3), 
    % otherwise it performs the Shapiro-Wilk test.
    [H(i), pValue(i), W(i)] = swtest(returns_2023(:,i));
    fprintf('%s --> H%d\n', names{i}, H(i));
end
fprintf('\n')
%% Histograms
figure;
for i = 1:16
    subplot(4,4,i);
    groupIdx = find(cellfun(@(v) ismember(names(i), v), groups));
    colors = ["#77AC30", "#EDB120", "#A2142F", "#4DBEEE"];
    yy = returns_2023(:,i);
    histogram(yy, 'FaceColor',colors(groupIdx));
    title(names(i));
end

%% Scatter
% for i = 1:16
%     for j = 1:16
%     subplot(16,16, (i-1) * 16 + j);
%     xx = returns_2023(:,i);
%     yy = returns_2023(:,j);
%     scatter(xx, yy)
%     end
% end
%% Corr plot
V = corr(returns_2023);
figure;
heatmap(V);
% Set custom x-axis and y-axis labels
ax = gca; % Get current axes
ax.XDisplayLabels = names; % Set x-axis labels
ax.YDisplayLabels = names; % Set y-axis labels
colormap summer;