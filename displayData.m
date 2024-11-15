clc
clear all
close all
warning off

%% PART A 

% Load data from Excel files
prices = readtable('prices.xlsx');
capitalizations = readtable('capitalizations.xlsx');

% Names
names = capitalizations.Properties.VariableNames(2:end);

% Extract dates and data
dates = prices{:,1}; % First column contains dates
prices_data = prices{:,2:end}; % Data starts from the second column

% Convert dates to MATLAB date format if needed
dates = datetime(dates);

% Filter prices for 2023
start_date = datetime(2023,1,1);
end_date = datetime(2023,12,31);
prices_2023 = prices_data(dates >= start_date & dates <= end_date, :);

% Calculate daily returns for each index in 2023
returns_2023 = diff(log(prices_2023));

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
for i =1:16
    V(i,i) = 0;
end
figure;
heatmap(V);
colormap summer;
%%




