function [] = disp_weights(table)
%DISP_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here
hexColors = {
    '#AF1740', '#1F4529', '#FFC436', '#526E48', '#DE7C7D', ...
    '#C2FFC7', '#F7E987', '#CC2B52', '#F8DE22', '#9EDF9C', ...
    '#62825D', '#133E87', '#4A628A', '#7AB2D3', '#B9E5E8', ...
    '#DFF2EB'};
colors = hexToRGB(hexColors);

ws = table2array(table);
names = table.Properties.VariableNames;
asset_names = table.Properties.RowNames;
n_assets = length(asset_names);
n_ptfs = size(ws,2);
pie_weights_fig = figure;
set(pie_weights_fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
for i=1:n_ptfs
    subplot(2,round(n_ptfs/2),i)
    labels = [];
    for j=1:n_assets
        if ws(j,i)< 4/100
            labels = [labels " "];
        else
            labels = [labels {sprintf('%.1f%%', ws(j,i)*100)}];
        end
    end
    pie(ws(:,i),labels);
    colormap(colors)
    title(names(i))
end
leg = legend(asset_names);
axis off;
set(leg,...
    'Position',[0.36092470716097 ...
                0.484668335419274...
                0.27580998671658 ...
                0.119524405506884],...
    'NumColumns',2);
sgtitle("Portfolio Weights");

% Alternative
figure;
cmap = colors; % Generate colormap with colors for each asset
b = bar(ws', 'stacked'); % Create stacked bar plot
b
% Assign consistent colors for each asset
for j = 1:n_assets
    b(j).FaceColor = 'flat'; % Enable individual bar coloring
    b(j).CData = repmat(colors(j,:),size(n_ptfs,1),1); % Assign consistent color to the asset
end

legend(asset_names, 'Location', 'bestoutside'); % Add legend for the assets
xlabel('Portfolio Index');
xticklabels(names)
ylabel('Weight');
title('Portfolio Weights Distribution');

end

