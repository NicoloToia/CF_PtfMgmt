function [] = disp_weights(table)
%DISP_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here
ws = table2array(table);
names = table.Properties.VariableNames;
asset_names = table.Properties.RowNames;
n_assets = length(asset_names);
n_ptfs = size(ws,2);
figure;
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
cmap = parula(length(n_assets)); % Generate colormap with colors for each asset
b = bar(ws', 'stacked'); % Create stacked bar plot
b
% Assign consistent colors for each asset
for j = 1:n_assets
    b(j).FaceColor = 'flat'; % Enable individual bar coloring
    b(j).CData = repmat(cmap(j), size(ws, 2), 1); % Assign consistent color to the asset
end

legend(asset_names, 'Location', 'bestoutside'); % Add legend for the assets
xlabel('Portfolio Index');
ylabel('Weight');
title('Portfolio Weights Distribution');

end

