function [] = disp_weights(table, desiredOrder, market)

    % This function displays the weights of the portfolios in a pie chart
    % and a stacked bar plot. The pie chart shows the weights of the assets
    % in each portfolio, while the stacked bar plot groups the assets into
    % sectors and factors, and shows the weights of the assets in each
    % portfolio.
    % INPUTS
    % table:                table of portfolio weights
    % desiredOrder:         desired order of the assets
    % market:               market data

    % Define colors for the pie chart and the stacked bar plot
    hexColors = {
        '#AF1740', '#1F4529', '#FFC436', '#526E48', '#DE7C7D', ...
        '#C2FFC7', '#F7E987', '#CC2B52', '#F8DE22', '#9EDF9C', ...
        '#62825D', '#133E87', '#4A628A', '#7AB2D3', '#B9E5E8', ...
        '#DFF2EB'};
    colors = hexToRGB(hexColors);

    %% Pie chart
    ws = table2array(table);
    names = table.Properties.VariableNames;
    asset_names = table.Properties.RowNames;
    n_assets = length(asset_names);
    n_ptfs = size(ws,2);
    pie_weights_fig = figure;
    set(pie_weights_fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    [~, order] = ismember(desiredOrder, asset_names);
    ws = ws(order, :); % Reorder rows of weights matrix
    colors = colors(order, :); % Reorder colors accordingly
    asset_names = asset_names(order); % Reorder asset names

    % % Define groups for sectors and factors
    [group_names, group_members] = build_groups(market);

    group_colors = [
        mean(colors(ismember(asset_names, group_members{1}), :), 1); % Cyclical
        mean(colors(ismember(asset_names, group_members{2}), :), 1); % Defensive
        mean(colors(ismember(asset_names, group_members{3}), :), 1); % Sensitive
        mean(colors(ismember(asset_names, group_members{4}), :), 1)  % Factors
    ];
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

    %% Stacked bar plot
    % Create stacked bar plot
    figure;
    b = bar(ws', 'stacked'); % Create stacked bar plot

    % Assign consistent colors for each asset
    for j = 1:n_assets
        b(j).FaceColor = 'flat'; % Enable individual bar coloring
        b(j).CData = repmat(colors(j,:), size(ws, 2), 1); % Assign consistent color to the asset
    end

    % Primary legend for individual assets
    legend(asset_names, 'Location', 'bestoutside'); % Add legend for the assets

    % Secondary legend clustering the assets into groups
    x_offset = 0.83;
    y_start = 0.55;
    y_step = 0.05;

    for i = 1:length(group_names)
        % Add a rectangle with the group color
        annotation('rectangle', 'Position', [x_offset, y_start - (i - 1) * y_step, 0.02, 0.02], ...
                'FaceColor', group_colors(i, :), 'EdgeColor', 'none');
        % Add the group name
        annotation('textbox', 'Position', [x_offset + 0.03, y_start - (i - 1) * y_step, 0.1, 0.02], ...
                'String', group_names{i}, 'LineStyle', 'none', 'FontSize', 10);
    end

    % Add axis labels and title
    xticks(1:size(ws, 2));
    xticklabels(names); 
    xlabel('Portfolio Index');
    ylabel('Weight');
    title('Portfolio Weights Distribution');

end