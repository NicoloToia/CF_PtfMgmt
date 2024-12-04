function [] = disp_pie_weights(weightsTable)
% This function displays pie charts of the asset allocations for each portfolio.

% Get the number of portfolios
numPortfolios = numel(weightsTable.Properties.VariableNames);

% Set two different windows for the pie charts of the portfolios
halfPortfolios = ceil(numPortfolios / 2);

% get the names of the portfolios and assets
portfolioNames = weightsTable.Properties.VariableNames; 
assetNames = weightsTable.Properties.RowNames;          

%%
    % First window
    figure('Name', 'Portfolio Allocations');
    rows = ceil(sqrt(halfPortfolios));
    cols = ceil(halfPortfolios / rows);

    for i = 1:halfPortfolios
        % Extract the weights for the current portfolio
        portfolioWeights = weightsTable{:, portfolioNames{i}};
        
        % Create a subplot
        subplot(rows, cols, i);
        pie(portfolioWeights, string(round(portfolioWeights * 100, 1)) + "%");
        title(['Asset Allocation ', portfolioNames{i}]);
    end

    % Add a shared vertical legend for the first window
    % Positioned to the right of the window
    legend(assetNames, 'Location', 'eastoutside', 'Orientation', 'vertical');
    sgtitle('Portfolio Allocations'); % General title of the window

%%
    % Second window
    figure('Name', 'Portfolio Allocations');
    rows = ceil(sqrt(numPortfolios - halfPortfolios));
    cols = ceil((numPortfolios - halfPortfolios) / rows);

    for i = halfPortfolios+1:numPortfolios
        % Extract the weights for the current portfolio
        portfolioWeights = weightsTable{:, portfolioNames{i}};
        
        % Create a subplot
        subplot(rows, cols, i - halfPortfolios);
        pie(portfolioWeights, string(round(portfolioWeights * 100, 1)) + "%");
        title(['Asset Allocation  ', portfolioNames{i}]);
    end

    % Add a shared vertical legend for the second window
    % Positioned to the right of the window
    legend(assetNames, 'Location', 'eastoutside', 'Orientation', 'vertical');
    sgtitle('Portfolio Allocations'); % General title of the window

end