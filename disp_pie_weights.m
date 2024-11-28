function [] = disp_pie_weights(weightsTable)

% Ottieni il numero di portafogli
numPortfolios = numel(weightsTable.Properties.VariableNames);

% Dividi il numero di portafogli per due finestre
halfPortfolios = ceil(numPortfolios / 2);

% Ottieni i nomi dei portafogli e degli asset
portfolioNames = weightsTable.Properties.VariableNames; % Nomi dei portafogli
assetNames = weightsTable.Properties.RowNames;          % Nomi degli asset

% Prima finestra
figure('Name', 'Portfolio Allocations');
rows = ceil(sqrt(halfPortfolios));
cols = ceil(halfPortfolios / rows);

for i = 1:halfPortfolios
    % Estrai i pesi per il portafoglio corrente
    portfolioWeights = weightsTable{:, portfolioNames{i}};
    
    % Crea un subplot
    subplot(rows, cols, i);
    pie(portfolioWeights, string(round(portfolioWeights * 100, 1)) + "%");
    title(['Asset Allocation ', portfolioNames{i}]);
end

% Aggiungi una legenda verticale condivisa per la prima finestra
% Posizionata a destra della finestra
legend(assetNames, 'Location', 'eastoutside', 'Orientation', 'vertical');
sgtitle('Portfolio Allocations'); % Titolo generale della finestra

% Seconda finestra
figure('Name', 'Portfolio Allocations');
rows = ceil(sqrt(numPortfolios - halfPortfolios));
cols = ceil((numPortfolios - halfPortfolios) / rows);

for i = halfPortfolios+1:numPortfolios
    % Estrai i pesi per il portafoglio corrente
    portfolioWeights = weightsTable{:, portfolioNames{i}};
    
    % Crea un subplot
    subplot(rows, cols, i - halfPortfolios);
    pie(portfolioWeights, string(round(portfolioWeights * 100, 1)) + "%");
    title(['Asset Allocation  ', portfolioNames{i}]);
end

% Aggiungi una legenda verticale condivisa per la seconda finestra
% Posizionata a destra della finestra
legend(assetNames, 'Location', 'eastoutside', 'Orientation', 'vertical');
sgtitle('Portfolio Allocations'); % Titolo generale della finestra