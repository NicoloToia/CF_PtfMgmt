function [] = Equally_weighted_pf(metricsTable)
    % Verifica che la tabella abbia almeno 3 righe
    if size(metricsTable, 1) < 3
        error('metricsTable deve contenere almeno 3 righe');
    end

    % Nomi delle etichette dei portafogli
    portfolioNames = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'M', 'N', 'P', 'Q', 'EW', 'CAPS'};
    
    % Verifica che il numero di colonne corrisponda al numero di etichette
    if size(metricsTable{1, :}, 2) ~= length(portfolioNames)
        error('Il numero di colonne di metricsTable deve essere %d per corrispondere ai nomi dei portafogli', length(portfolioNames));
    end

    % Genera una palette di colori per le barre
    numPortfolios = length(portfolioNames);
    colors = lines(numPortfolios); % 'lines' colormap for distinct colors
    
    % Crea una nuova figura
    figure;

    % Primo grafico: Expected Return
    subplot(2,2,1);
    b1 = bar(metricsTable{1,:}, 'FaceColor', 'flat'); 
    for i = 1:numPortfolios
        b1.CData(i,:) = colors(i,:); % Assegna un colore diverso a ciascuna barra
    end
    % Aggiungi la linea orizzontale corrispondente a EW
    hold on;
    plot([0:17],ones(18,1)*metricsTable{1,15} ,'r--', 'LineWidth', 0.5); % Linea rossa tratteggiata
    hold off;
    title('Expected Return');
    xlabel('Portfolios');
    ylabel('Value');
    xticks(1:numPortfolios);
    xticklabels(portfolioNames);
    xtickangle(45); 

    % Secondo grafico: Volatility
    subplot(2,2,2);
    b2 = bar(metricsTable{2,:}, 'FaceColor', 'flat'); 
    for i = 1:numPortfolios
        b2.CData(i,:) = colors(i,:); % Assegna un colore diverso a ciascuna barra
    end
    % Aggiungi la linea orizzontale corrispondente a EW
    hold on;
    plot([0:17],ones(18,1)*metricsTable{2,15}, 'r--', 'LineWidth', 0.5); % Linea rossa tratteggiata
    hold off;
    title('Volatility');
    xlabel('Portfolios');
    ylabel('Value');
    xticks(1:numPortfolios);
    xticklabels(portfolioNames);
    xtickangle(45); 

    % Terzo grafico: Sharpe Ratio
    subplot(2,2,3);
    b3 = bar(metricsTable{3,:}, 'FaceColor', 'flat'); 
    for i = 1:numPortfolios
        b3.CData(i,:) = colors(i,:); % Assegna un colore diverso a ciascuna barra
    end
    % Aggiungi la linea orizzontale corrispondente a EW
    hold on;
    plot([0:17],ones(18,1)*metricsTable{3,15}, 'r--', 'LineWidth', 0.5); % Linea rossa tratteggiata
    hold off;
    title('Sharpe Ratio');
    xlabel('Portfolios');
    ylabel('Value');
    xticks(1:numPortfolios);
    xticklabels(portfolioNames);
    xtickangle(45); 

    % Legenda dei portafogli (posizionata nel quarto quadrante)
    subplot(2,2,4);
    axis off; 
    hold on;
    for i = 1:numPortfolios
        plot(0.5, 1 - i * 0.05, 's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:)); 
        text(0.6, 1 - i * 0.05, portfolioNames{i}, 'FontSize', 10, 'HorizontalAlignment', 'left'); 
    end
    title('Portfolio Legend');
    hold off;
end
