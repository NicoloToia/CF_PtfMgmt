function [] =  print_portfolio(weights,names,returns,risk,SR) 


disp('Asset Name                Weight')
disp('-------------------------------------------')
for i = 1:length(weights)
    fprintf('%-25s %.4f\n', names{i}, weights(i));
end
disp('-------------------------------------------')
fprintf('%-25s %.4f\n', 'Expected Return', returns);
fprintf('%-25s %.4f\n', 'Volatility', risk);
fprintf('%-25s %.4f\n', 'Sharpe Ratio', SR);
fprintf('%-25s %.4f\n', 'Sum of weights', sum(weights));
disp('  ')

end