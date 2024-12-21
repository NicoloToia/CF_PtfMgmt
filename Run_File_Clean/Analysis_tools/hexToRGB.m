function rgb = hexToRGB(hexArray)
% Convert hexadecimal color codes to RGB triplets
%
% INPUT:
% hexArray:     Cell array of hexadecimal colour codes
%
% OUTPUT:
% rgb:          Array of rgb colour codes

n = length(hexArray);   % Dimension of input array
rgb = zeros(n, 3);      % Preallocate for speed

% Convert to rgb each element
for i = 1:n
    hex = hexArray{i};
    if startsWith(hex, '#')
        hex = hex(2:end); % Remove the '#' character
    end
    rgb(i, :) = sscanf(hex, '%2x%2x%2x', [1 3]) / 255; % Convert to RGB
end

end

