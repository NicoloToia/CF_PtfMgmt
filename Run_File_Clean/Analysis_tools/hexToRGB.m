function rgb = hexToRGB(hexArray)
    n = length(hexArray);
    rgb = zeros(n, 3); % Preallocate for speed
    for i = 1:n
        hex = hexArray{i};
        if startsWith(hex, '#')
            hex = hex(2:end); % Remove the '#' character
        end
        rgb(i, :) = sscanf(hex, '%2x%2x%2x', [1 3]) / 255; % Convert to RGB
    end
end

