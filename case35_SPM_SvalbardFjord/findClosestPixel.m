function [row, column, closest_pixel] = findClosestPixel(lon, lat, lonref, latref)

% Subtract refernce points from grids
minDist = abs(lon - lonref) + abs(lat - latref);

% Find the index of the closest point
closest_pixel = find(minDist == min(abs(minDist(:))));

% return the row and column
[row, column] = find(minDist == min(abs(minDist(:))));

% site = [x, y];
% 
% row = site(1);
% column = site(2);

end
