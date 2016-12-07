function [combinedData] = myColorRatio(red,green,blue)
% myColorRatio calculates the ratio of each color to the total for ternary
% diagrams.

total = red + green + blue;
R = red./total;
G = green./total;
B = blue./total;

combinedData = [R, G, B];
end

