function [I,S] = generatePhaseDomain(maxS,maxI,numGrid,sumConstraint)
% Arguments are:
% - Maximum value of S
% - Maximum value of I
% - Number of points in the grid
% - The sum of S(0)+I(0) which remains constant

% Outputs I and S as an empty grid to be filled in

sValues = linspace(0,maxS,numGrid);
iValues = linspace(0,maxI,numGrid);
[S,I] = meshgrid(sValues,iValues);
maxConstraint = sumConstraint-(S+I);
out = maxConstraint<0;
I(out) = nan;
S(out) = nan;
end