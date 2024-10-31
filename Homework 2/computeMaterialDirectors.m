function [m1,m2] = computeMaterialDirectors(a1, a2, theta)
% Inputs:
% a1 is a matrix of size ne x 3. This contains the first
% reference director on each edge.
% a2 is a matrix of size ne x 3.
% theta is a vector of size ne ( theta = q(4:4:end) )
%
% Outputs:
% m1 and m2 are matrices of size ne x 3. Each column of
% these matrices contain the material directors on each
% edge.

ne = length(theta);
m1 = zeros(ne,3);
m2 = zeros(ne,3);

for c=1:ne % loop over edges
    cs = cos(theta(c));
    ss = sin(theta(c));
    m1(c,:) = cs * a1(c,:) + ss * a2(c,:);
    m2(c,:) = - ss * a1(c,:) + cs * a2(c,:);
end

end
