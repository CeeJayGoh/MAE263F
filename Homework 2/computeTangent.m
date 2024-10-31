function tangent = computeTangent(q)
% tangent is a matrix of size ne x 3
% q is a vector of size 4*nv-1
nv = (length(q)+1)/4;
ne = nv - 1;
tangent = zeros(ne, 3);
% Each column of "tangent" matrix is a 3-element
% vector representing the tangent on each edge

for c = 1:ne
    xc = q(4*c-3: 4*c-1);
    xcp1 = q(4*c+1: 4*c+3);
    edge = xcp1 - xc;
    tangent(c,:) = edge / norm(edge);
end

end
