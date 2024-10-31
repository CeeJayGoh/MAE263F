function kappa = getkappa( q, m1, m2 )
% Inputs:
% q: DOF vector
% m1 is the first material director (size is ne x 3)
% m2 is the second material director (size is ne x 3)
% Output
% kappa is a matrix of size nv x 2. Each node has kappa_1 and kappa_2.
nv = ( length(q) + 1) / 4; % length(q) = 4*nv-1
ne = nv - 1;
kappa = zeros(nv, 2);

for c=2:ne % All internal nodes except first and last
    node0 = q(4*c-7:4*c-5);
    node1 = q(4*c-3:4*c-1);
    node2 = q(4*c+1:4*c+3);
    m1e = m1(c-1, :); % m1 vector of c-1-th edge
    m2e = m2(c-1, :); % m2 vector of c-1-th edge
    m1f = m1(c, :); % m1 vector of c-th edge
    m2f = m2(c, :); % m2 vector of c-th edge
    kappaLocal = computekappa(node0, node1, node2, m1e, m2e, m1f, m2f );
    kappa(c,1) = kappaLocal(1);
    kappa(c,2) = kappaLocal(2);
end

end
