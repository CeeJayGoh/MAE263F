function refTwist = computeRefTwist(a1, tangent, refTwist)
% a1 is a matrix of size ne x 3. Each column of a1
% contains the first time-parallel reference director
% on each edge.
%
% tangent is a matrix of size ne x 3. It contains all
% the tangent vectors on the edges.
%
% refTwist is a vector of size nv (one reference twist
% for each node). This (optionally) is provided as an
% input (guess solution).

[ne, ~] = size(a1); % number of edges
nv = ne + 1;

% refTwist = zeros(nv, 1);

for c=2:ne % all internal nodes
    u0 = a1(c-1,:); % a1 vector of previous edge
    u1 = a1(c, :); % a1 vector of current edge
    t0 = tangent(c-1,:); % tangent of previous edge
    t1 = tangent(c,:); % tangent of current edge
    ut = parallel_transport(u0, t0, t1);
    % ut and u1 are the same?
    % Method 1: Okay? But 2*pi issue?
    % refTwist(c) = signedAngle(ut, u1, t1);
    % Method 2
    ut = rotateAxisAngle( ut, t1, refTwist(c) );
    refTwist(c) = refTwist(c) + signedAngle(ut, u1, t1);
end

end
