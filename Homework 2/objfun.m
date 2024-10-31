function [q, u, a1, a2] = objfun(q0, u, a1, a2, freeIndex, tol, refTwist)
% Global
global Fg M dt

% Guess
q = q0;
iter = 1;
error = 10 * tol;

while error > tol
    % Compute reference frame
    [a1Iterate, a2Iterate] = computeTimeParallel(a1, q0, q);

    % Compute reference twist
    tangent = computeTangent(q);
    refTwist_iterate = computeRefTwist(a1Iterate, tangent, refTwist);

    % Compute material frame
    theta = q(4:4:end);
    [m1, m2] = computeMaterialDirectors( a1Iterate, a2Iterate, theta);

    % Compute force and Jacobian
    [Fb, Jb] = getFb(q, m1, m2); % Bending
    [Ft, Jt] = getFt(q, refTwist_iterate); % Twisting 
    [Fs, Js] = getFs(q); % Stretching

    % Equations of motion
    f = M/dt* ((q-q0)/dt - u) - Fb - Ft - Fs - Fg;
    J = M/dt^2 - Jb - Jt - Js;
    f_free = f(freeIndex);
    J_free = J(freeIndex, freeIndex);

    % Newton's update
    dq_free = J_free \ f_free ;
    q(freeIndex) = q(freeIndex) - dq_free;

    % Error
    error = sum( abs( f_free ) );
    fprintf('Iter = %d, error=%f\n', iter, error);
    iter = iter + 1;
end

u = (q - q0) / dt;
a1 = a1Iterate;
a2 = a2Iterate;
end