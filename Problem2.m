%% HW1 Assignment 2
clear; clc; close all; 

% N_vals = [11, 21, 61, 101, 201];  % Different numbers of nodes to test
% terminal_velocities = zeros(size(N_vals));  % Store terminal velocities for each N

dt_vals = [1, 1e-3, 1e-6];  % Different time steps to test
terminal_velocities_dt = zeros(size(dt_vals));  % Store terminal velocities for each dt

for idx = 1:length(dt_vals) % Loops different N values or dt values
    % N = N_vals(idx);
    N = 21;  % Update number of nodes 

    % Number of nodes
    ndof = N * 2; % number of degrees of freedom
    dt = 0.01; % second - Time step size
    RodLength = 0.1; %
    deltaL = RodLength / (N-1);
    
    % Radii of spheres
    R = zeros(N,1); % Vector of size N - Radius of N nodes
    R(:) = deltaL/10;
    midNode = (N+1)/2;
    R(midNode) = 0.025;
    
    % Density
    rho_metal = 7000; % kg/m^3
    rho_f = 1000; % fluid
    rho = rho_metal - rho_f;
    r0 = 0.001; % meter - rod radius
    Y = 1e9; % Young's modulus (Y instead of E for clarity)
    g = 9.8; % m/s^2 - gravity
    visc = 1000; % pa-s
    totalTime = 50; % second - total simulation time
    
    % Utility parameter
    ne = N - 1; % number of edges
    EI = Y * pi * r0^4 / 4; % Nm^2 - bending stiffness
    EA = Y * pi * r0^2; % Newton
    
    % Geometry - initial configuration
    nodes = zeros(N,2);
    
    for c=1:N % Loop over all the nodes
        nodes(c,1) = (c-1) * deltaL; % x coordinates
        nodes(c,2) = 0;
    end
    
    % Mass, M
    M = zeros(ndof, ndof);
    
    for k=1:N
        M(2*k-1, 2*k-1) = 4/3*pi*R(k)^3*rho_metal; % Mass for x_k
        M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
    end
    
    % Viscous damping matrix, C
    C = zeros(ndof,ndof);
    
    for k=1:N
        C(2*k-1, 2*k-1) = 6 * pi * visc * R(k);
        C(2*k, 2*k) = C(2*k-1, 2*k-1);
    end
    
    % Weight vector, W
    W = zeros(ndof, 1);
    
    for k=1:N
        W(2*k-1) = 0; % weight along x is zero
        W(2*k) = -4/3*pi*R(k)^3*rho*g;
    end
    
    % Initial DOF
    q0 = zeros(ndof, 1);
    
    for c=1:N % loop over nodes
        q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
        q0( 2*c ) = nodes(c,2); % y1, y2, y3
    end
    
    u0 = zeros(ndof, 1); % old velocity (initial velocity)
    
    % tolerance
    tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected
    
    % Time marching scheme
    Nsteps = round(totalTime/dt);
    
    % Storage for y-velocity of the middle node
    all_mid_v = zeros(Nsteps, 1);
    all_mid_q = zeros(Nsteps, 1);
    
    for c = 2:Nsteps
        fprintf('Time = %f\n', (c-1) * dt);
        
        % Guess
        q = q0; % New DOFs are initialized to be equal to old DOFs
        
        % Newton Raphson
        err = 10 * tol;
        
        while err > tol
            f = M / dt * ( (q-q0)/dt - u0 );
            J = M / dt^2;
            %
            % Elastic forces
            %
            % Linear spring
            for k=1:N-1
                xk = q(2*k-1);
                yk = q(2*k);
                xkp1 = q(2*k+1);
                ykp1 = q(2*k+2);
                l_k = deltaL;
                dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
                dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
                ind = [2*k-1, 2*k, 2*k+1, 2*k+2];
                f(ind) = f(ind) + dF;
                J(ind,ind) = J(ind,ind) + dJ;
            end
    
            % Bending spring
            for k=2:N-1
                xkm1 = q(2*k-3);
                ykm1 = q(2*k-2);
                xk = q(2*k-1);
                yk = q(2*k);
                xkp1 = q(2*k+1);
                ykp1 = q(2*k+2);
                curvature0 = 0;
                l_k = deltaL;
                dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
                dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
                ind = [2*k-3, 2*k-2, 2*k-1, 2*k, 2*k+1, 2*k+2];
                f(ind) = f(ind) + dF;
                J(ind, ind) = J(ind, ind) + dJ;
            end
    
        % Viscous force
        f = f + C * (q-q0) / dt;
        J = J + C / dt;
    
        % Weight
        f = f - W;
    
        % Update
        q = q - J \ f;
        err = sum ( abs(f) );
    
        end
    
        % New velocity
        u = (q - q0) / dt;
    
        % Store some information
        all_mid_v(c) = u(2*midNode);
        all_mid_q(c) = q(2*midNode);
    
        % % Plot
        % figure(1);
        % plot( q(1:2:end), q(2:2:end), 'ro-');
        % axis equal
        % xlabel('x [meter]');
        % ylabel('y [meter]');
        % drawnow

        % Update (new becomes old)
        q0 = q;
        u0 = u;
    end
    % 
    % % Plot middle node downward velocity
    % figure(2);
    % timeArray = (1:Nsteps) * dt;
    % plot(timeArray, all_mid_v, 'k-');
    % xlabel('Time, t [sec]');
    % ylabel('Velocity (vertical) of middle node, v [m/s]');
    % 
    % % Plot middle node position vs time
    % figure(3)
    % plot(timeArray, all_mid_q, 'k-');
    % title('Position of middle node vs Time')
    % xlabel('Time, t [sec]');
    % ylabel('Position of middle node,  [m]');
    % 
    % Capture the terminal velocity after simulation
    terminal_velocities(idx) = all_mid_v(end);  % Assuming last velocity is terminal

end

% % Plot terminal velocity vs. number of nodes
% figure;
% plot(N_vals, terminal_velocities, 'o-');
% xlabel('Number of Nodes (N)');
% ylabel('Terminal Velocity [m/s]');
% title('Terminal Velocity vs. Number of Nodes');

% Plot terminal velocity vs. time step size
figure;
plot(dt_vals, terminal_velocities_dt, 'o-');
xlabel('Time Step Size (\Delta t)');
ylabel('Terminal Velocity [m/s]');
title('Terminal Velocity vs. Time Step Size');