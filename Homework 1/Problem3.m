%% HW1 Assignment 3
clear; clc; close all; 

%% Inputs
% (1) Define the number of DOF related quantities

nv = 50; % Number of nodes or vertices
ne = nv - 1; % Number of edges
ndof = 2 * nv; % Number of DOF

% (2) Define time step related quantities
dt = 0.01; % Time step size in seconds
totalTime = 1; % Simulation time in seconds
Nsteps = round(totalTime / dt);

% (3) Physical parameters
% (3a) Geometry
L = 1; % Length of the beam in meters
outerRadius = 0.013; % Outer radius in meters
innerRadius = 0.011; % Inner radius in meters

% (3b) Material parameters
E = 70e9; % Young's modulus in Pa (GPa converted to Pa)
rho = 2700; % Density of aluminum in kg/m^3
P = -2000; % Point load applied at 0.75 m from the left-hand edge in N

% Moment of inertia of the cross-section
I = pi/4 * (outerRadius^4 - innerRadius^4); % Moment of inertia

%% Stiffness Variables
EA = E * pi * (outerRadius^2 - innerRadius^2); % Axial stiffness
EI = E * I; % Bending stiffness

%% Mass Matrix
% Compute mass per node based on density and geometry
m = pi * (outerRadius^2 - innerRadius^2) * L * rho / ne; % Mass per element
massVector = zeros(ndof, 1);

for c = 1:nv
    massVector(2 * c - 1:2 * c) = m / 2; % Divide mass evenly for each node
end

M = diag(massVector);

%% Initial DOF Vector
nodes = zeros(nv, 2);
deltaL = L / (nv - 1); % Length of each segment

for c = 1:nv
    nodes(c, 1) = (c - 1) * deltaL; % x-coordinates (initial positions along the beam)
    nodes(c, 2) = 0; % y-coordinates (initially all zero)
end

q0 = zeros(ndof, 1); % Initial condition

for c = 1:nv
    q0(2 * c - 1:2 * c) = nodes(c, :);
end

u0 = zeros(ndof, 1); % Initial velocity

%% Define Fixed and Free DOFs
% Constraints: First node fixed in both x and y, last node fixed in y
fixedIndex = [1, 2, ndof];
freeIndex = setdiff(1:ndof, fixedIndex);

%% External Force Vector
F_ext = zeros(ndof, 1);
% Apply point load at the node closest to 0.75 m from the left end
[~, loadNode] = min(abs(nodes(:, 1) - 0.75));
F_ext(2 * loadNode) = -P; % Apply load in the negative y direction

%% Tolerance for Newton-Raphson
tol = 1e-6;

%% Time Stepping Scheme
ctime = 0; % Current time
all_mid_q = zeros(Nsteps, 1); % To store vertical displacement of the middle node
ymax = zeros(Nsteps, 1); % To store maximum displacement

for timeStep = 1:Nsteps
    fprintf('Current time = %f\n', ctime);
    
    % Newton-Raphson Iteration
    err = 10 * tol;
    q = q0;
    while err > tol
        % Compute the residual force
        f = M / dt * ((q - q0) / dt - u0);
        J = M / dt^2;
        
        % Elastic Forces (Axial and Bending)
        for k = 1:ne
            % Axial spring between nodes k and k+1
            xk = q(2 * k - 1);
            yk = q(2 * k);
            xkp1 = q(2 * (k + 1) - 1);
            ykp1 = q(2 * (k + 1));
            l_k = deltaL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            ind = [2 * k - 1, 2 * k, 2 * (k + 1) - 1, 2 * (k + 1)];
            f(ind) = f(ind) + dF;
            J(ind, ind) = J(ind, ind) + dJ;
        end
        
        % Bending Spring Forces
        for k = 2:ne
            xkm1 = q(2 * (k - 1) - 1);
            ykm1 = q(2 * (k - 1));
            xk = q(2 * k - 1);
            yk = q(2 * k);
            xkp1 = q(2 * (k + 1) - 1);
            ykp1 = q(2 * (k + 1));
            curvature0 = 0;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            ind = [2 * (k - 1) - 1, 2 * (k - 1), 2 * k - 1, 2 * k, 2 * (k + 1) - 1, 2 * (k + 1)];
            f(ind) = f(ind) + dF;
            J(ind, ind) = J(ind, ind) + dJ;
        end
        
        % Apply external forces
        f = f + F_ext;
        
        % Update for fixed DOFs
        f(fixedIndex) = 0;
        J(fixedIndex, :) = 0;
        J(:, fixedIndex) = 0;
        J(fixedIndex, fixedIndex) = eye(length(fixedIndex));
        
        % Newton-Raphson update
        delta_q = -J \ f;
        q = q + delta_q;
        
        % Calculate error
        err = norm(f(freeIndex));

        % Plot beam shape
        figure(1);
        plot(q(1:2:end), q(2:2:end), 'b-o');
        axis equal;
        title('Beam Deformation vs Time');
        xlabel('x [m]');
        ylabel('y [m]');
    end
    
    % Update velocity
    u = (q - q0) / dt;
    
    % Store maximum displacement and middle node displacement
    ymax(timeStep) = max(abs(q(2:2:end)));
    all_mid_q(timeStep) = q(2 * round(nv / 2));
    
    % Update positions and velocities for the next step
    q0 = q;
    u0 = u;
    
    % Update current time
    ctime = ctime + dt;
end

%% Plot Maximum Vertical Displacement Over Time
figure;
timeArray = (1:Nsteps) * dt;
plot(timeArray, ymax, 'r-');
title('Maximum Vertical Displacement vs Time');
xlabel('Time [s]');
ylabel('Maximum Vertical Displacement [m]');

%% Plot Middle Node Vertical Displacement Over Time
figure;
plot(timeArray, all_mid_q, 'b-');
title('Middle Node Vertical Displacement vs Time');
xlabel('Time [s]');
ylabel('Vertical Displacement of Middle Node [m]');

%% Euler Beam Theory
d = 0.75; % Applied load position
c = L-d;
ymax_euler = (P*c*(L^2-c^2)^1.5)/(9*sqrt(3)*EI*L);

%% Plot P vs ymax for Different Loads
P_values = -2000:-1000:-20000; % Range of point loads
ymax_simulated = zeros(length(P_values), 1);
ymax_theory = zeros(length(P_values), 1);

for i = 1:length(P_values)
    P_current = P_values(i);
    F_ext(2 * loadNode) = -P_current;
    
    % Re-run simulation for the current load (using last time step result as initial condition)
    q = q0;
    for timeStep = 1:Nsteps
        err = 10 * tol;
        while err > tol
            f = M / dt * ((q - q0) / dt - u0);
            J = M / dt^2;
            
            for k = 1:ne
                xk = q(2 * k - 1);
                yk = q(2 * k);
                xkp1 = q(2 * (k + 1) - 1);
                ykp1 = q(2 * (k + 1));
                l_k = deltaL;
                dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
                dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
                ind = [2 * k - 1, 2 * k, 2 * (k + 1) - 1, 2 * (k + 1)];
                f(ind) = f(ind) + dF;
                J(ind, ind) = J(ind, ind) + dJ;
            end
            
            for k = 2:ne
                xkm1 = q(2 * (k - 1) - 1);
                ykm1 = q(2 * (k - 1));
                xk = q(2 * k - 1);
                yk = q(2 * k);
                xkp1 = q(2 * (k + 1) - 1);
                ykp1 = q(2 * (k + 1));
                curvature0 = 0;
                dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
                dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
                ind = [2 * (k - 1) - 1, 2 * (k - 1), 2 * k - 1, 2 * k, 2 * (k + 1) - 1, 2 * (k + 1)];
                f(ind) = f(ind) + dF;
                J(ind, ind) = J(ind, ind) + dJ;
            end
            
            f = f + F_ext;
            f(fixedIndex) = 0;
            J(fixedIndex, :) = 0;
            J(:, fixedIndex) = 0;
            J(fixedIndex, fixedIndex) = eye(length(fixedIndex));
            
            delta_q = -J \ f;
            q = q + delta_q;
            err = norm(f(freeIndex));
        end
    end
    ymax_simulated(i) = umax(abs(q(2:2:end)));
    ymax_theory(i) = (P_current * c * (L^2 - c^2)^1.5) / (9 * sqrt(3) * EI * L);
end

figure;
plot(P_values, ymax_simulated, 'b-o', P_values, ymax_theory, 'r--');
legend('Simulated', 'Euler Beam Theory');
title('P vs ymax Comparison');
xlabel('Point Load P [N]');
ylabel('Maximum Vertical Displacement ymax [m]');
