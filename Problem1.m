%% HW1 Assignment 1
clear; clc; close all; 

N = 6;
RodLength = 0.1; % l = 1 m
deltaL = RodLength / (N-1);

% Radii of spheres
R1 = 0.005; % meter
R2 = 0.025; % meter
R3 = 0.005; % meter

% Density
rho_metal = 7000; % kg/m^3
rho_f = 1000; % fluid
rho = rho_metal - rho_f;
r0 = 0.001; % meter - rod radius
Y = 1e9; % Young's modulus (Y instead of E for clarity)
g = 9.8; % m/s^2 - gravity
visc = 1000; % pa-s
totalTime = 10; % second - total simulation time

% Utility parameter
EI = Y * pi * r0^4 / 4; % Nm^2 - bending stiffness
EA = Y * pi * r0^2; % Newton

% Geometry - initial configuration
nodes = zeros(N,2);

for c=1:N/2 % Loop over all the nodes
    nodes(c,1) = (c-1) * deltaL; % x coordinates
    nodes(c,2) = 0;
end

% Mass, M
M1 = 4/3*pi*R1^3*rho_metal;
M2 = 4/3*pi*R2^3*rho_metal;
M3 = 4/3*pi*R3^3*rho_metal;
mvec = [M1 M1 M2 M2 M3 M3];
M = diag(mvec);

% Viscous damping matrix, C
C1 = 6 * pi * visc * R1;
C2 = 6 * pi * visc * R2;
C3 = 6 * pi * visc * R3;
cvec = [C1 C1 C2 C2 C3 C3];
C = diag(cvec);

% Weight vector, W
W = zeros(N, 1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;

% Initial DOF
q0 = zeros(N, 1);

for c=1:N/2 % loop over nodes
    q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
    q0( 2*c ) = nodes(c,2); % y1, y2, y3
end

u0 = zeros(N, 1); % old velocity (initial velocity)

%% Prompt for implicit or Explicit 
method = input("Choose between implicit (1) or explicit (2): ");

if method == 1
    %Implicit Approach
    dt = 1e-2; % second - Time step size
    
    % Time marching scheme
    Nsteps = round(totalTime/dt);
    
    for c = 2:Nsteps
        fprintf('Time = %f\n', (c-1) * dt);
    
        % Guess
        q = q0; % New DOFs are initialized to be equal to old DOFs
    
        % Newton Raphson
        tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected
        err = 10 * tol;
    
        while err > tol
            f = M / dt * ( (q-q0)/dt - u0 );
            J = M / dt^2;
            
            % Elastic forces
            
            % Linear spring between nodes 1 and 2
            xk = q(1);
            yk = q(2);
            xkp1 = q(3);
            ykp1 = q(4);
            l_k = deltaL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            f(1:4) = f(1:4) + dF;
            J(1:4,1:4) = J(1:4,1:4) + dJ;
    
            % Linear spring between nodes 2 and 3
            xk = q(3);
            yk = q(4);
            xkp1 = q(5);
            ykp1 = q(6);
            l_k = deltaL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            f(3:6) = f(3:6) + dF;
            J(3:6, 3:6) = J(3:6, 3:6) + dJ;
    
            % Bending spring at node 2
            xkm1 = q(1);
            ykm1 = q(2);
            xk = q(3);
            yk = q(4);
            xkp1 = q(5);
            ykp1 = q(6);
            curvature0 = 0;
            l_k = deltaL;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            f(1:6) = f(1:6) + dF;
            J(1:6, 1:6) = J(1:6, 1:6) + dJ;
    
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
        all_mid_v(c) = u(4); %y velocity of middle sphere
        all_mid_q(c,1) = q(3); %x position of middle sphere
        all_mid_q(c,2) = q(4); %y positions

        % Plot
        figure(1);
        plot( q(1:2:end), q(2:2:end), 'ro-');
        axis equal
        title('Position of Spheres')
        xlabel('x [meter]');
        ylabel('y [meter]');
        drawnow
        % Update (new becomes old)
        q0 = q;
        u0 = u;
    
        % Store shape of structure
        if (c-1)*dt == 0.01
            q1 = q;
        elseif (c-1)*dt == 0.05
            q2 = q;
        elseif (c-1)*dt == 0.1
            q3 = q;
        elseif (c-1)*dt == 1
            q4 = q;
        elseif (c-1)*dt == 9.99
            q5 = q;
        end

    end
    % Plot middle node downward velocity
    figure(2);
    timeArray = (1:Nsteps) * dt;
    plot(timeArray, all_mid_v, 'k-');
    title('Velocity vs Time')
    xlabel('Time, t [sec]');
    ylabel('Velocity (vertical) of middle node, v [m/s]');
    
    % Plot middle node position vs time
    figure(3)
    plot(timeArray, all_mid_q(:,2), 'k-');
    title('Position of middle node vs Time')
    xlabel('Time, t [sec]');
    ylabel('Position of middle node,  [m]');

else
    % Explicit Approach
    dt = 1e-5; % second - Time step size
    Nsteps = round(totalTime / dt);

    % Initial Condition
    q = q0; % New DOFs are initialized to be equal to old DOFs

    for c = 2:Nsteps
        fprintf('Time = %f\n', (c - 1) * dt);

        % Compute elastic forces
        f = zeros(N, 1); % Initialize force vector to zero

        % Linear spring between nodes 1 and 2
        xk = q0(1);  % x position of node 1
        yk = q0(2);  % y position of node 1
        xkp1 = q0(3);  % x position of node 2
        ykp1 = q0(4);  % y position of node 2
        l_k = deltaL;  % Natural length of the spring
        dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);  % Gradient of stretching energy
        f(1:4) = f(1:4) + dF;  % Add elastic forces

        % Linear spring between nodes 2 and 3
        xk = q0(3);  % x position of node 2
        yk = q0(4);  % y position of node 2
        xkp1 = q0(5);  % x position of node 3
        ykp1 = q0(6);  % y position of node 3
        dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);  % Gradient of stretching energy
        f(3:6) = f(3:6) + dF;  % Add elastic forces

        % Bending spring at node 2
        xkm1 = q0(1);  % x position of node 1
        ykm1 = q0(2);  % y position of node 1
        xk = q0(3);  % x position of node 2
        yk = q0(4);  % y position of node 2
        xkp1 = q0(5);  % x position of node 3
        ykp1 = q0(6);  % y position of node 3
        curvature0 = 0;  % Initial curvature (zero for a straight beam)
        dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);  % Gradient of bending energy
        f(1:6) = f(1:6) + dF;  % Add bending forces

        % Viscous force (damping)
        f = f - C * (q-q0) / dt;  % Damping force is proportional to the velocity

        % Gravitational force (downward weight)
        f = f + W;  % Subtract weight (gravitational force)

        % Explicit velocity update
        u = u0 + (f ./ diag(M)) * dt;  % Update velocity

        % Explicit position update
        q = q0 + u * dt;  % Update position based on velocity

        % Store y velocity of the middle node (node 2)
        all_mid_v(c) = u(4);  % y-velocity of the middle sphere (node 2)
        all_mid_q(c, 1) = q(3);  % x position of the middle sphere (node 2)
        all_mid_q(c, 2) = q(4);  % y position of the middle sphere (node 2)

        % Plot positions at each time step
        figure(1);
        plot(q(1:2:end), q(2:2:end), 'ro-');  % Plot x and y coordinates
        axis equal;
        title('Position of Spheres (Explicit)');
        xlabel('x [meter]');
        ylabel('y [meter]');
        drawnow;

        % Update positions and velocities for the next step
        q0 = q;  % Current positions become the new initial positions
        u0 = u;  % Current velocities become the new initial velocities

    end

    % Plot middle node downward velocity
    figure(2);
    timeArray = (1:Nsteps) * dt;
    plot(timeArray, all_mid_v, 'k-');
    title('Velocity vs Time (Explicit)');
    xlabel('Time, t [sec]');
    ylabel('Velocity (vertical) of middle node, v [m/s]');

    % Plot middle node position vs time
    figure(3);
    plot(timeArray, all_mid_q(:, 2), 'k-');
    title('Position of middle node vs Time (Explicit)');
    xlabel('Time, t [sec]');
    ylabel('Position of middle node, [m]');
end

