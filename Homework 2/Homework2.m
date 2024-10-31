%% Chapter 7: Discrete Elastic Rods Assignment

clear; clc; close all;

%% Global variables
global Fg M dt
global kappaBar EI GJ voronoiLength
global EA refLen

%% Inputs
% (1) Define the number of DOF related quantities
nv = 50; % number of nodes or vertices
ne = nv - 1; % number of edges
ndof = 3 * nv + ne; % number of DOF

% (2) Define time step related quantities
dt = 0.01; % time step size in seconds
totalTime = 5; % simulation time in seconds
Nsteps = round(totalTime / dt);

% (3) Physical parameters
% (3a) Geometry
RodLength = 0.2; % total length of the rod in meters
natR = 0.02; % natural radius in meters (initial curvature)
r0 = 0.001; % cross-sectional radius in meters

% (3b) Material parameters
Y = 10e6; % Young's modulus in Pascals
G = Y / 3; % Shear modulus for incompressible material
rho = 1000; % density in kg/m^3

g = [0; 0; -9.81]; % gravitational acceleration in m/s^2

%% Stiffness variables
EI = Y * pi * r0^4 / 4; % bending stiffness
GJ = G * pi * r0^4 / 2; % twisting stiffness
EA = Y * pi * r0^2; % stretching stiffness

%% Tolerance
tol = EI / RodLength^2 * 1e-6;

%% Mass Matrix
totalM = pi * r0^2 * RodLength * rho; % total mass in kg
dm = totalM / ne; % mass per edge
massVector = zeros(ndof, 1);

for c = 1:nv
    ind = [4 * c - 3; 4 * c - 2; 4 * c - 1]; % c-th node
    if c == 1 || c == nv
        massVector(ind) = dm / 2; % for first and last nodes
    else
        massVector(ind) = dm; % internal nodes
    end
end

for c = 1:ne
    ind = 4 * c;
    massVector(ind) = 1 / 2 * dm * r0^2; % moment of inertia for twist
end

M = diag(massVector);

%% Initial DOF vector
nodes = zeros(nv, 3);
dTheta = (RodLength / natR) * (1 / ne);

for c = 1:nv
    nodes(c, 1) = natR * cos((c - 1) * dTheta);
    nodes(c, 2) = natR * sin((c - 1) * dTheta);
    nodes(c, 3) = 0;
end

q0 = zeros(ndof, 1); % Initial condition

for c = 1:nv
    ind = [4 * c - 3; 4 * c - 2; 4 * c - 1]; % c-th node
    q0(ind) = nodes(c, :);
end

u = zeros(ndof, 1); % Initial velocity

%% Reference length for each edge
refLen = zeros(ne, 1);
for c = 1:ne
    dx = nodes(c + 1, :) - nodes(c, :); % dx = x_{c+1} - x_c
    refLen(c) = norm(dx);
end

%% Voronoi length
voronoiLength = zeros(nv, 1);
for c = 1:nv
    if c == 1
        voronoiLength(c) = 1 / 2 * refLen(c);
    elseif c == nv
        voronoiLength(c) = 1 / 2 * refLen(c - 1);
    else
        voronoiLength(c) = 1 / 2 * refLen(c - 1) + 1 / 2 * refLen(c);
    end
end

%% Reference frame (At t=0, initialize using space parallel transport)
a1 = zeros(ne, 3); % First reference director for all the edges
a2 = zeros(ne, 3); % Second reference director for all the edges
tangent = computeTangent(q0); % Tangent for all the edges

% Compute a1 for the first edge
t0 = tangent(1, :);
t1 = [0; 0; -1]; % arbitrary
a1Tmp = cross(t0, t1);

if abs(a1Tmp) < 1e-6
    t1 = [0; 1; 0]; % arbitrary
    a1Tmp = cross(t0, t1);
end

a1(1, :) = a1Tmp / norm(a1Tmp);
a2(1, :) = cross(tangent(1, :), a1(1, :));

% Done with the first edge
for c = 2:ne
    t0 = tangent(c - 1, :);
    t1 = tangent(c, :);
    a1_0 = a1(c - 1, :);
    a1_1 = parallel_transport(a1_0, t0, t1);
    a1(c, :) = a1_1 / norm(a1_1);
    a2(c, :) = cross(t1, a1(c, :));
end

%% Material frame
theta = q0(4:4:end); % Vector of size ne
[m1, m2] = computeMaterialDirectors(a1, a2, theta); % Compute material frame

%% Reference twist
refTwist = zeros(nv, 1); % initialize zero for first node
refTwist = computeRefTwist(a1, tangent, refTwist);

%% Natural curvature
kappaBar = getkappa(q0, m1, m2); % Natural curvature at each node

%% Gravity
Fg = zeros(ndof, 1);
for c = 1:nv
    ind = [4 * c - 3; 4 * c - 2; 4 * c - 1]; % c-th node
    Fg(ind) = massVector(ind) .* g;
end

%% Fixed and free DOFs
fixedIndex = 1:7; % Clamped boundary conditions
freeIndex = 8:ndof;

%% Time stepping scheme
ctime = 0; % Current time
endZ = zeros(Nsteps, 1); % z-coordinate of the last node with time

for timeStep = 1:Nsteps
    fprintf('Current time = %f\n', ctime);
    [q, u, a1, a2] = objfun(q0, u, a1, a2, freeIndex, tol, refTwist);
    ctime = ctime + dt;

    % Update q0 (old position)
    q0 = q;
    % Store endZ
    endZ(timeStep) = q(end);

    % Plot
    if mod(timeStep, 100) == 0
        theta = q(4:4:end);
        [m1, m2] = computeMaterialDirectors(a1, a2, theta);
        plotrod(q, a1, a2, m1, m2, ctime);
    end
end

% Visualization
figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, endZ, 'ro-');
xlabel('Time, t [sec]');
ylabel('z-coordinate of last node, \delta_z [m]');
