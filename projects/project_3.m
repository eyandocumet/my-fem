% Eyan Documet
clear all;
close all;
clc;

x0 = 0;      % Left endpoint location
L = 1;       % right endpoint location
Efunc = 0.4; % Youngs modulus function
BC0 = 0;      % u(x0) = 0
BCL = -6;      % u(L) = -6
res = 0;     % Number of resample points
p = 1;       % Shape function order
Ne = 20;     % 20 elements to start
h = (L-x0)/Ne * ones(Ne,1); % Array of initial element sizes
TOL = 0.05; % Error tolerance in energy norm

% Array to see if BC is Dirichlet (1) or Neumann (0)
% First entry (BCType(1)) is for left boundary
% Second entry (BCType(2)) is for right boundary
BCType = [1 1]; % Both Endpoints are Dirichlet for HW3

force = @(x) (Efunc .* 6 .* pi) .* (6 .* pi .* (x.^4 + 1) .* cos(3 .* pi .* x) .* sin(pi .* x.^3) + ...
             (4 .* x) .* (3 .* pi .* x .* sin(3 .* pi .* x) - cos(3 .* pi .* x)) .* cos(pi .* x.^3) + ...
             6 .* sin(3 .* pi .* x) + 9 .* pi .* x .* cos(3 .* pi .* x)); % Done by hand

uTrue = @(x) (4 .* sin(pi .* x.^3) + 6 .* x) .* cos(3 .* pi .* x); % Given

duTrue = @(x) cos(3 .* pi .* x) .* (12 .* pi .* x.^2 .* cos(pi .* x.^3) + 6) - ...
              3 .* pi .* sin(3 .* pi .* x) .* (4 .* sin(pi .* x.^3) + 6 .* x); % Done by hand


%% Problem 4(a-c), adaptive regime
disp("Solving Homework 3...");
tic;
% Initial mesh with Ne = 20
[xglobe, Nn, conn] = Mesh1D(p, Ne, x0, h);

% Solve FEM with initial Ne and get intial error
[xh, uN, error, PE] = myFEM1D(p, Ne, Nn, conn, xglobe, force, Efunc, BC0, BCL, BCType, duTrue, res);

NeVec = Ne; % Bookkeeping
PEVec = PE; % Ibid.
while any(error >= TOL) % Iterate until all nodes are compliant
    % Remesh
    [xglobe, conn, Ne, Nn] = addNodes(xglobe, error, TOL, p);

    % Resolve
    [xh, uN, error, PE] = myFEM1D(p, Ne, Nn, conn, xglobe, force, Efunc, BC0, BCL, BCType, duTrue, res);
    NeVec = [NeVec, Ne]; % Bookkeeping
    PEVec = [PEVec, PE]; % Ibid.
end

% Plotting -- Adaptive Meshing
figure;
hold on;

xtrue = linspace(0,1,1000); 

scatter(xh,zeros(1,numel(xh)),64,'+k','LineWidth',1.5) % 1D Mesh

plot(xtrue, uTrue(xtrue), 'LineWidth', 2); % True Solution

plot(xh, uN,'o--','LineWidth',2); % Approx. Solution

xlim([0,1]);

legend(["1D Mesh, Nodes", "True Solution - u", "FEM Solution - u^N"],'Location','northwest');

title("FEM h-Refinement, Adaptive Mesh");
ylabel('Displacement');
xlabel('Initial Position');
grid on;
hold off;

set(gcf, 'Color', 'w');  % Set background to white
saveas(gcf, 'adaptive_mesh.png', 'png');  % Save as PNG
print(gcf, 'adaptive_mesh.png', '-dpng', '-r300');  % Save with 300 dpi
close;

%% Problem 4(a-c), uniform regime

errFlag = false;
NeUniform = 20; % coarse initial mesh
NeVecUniform = [];
PEVecUniform = [];
while ~errFlag
    % Mesh and solve FEM problem, like Project 1
    h = 1/NeUniform * ones(NeUniform,1);
    
    [xglobeUni, Nn, conn] = Mesh1D(p, NeUniform, x0, h);
    
    [xh_uni, uN_uni, error_uni, PE] = myFEM1D(p, NeUniform, Nn, conn, xglobeUni, force, Efunc, BC0, BCL, BCType, duTrue, res);
   
    if max(error_uni) <= 0.05 % If error is above threshold, increment Ne
        errFlag = true;
    else
        NeUniform = NeUniform + 1;
          
    end
    NeVecUniform = [NeVecUniform, NeUniform]; % append new Ne to NeVec
    PEVecUniform = [PEVecUniform, PE]; % append new PE to PEVec
end

% Plotting -- Uniform mesh refinement
figure;
hold on;

xtrue = linspace(0,1,1000); 

scatter(xh_uni,zeros(1,numel(xh_uni)),64,'+k','LineWidth',1.5) % 1D Mesh

plot(xtrue, uTrue(xtrue), 'LineWidth', 2); % True Solution

plot(xh_uni, uN_uni,'o--','LineWidth',2); % Approx. Solution

xlim([0,1]);

legend(["1D Mesh, Nodes", "True Solution - u", "FEM Solution - u^N"],'Location','northwest');

title("FEM h-Refinement, Uniform Mesh");
ylabel('Displacement');
xlabel('Initial Position');
grid on;

set(gcf, 'Color', 'w');  % Set background to white
saveas(gcf, 'uniform_mesh.png', 'png');  % Save as PNG
print(gcf, 'uniform_mesh.png', '-dpng', '-r300');  % Save with 300 dpi
hold off;

close;

% Plotting -- PE vs. Ne
figure;
hold on;

plot(NeVec, PEVec, '-o', 'LineWidth',2);


plot(NeVecUniform, PEVecUniform, '-o', 'LineWidth', 2);
title("Potential Energy of an FEM Approximation as a Function of Element Count")
subtitle("Comparing Adaptive and Uniform Meshing")
xlabel('Number of Elements, N_e');
ylabel('Potential Energy, PE'); 
grid on;

set(gcf, 'Color', 'w');  % Set background to white
saveas(gcf, 'pe_ne.png', 'png');  % Save as PNG
print(gcf, 'pe_ne.png', '-dpng', '-r300');  % Save with 300 dpi
hold off;

close;

%% 4d
% Adaptive regime
XI = (xglobe(1:end-1) + xglobe(2:end)) / 2;
AI = error;

% Uniform regime
XI_uni = (xglobeUni(1:end-1) + xglobeUni(2:end)) / 2;
AI_uni = error_uni;

% Plotting -- XI and AI
figure;
hold on;
plot(XI, AI, '-o', 'LineWidth', 2);
plot(XI_uni, AI_uni, '-x', 'LineWidth', 2);

xlabel('Middle Position X_I');
ylabel('Error A_I');
legend('Adaptive Mesh', 'Uniform Mesh', 'Location','northwest');
grid on;
title('Error, AI, vs. Middle Position, XI, for Adaptive and Uniform FEM Analysis');

set(gcf, 'Color', 'w');  % Set background to white
saveas(gcf, 'error_trend.png', 'png');  % Save as PNG
print(gcf, 'error_trend.png', '-dpng', '-r300');  % Save with 300 dpi
hold off;

close;

% Bookkeeping -- Minimum elements in both regimes.
disp("min(Ne,adaptive) =");
disp(max(NeVec));

disp("min(Ne,uniform) =");
disp(max(NeVecUniform));

disp("Done!");
toc;

%% Mesh1D - 
% Generates a 1D mesh along a domain.
function [xglobe, Nn, conn] = Mesh1D(p, Ne, x0, h)

    % Assigning nodal spacings
    h_node = repelem(h, p) / p; 

    % Mesh vector
    xglobe = x0 + [0, cumsum(h_node)']; 
    [~, Nn] = size(xglobe);
    
    % Connectivity Matrix (Each row represents which global nodes correspond to the local node)
    conn = (1:p:p*(Ne))' + (0:p); 
end

%% evalShape - 
% Determines the shape functions & their derivatives.
function [ShapeFunc, ShapeDer] = evalShape(p, pts)
    switch p
        case 1 % Linear shape functions
            ShapeFunc = [(1-pts)./2, (1+pts)./2]; % Eq 3.27
            ShapeDer = [-1/2, 1/2].*ones(size(pts));
        
        case 2 % Quadratic Shape functions
            ShapeFunc = [pts.*(pts-1)/2, ...
                1-pts.^2, ...
                pts.*(pts+1)/2]; % Sec. 3.10

            ShapeDer = [pts - 1/2, ...
                        -2 .* pts, pts + 1/2];
            
        case 3 % Cubic Shape Functions
            ShapeFunc = [-9/16*(pts+1/3).*(pts-1/3).*(pts-1),...
                        27/16*(pts+1).*(pts-1/3).*(pts-1),...
                        -27/16*(pts+1).*(pts+1/3).*(pts-1),...
                        9/16*(pts+1).*(pts+1/3).*(pts-1/3)];

           ShapeDer = [-(27*pts.^2 - 18*pts - 1)/16,...
                      9*(9*pts.^2 - 2*pts - 3)/16,...
                      -9*(9*pts.^2 + 2*pts - 3)/16,...
                      (27*pts.^2 + 18*pts - 1)/16]; % Given
    end
end

function [xh, uN, error, PE] = myFEM1D(p, Ne, Nn, conn, xglobe, force, Efunc, BC0, BCL, BCType, duTrue, res)
    % Defining weights and Gauss points
    [wts, pts] = myGauss(p);

    % Evaluating shape functions and their derivatives
    [ShapeFunc, ShapeDer] = evalShape(p,pts);

    % Initializing stiffness Matrix
    K = zeros(Nn,Nn);

    % Initializing FEM solution vector
    uN = zeros(Nn,1);

    % Initializing Forcing vector
    R = zeros(Nn,1);

    % Iterating stiffness matrix
    for e = 1:Ne
        % Extracting element id from conn matrix
        id = conn(e, :); % Whole row-elem
        for gauss_pt = 1:numel(pts)
            % Jacobian
            J = xglobe(id) * ShapeDer(gauss_pt,:)';

            % Evaluating elemental stiffness matrix
            Ke = wts(gauss_pt).*ShapeDer(gauss_pt,:)'.*(1./J).*Efunc*ShapeDer(gauss_pt,:);
            K(id, id) = K(id, id) + Ke;

            % Mapping local to global coordinates
            x_zeta = xglobe(id) * ShapeFunc(gauss_pt,:)'; % Eq. 3.26

            % Evaluating forcing function & updating elemental loading term.
            Re = wts(gauss_pt)*ShapeFunc(gauss_pt,:)*force(x_zeta)*J;
            R(id) = R(id) + Re'; 
        end
    end

    KPE = K;
    RPE = R;

    % Boundary conditions
    if BCType(1) % Left Dirichlet BC
    uN(1) = BC0; % Set Dirichlet condition at the left

    % Adjust loading terms 
    switch p
        case 1 % Linear
            R(2) = R(2) - K(2,1) * BC0; % Modify the 2nd entry of R
        case 2 % Quadratic
            R(2) = R(2) - K(2,1) * BC0; % Modify R(2) for the 1st Dirichlet condition
            R(3) = R(3) - K(3,1) * BC0; % Modify R(3) as well for quadratic shape
        case 3 % Cubic
            R(2) = R(2) - K(2,1) * BC0; % Modify R(2) for cubic shape
            R(3) = R(3) - K(3,1) * BC0; % Modify R(3)
            R(4) = R(4) - K(4,1) * BC0; % Modify R(4)
        end
    else % Left Neumann BC
        R(1) = R(1) + BC0; % Neumann boundary modifies the load vector
    end

    if BCType(2) % Right Dirichlet BC
        uN(Nn) = BCL; % Set Dirichlet condition at the right

        % Adjust second to last loading terms
        switch p
            case 1 % Linear
            R(Nn-1) = R(Nn-1) - K(Nn-1, Nn) * BCL; % Modify R(Nn-1)
            case 2 % Quadratic
            R(Nn-1) = R(Nn-1) - K(Nn-1, Nn) * BCL; % Modify R(Nn-1)
            R(Nn-2) = R(Nn-2) - K(Nn-2, Nn) * BCL; % Modify R(Nn-2)
            case 3 % Cubic
            R(Nn-1) = R(Nn-1) - K(Nn-1, Nn) * BCL; % Modify R(Nn-1)
            R(Nn-2) = R(Nn-2) - K(Nn-2, Nn) * BCL; % Modify R(Nn-2)
            R(Nn-3) = R(Nn-3) - K(Nn-3, Nn) * BCL; % Modify R(Nn-3)
            end
    else % Right Neumann BC
        R(Nn) = R(Nn) + BCL; % Neumann boundary modifies the load vector
    end
    
    % Solving for uN (Removing the boundary conditions)
    freeNodes = 1:Nn;
    if BCType(1)
        freeNodes(1) = [];
    end
    if BCType(2)
        freeNodes(end) = [];
    end
    
    uN(freeNodes) = K(freeNodes, freeNodes) \ R(freeNodes);
    
    % Calculating PE
    PE = (1/2) * (uN' * KPE * uN) - (uN' * RPE);
    
    % Evaluating Error
    % Initialize error numerator and denominator
    errNum = zeros(Ne,1);
    errDen = 0;

    for e = 1:Ne
        id = conn(e,:); % Extract element ID
        % Loop through Gauss points, eval error
        for q = 1:numel(pts)

            x_zeta = xglobe(id) * ShapeFunc(q, :)';  % Eq. 3.26

            J = xglobe(id) * ShapeDer(q,:)'; 

            duN = uN(id)'*ShapeDer(q,:)'*(1/J); 

            errNum(e) = errNum(e) + wts(q)*(duTrue(x_zeta)-duN)*Efunc*(duTrue(x_zeta)-duN)*J;

            errDen = errDen + wts(q)*(duTrue(x_zeta))*Efunc*(duTrue(x_zeta))*J;
        end
        errNum(e) = errNum(e)/diff(xglobe(id)); % new
    end

    errDen = errDen/(xglobe(end)-xglobe(1));

    % Final error vector
    error = sqrt(errNum./errDen);
    
    % Resampling
    if res
        xh = zeros(res*Ne, 1);      % Initializing Positions of resample points
        uNres = zeros(res*Ne, 1);   % Initializing uN placeholder
        S = linspace(-1,1,res)';    % Sample points in the reference element ([-1, 1])
        [ResShapeFunc, ~] = evalShape(p,S); % ShapeFunc evaluated at res points
        
        n = 1; % Counter for indexing xh and uNRes through all res*Ne sample points
        for e = 1:Ne
            id = conn(e, :);
            for i = 1:res                        
                xh(n) = xglobe(id) * ResShapeFunc(i,:)'; % x(zeta)
                uNres(n) = uN(id)' * ResShapeFunc(i,:)'; % evaluating uh at sample point
                n = n + 1;
            end
        end
        uN = uNres;
    else
        xh = xglobe;
    end

end

function [xglobeNew, connNew, NeNew, NnNew] = addNodes(xglobe,error,TOL,p)

    naughty = find(error >= TOL); % Finding bad elements indices
    newNodes = zeros(size(naughty'));

    for i = 1:length(naughty)
        e = naughty(i);
        % Getting the midpoint of the "naughty" element, adding to "new
        % nodes" array
        bad_L = xglobe(e);
        bad_R = xglobe(e + 1);
        newNode = (bad_L + bad_R) / 2;
        newNodes(i) = newNode;
    end
    xglobeNew = sort([xglobe, newNodes]); % Create new xglobe and sort
    NeNew = length(xglobeNew) - 1; % See Mesh1D
    NnNew = length(xglobeNew); % See Mesh1D
    connNew = (1:p:p*(NeNew))' + (0:p); % See Mesh1D
end

%% DO NOT MODIFY BELOW -- use myGauss for myFEM1D
function [wts,pts] = myGauss(p)
    
    ptsNeed = ceil((p+1)/2);
    
    switch ptsNeed + 2
        case 1
            wts = 2;
            pts = 0;
        case 2
            wts = [1; 1];
            pts = [-0.5773502691896257; 0.5773502691896257];
        case 3
            wts = [0.8888888888888888; 0.5555555555555556; ...
                0.5555555555555556];
            pts = [0; -0.7745966692414834; 0.7745966692414834];
        case 4
            wts = [0.6521451548625461; 0.6521451548625461; ...
                0.3478548451374538; 0.3478548451374538];
            pts = [-0.3399810435848563; 0.3399810435848563; ...
                -0.8611363115940526; 0.8611363115940526];
        case 5
            wts = [0.5688888888888889; 0.4786286704993665;...
                0.4786286704993665; 0.2369268850561891; ...
                0.2369268850561891];
            pts = [0; -0.5384693101056831;  0.5384693101056831;...
                -0.9061798459386640; 0.9061798459386640];
    end
end