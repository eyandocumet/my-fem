clear
clc

% Vector of elements we want to evaluate at.
NeVec = [4,16,64,256,1024];


% Vector of k constants to evaluate at.
kVec = [1,2,4,16,32];

% Parameters --
x0 = 0;         % Left Boundary
L = 1;          % Domain Length
p = 1;          % Element polynomial order (1 = Linear)
Efunc = 0.2;    % Stiffness Function
BC0 = 0;        % u(0) = 0
BCL = 1;        % u(L) = 1

% Array to see if BC is Dirichlet (1) or Neumann (0)
BCType = [1 1];
% First entry (BCType(1)) is for left boundary
% Second entry (BCType(2)) is for right boundary

% Equations --
force = @(x, k) -(k^2 .* sin((2 * pi * k * x) / L) + 2 * x.^2); % Given

uTrue = @(x, k) (-L^2) / (4 * Efunc * pi^2) .* sin(2 * pi * k * x / L) + ...
    (x.^4) / (6 * Efunc) + ((1 / L) + (L * sin(2 * pi * k) / (4 * Efunc * pi^2)) - ...
    (L^3 / (6 * Efunc))) .* x; % Solved by hand (See #1)

duTrue = @(x, k) -(L * k .* cos((2 * pi * k * x) / L)) / (2 * pi * Efunc) + ...
    (2 * x.^3) / (3 * Efunc) + ...
    (L * sin(2 * pi * k)) / (4 * pi^2 * Efunc) - ...
    L^3 / (6 * Efunc) + ...
    1 / L; % Solved w/ Symbolic Toolbox

disp("Solving HW #1...")
%% P3a
disp("(a) What is the minimum number of elements for various values of k?")
tic;
minNe = zeros(numel(kVec),1);
for i = 1:numel(kVec)   % loop through stiffness values
    k_i = kVec(i);

    errFlag = false;
    Ne = 4; % coarse initial mesh
    while ~errFlag
        % mesh and solve FEM problem
        h = 1/Ne * ones(Ne,1);
        [xglobe, Nn, conn] = Mesh1D(p, Ne, x0, h);
        [uN, error] = myFEM1D(p, Ne, Nn, conn, xglobe, @(x) force(x, k_i), Efunc, BC0, BCL, BCType, @(x) duTrue(x, k_i) );

        if error <= 0.05    % if error is above threshold, increment element count and try again
            minNe(i) = Ne;
            errFlag = true;
        else
            Ne = Ne + 1;
        end  
    end
end
toc;

min_elements_table = table(kVec', minNe, 'VariableNames', {'k', 'Min_Ne'})

%% P3b
disp("(b) Plot the numerical solutions for N = 4, 16, 64, 256, 1024, for each k, along with the true solution (make a plot u vs x for each k, with the true solution and the numerical solutions visible and discuss).")
tic;
close all
markers = {'o','+','*','s','d','v','>','h'};
% List a bunch of colors; like the markers, they 
% will be selected circularly. 
colors = {'b','g','r','k','m','c'};
% Same with line styles
linestyle = {'-','--','-.',':'};
% this function will do the circular selection
% Example:  getprop(colors, 7) = 'b'
getFirst = @(v)v{1}; 
getprop = @(options, idx)getFirst(circshift(options,-idx+1));


for j = 1:numel(kVec)
figure
hold on;

% Iterating k
k_j = kVec(j);

    for i = 1:numel(NeVec)
        Ne = NeVec(i);
        h = 1/Ne * ones(Ne,1);

        [xglobe, Nn, conn] = Mesh1D(p, Ne, x0, h);

        [uN, error] = myFEM1D(p, Ne, Nn, conn, xglobe, @(x) force(x, k_j), Efunc, BC0, BCL, BCType, @(x) duTrue(x, k_j) );

        plot(xglobe, uN,...
            'color',getprop(colors,i),...
            'Marker',getprop(markers,i),...
            'MarkerSize', 2,...
            'linestyle',getprop(linestyle,i),...
            'LineWidth',1,...
            'DisplayName', ['Ne = ', num2str(Ne)]);
    end

    xTrue = 0:0.001:1;
    plot(xTrue, uTrue(xTrue, k_j), 'DisplayName', 'True Sol','LineWidth',1);

    title(['k = ' num2str(k_j)], FontSize=24)
    xlabel('x');
    ylabel('uN');
    legend('Location','best')
    grid on;

    set(gcf, 'Color', 'w');  % Set background to white
    saveas(gcf, ['plot_k_' num2str(k_j) '.png'], 'png');  % Save as PNG
    print(gcf, ['plot_k_' num2str(k_j) '.png'], '-dpng', '-r300');  % Save with 300 dpi

end
toc;
disp('Saved to home directory.');

%% P3c - Not computationally efficient to run this again, but makes code more legible. In future projects, will probably refactor.
disp("Plot e^N (error) vs 1/N for each k (log-log)")
tic;
close all

figure

for i = 1:numel(kVec)
    k_i = kVec(i);  % Corrected this line

    error_k_i = zeros(1, numel(NeVec));

    for j = 1:numel(NeVec)
        Ne = NeVec(j);
        h = 1/Ne * ones(Ne, 1);

        [xglobe, Nn, conn] = Mesh1D(p, Ne, x0, h);

        [~, error] = myFEM1D(p, Ne, Nn, conn, xglobe, @(x) force(x, k_i), Efunc, BC0, BCL, BCType, @(x) duTrue(x, k_i));
        error_k_i(j) = error;  % Corrected this line
    end

    plot(1./NeVec, error_k_i, ...
        'color', getprop(colors, i), ...
        'Marker', getprop(markers, i), ...
        'MarkerSize', 5, ...
        'linestyle', getprop(linestyle, i), ...
        'LineWidth', 2, ...
        'DisplayName', ['k = ', num2str(k_i)]);

    yscale log;
    xscale log;


    hold on;
end

title('e^N (error) vs. 1/N for Different k Values')
xlabel('$\frac{1}{N}$', 'Interpreter', 'latex');
ylabel('$e^{N}$', 'Interpreter', 'latex');
legend('Location', 'best')
grid on;



set(gcf, 'Color', 'w');  % Set background to white



saveas(gcf, ['eN_vs_1_N.png'], 'png');  % Save as PNG
print(gcf, ['eN_vs_1_N.png'], '-dpng', '-r300');  % Save with 300 dpi

toc;
close all

disp('See chart.');

%% Mesh1D - 
% Generates a 1D mesh along a domain.
function [xglobe, Nn, conn] = Mesh1D(p, Ne, x0, h)
    Nn = p*Ne+1; % Number of nodes.
    Nne = p + 1; % Number of nodes per element

    % Initializing real domain positions
    xglobe = x0*ones(1,Nn);
    for node = 2:Nn
        xglobe(node) = xglobe(node-1) + h(node-1);
    end
    
    % Initializing connectivity matrix
    % <!> NOTE: Temporary for project 1. Will need to upgrade later. <!>
    conn = [(1:Nn-1)' (2:Nn)'];
end
%% myFEM1D - 
% Performs FEA on a 1D stick.
function [uN, error] = myFEM1D(p, Ne, Nn, conn, xglobe, force, Efunc, BC0, BCL, BCType, duTrue)
    % Defining weights and Gauss points
    [wts, pts] = myGauss(p);

    % Evaluating shape functions and their derivatives
    [ShapeFunc, ShapeDer] = evalShape(p,pts);

    % Initializing stiffness Matrix
    K = zeros(Nn, Nn);

    % Initializing FEM solution vector & forcing vector
    uN = zeros(Nn, 1);
    R = zeros(Nn, 1);

    % And now, the FEM
    % Loop through each element - 
    for elem = 1:Ne
        % Extract nodal indices (id) from conn matrix
        id = conn(elem, :); % Whole row-elem
        for q = 1:numel(pts) % Looping through each gauss point
            % Evaluating Jacobian
            J = xglobe(id) * ShapeDer(q,:)';

            % Evaluating & "stamping in" elemental stiffness matrix (2 x 2)
            Ke = wts(q) * ShapeDer(q,:)' * (1./J) * Efunc * ShapeDer(q, :);
            K(id, id) = K(id, id) + Ke;

            % Map gauss point from local to global coordinates
            x_zeta = xglobe(id) * ShapeFunc(q,:)'; 
            
            % Evaluating and assembling forcing function for elemental loading terms
            Re = wts(q)*ShapeFunc(q,:)*force(x_zeta)*J;
            R(id) = R(id) + Re';            
        end
    end
    
    % Boundary conditions
    if BCType(1) % Left Dirichlet BC
        uN(1) = BC0;
        % Adjust second loading term 
        R(2) = R(2) - K(2,1)*BC0;
    else % Nothing else for now
    end

    if BCType(2) % Right Dirichlet BC
        uN(Nn) = BCL;
        % Adjust second to last loading term 
        R(Nn-1) = R(Nn-1) - K(Nn-1, Nn)*BCL;
    else % Nothing else for now
    end
    
    % Calculating uN (with removed BC terms)
    uN(2:end-1) = K(2:end-1, 2:end-1) \ R(2:end-1);
    
    % Evaluating Error
    % Initialize error numerator and denominator
    errNum = 0;
    errDen = 0;

    for e = 1:Ne
        % Extract element ID
        id = conn(e, :);

        % Loop through Gauss points
        for q = 1:numel(pts)

            J = xglobe(id) * ShapeDer(q,:)'; % See above

            % Derivative of numerical solution
            duN = uN(id)' * ShapeDer(q,:)' * (1/J); % Given form DQ 4

            % Map gauss point from local to true spatial coordinates
            x_zeta = xglobe(id) * ShapeFunc(q,:)'; % See above

            % Error numerator and denominator
            errNum = errNum + wts(q) * (duTrue(x_zeta)-duN) * Efunc * (duTrue(x_zeta)-duN) * J;
            errDen = errDen + wts(q) * (duTrue(x_zeta)) * Efunc * (duTrue(x_zeta)) * J;
        end
    end
    
    % Final error
    error = sqrt(errNum/errDen);
end

%% evalShape - 
% Determines the shape functions & their derivatives.
function [ShapeFunc, ShapeDer] = evalShape(p,pts)
    switch p
        case 1
        % Linear Shape Functions (HW#1)
        ShapeFunc = [(1-pts)./2, (1+pts)./2]; % Eq 3.27
        ShapeDer = [-1/2, 1/2].*ones(size(pts));
    end
    % In future projects, other cases will be included.
end

%% myGauss - 
% Gaussian Quadrature Function; Quick Integration
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

%% GLOSSARY:
% p (int)                       : Polynomial Order of Shape Fnt.
% Ne (int)                      : Desired # Elements in Mesh or Number of Elements
% x0 (double)                   : Domain Start
% h (Ne x 1 double)             : Spacing of Elements        
% xglobe (1 x Nn double)        : Global Coordinates of Nodes 
% Nn (int)                      : Number of nodes in domain
% conn (Ne x Nn int)            : Connectivity Matrix
% force (function handle)       : Loading Function
% Efunc (int <!>)               : Function for Material Stiffness (E)
% BC0 (double)                  : Boundary Condition at x = 0 (Left End)
% BCL (double)                  : Boundary Condition at x = L (Right End)
% BCType (int)                  : Type of Boundary Condition (0 for Dirichlet, 1 for Neumann)
% duTrue (function handle)      : True Derivative of Displacement (for error computation)
% uN (Nn x 1 double)            : Approximated Solution at Nodes
% error (double)                : Error between True and Approximated Solution
% pts (1 x N double)            : Gauss points.
% ShapeFunc (N x (p+1) double)  : Shape function values at the points.
% ShapeDer  (N x (p+1) double)  : Derivatives of shape functions at the points.
% wts (1 x ptsNeed double)      : Gauss weights for integration.
