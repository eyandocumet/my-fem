clear
clc

% Vector of elements we want to evaluate at.
ne_vec = [4,16,64,256,1024];

% Vector of k constants to evaluate at.
k_vec = [1,2,4,16,32];

% Parameters --
x0 = 0;         % Left Boundary
L = 1;          % Domain Length
p = 1;          % Element polynomial order (1 = Linear)
efunc = 0.2;    % Stiffness Function
bc0 = 0;        % u(0) = 0
bcl = 1;        % u(L) = 1

% Array to see if BC is Dirichlet (1) or Neumann (0)
bc_type = [1 1];
% First entry (bc_type(1)) is for left boundary
% Second entry (bc_type(2)) is for right boundary

% Equations --
force = @(x, k) -(k^2 .* sin((2 * pi * k * x) / L) + 2 * x.^2); % Given

u_true = @(x, k) (-L^2) / (4 * efunc * pi^2) .* sin(2 * pi * k * x / L) + ...
    (x.^4) / (6 * efunc) + ((1 / L) + (L * sin(2 * pi * k) / (4 * efunc * pi^2)) - ...
    (L^3 / (6 * efunc))) .* x; % Solved by hand (See #1)

du_true = @(x, k) -(L * k .* cos((2 * pi * k * x) / L)) / (2 * pi * efunc) + ...
    (2 * x.^3) / (3 * efunc) + ...
    (L * sin(2 * pi * k)) / (4 * pi^2 * efunc) - ...
    L^3 / (6 * efunc) + ...
    1 / L; % Solved w/ Symbolic Toolbox

disp("Solving HW #1...")
%% P3a
disp("(a) What is the minimum number of elements for various values of k?")
tic;
min_ne = zeros(numel(k_vec),1);
for i = 1:numel(k_vec)   % loop through stiffness values
    k_i = k_vec(i);

    err_flag = false;
    ne = 4; % coarse initial mesh
    while ~err_flag
        % mesh and solve FEM problem
        h = 1/ne * ones(ne,1);
        [x_globe, nn, conn] = mesh1d(p, ne, x0, h);
        [u_n, error] = my_fem1d(p, ne, nn, conn, x_globe, @(x) force(x, k_i), efunc, bc0, bcl, bc_type, @(x) du_true(x, k_i) );

        if error <= 0.05    % if error is above threshold, increment element count and try again
            min_ne(i) = ne;
            err_flag = true;
        else
            ne = ne + 1;
        end  
    end
end
toc;

min_elements_table = table(k_vec', min_ne, 'VariableNames', {'k', 'Min_Ne'})

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
get_first = @(v)v{1}; 
getprop = @(options, idx)get_first(circshift(options,-idx+1));


for j = 1:numel(k_vec)
figure
hold on;

% Iterating k
k_j = k_vec(j);

    for i = 1:numel(ne_vec)
        ne = ne_vec(i);
        h = 1/ne * ones(ne,1);

        [x_globe, nn, conn] = mesh1d(p, ne, x0, h);

        [u_n, error] = my_fem1d(p, ne, nn, conn, x_globe, @(x) force(x, k_j), efunc, bc0, bcl, bc_type, @(x) du_true(x, k_j) );

        plot(x_globe, u_n,...
            'color',getprop(colors,i),...
            'Marker',getprop(markers,i),...
            'MarkerSize', 2,...
            'linestyle',getprop(linestyle,i),...
            'LineWidth',1,...
            'DisplayName', ['Ne = ', num2str(ne)]);
    end

    x_true = 0:0.001:1;
    plot(x_true, u_true(x_true, k_j), 'DisplayName', 'True Sol','LineWidth',1);

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
disp('See chart.');

%% P3c - Not computationally efficient to run this again, but makes code more legible. In future projects, will probably refactor.
disp("Plot e^N (error) vs 1/N for each k (log-log)")
tic;
close all

figure

for i = 1:numel(k_vec)
    k_i = k_vec(i);  % Corrected this line

    error_k_i = zeros(1, numel(ne_vec));

    for j = 1:numel(ne_vec)
        ne = ne_vec(j);
        h = 1/ne * ones(ne, 1);

        [x_globe, nn, conn] = mesh1d(p, ne, x0, h);

        [~, error] = my_fem1d(p, ne, nn, conn, x_globe, @(x) force(x, k_i), efunc, bc0, bcl, bc_type, @(x) du_true(x, k_i));
        error_k_i(j) = error;  % Corrected this line
    end

    plot(1./ne_vec, error_k_i, ...
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

%% mesh1d - 
% Generates a 1D mesh along a domain.
function [x_globe, nn, conn] = mesh1d(p, ne, x0, h)
    nn = p*ne+1; % Number of nodes.
    nne = p + 1; % Number of nodes per element

    % Initializing real domain positions
    x_globe = x0*ones(1,nn);
    for node = 2:nn
        x_globe(node) = x_globe(node-1) + h(node-1);
    end
    
    % Initializing connectivity matrix
    % <!> NOTE: Temporary for project 1. Will need to upgrade later. <!>
    conn = [(1:nn-1)' (2:nn)'];
end
%% my_fem1d - 
% Performs FEA on a 1D stick.
function [u_n, error] = my_fem1d(p, ne, nn, conn, x_globe, force, efunc, bc0, bcl, bc_type, du_true)
    % Defining weights and Gauss points
    [wts, pts] = my_gauss(p);

    % Evaluating shape functions and their derivatives
    [shape_func, shape_der] = eval_shape(p,pts);

    % Initializing stiffness Matrix
    k = zeros(nn, nn);

    % Initializing FEM solution vector & forcing vector
    u_n = zeros(nn, 1);
    r = zeros(nn, 1);

    % And now, the FEM
    % Loop through each element - 
    for elem = 1:ne
        % Extract nodal indices (id) from conn matrix
        id = conn(elem, :); % Whole row-elem
        for q = 1:numel(pts) % Looping through each gauss point
            % Evaluating Jacobian
            j = x_globe(id) * shape_der(q,:)';

            % Evaluating & "stamping in" elemental stiffness matrix (2 x 2)
            ke = wts(q) * shape_der(q,:)' * (1./j) * efunc * shape_der(q, :);
            k(id, id) = k(id, id) + ke;

            % Map gauss point from local to global coordinates
            x_zeta = x_globe(id) * shape_func(q,:)'; 
            
            % Evaluating and assembling forcing function for elemental loading terms
            re = wts(q)*shape_func(q,:)*force(x_zeta)*j;
            r(id) = r(id) + re';            
        end
    end
    
    % Boundary conditions
    if bc_type(1) % Left Dirichlet BC
        u_n(1) = bc0;
        % Adjust second loading term 
        r(2) = r(2) - k(2,1)*bc0;
    else % Nothing else for now
    end

    if bc_type(2) % Right Dirichlet BC
        u_n(nn) = bcl;
        % Adjust second to last loading term 
        r(nn-1) = r(nn-1) - k(nn-1, nn)*bcl;
    else % Nothing else for now
    end
    
    % Calculating u_n (with removed BC terms)
    u_n(2:end-1) = k(2:end-1, 2:end-1) \ r(2:end-1);
    
    % Evaluating Error
    % Initialize error numerator and denominator
    err_num = 0;
    err_den = 0;

    for e = 1:ne
        % Extract element ID
        id = conn(e, :);

        % Loop through Gauss points
        for q = 1:numel(pts)

            j = x_globe(id) * shape_der(q,:)'; % See above

            % Derivative of numerical solution
            du_n = u_n(id)' * shape_der(q,:)' * (1/j); % Given form DQ 4

            % Map gauss point from local to true spatial coordinates
            x_zeta = x_globe(id) * shape_func(q,:)'; % See above

            % Error numerator and denominator
            err_num = err_num + wts(q) * (du_true(x_zeta)-du_n) * efunc * (du_true(x_zeta)-du_n) * j;
            err_den = err_den + wts(q) * (du_true(x_zeta)) * efunc * (du_true(x_zeta)) * j;
        end
    end
    
    % Final error
    error = sqrt(err_num/err_den);
end

%% eval_shape - 
% Determines the shape functions & their derivatives.
function [shape_func, shape_der] = eval_shape(p,pts)
    switch p
        case 1
        % Linear Shape Functions (HW#1)
        shape_func = [(1-pts)./2, (1+pts)./2]; % Eq 3.27
        shape_der = [-1/2, 1/2].*ones(size(pts));
        
        % case 2
        % Quadratic Shape Functions (HW#2)

        % case 3
        % Cubic Shape Functions (HW#2)
    end
end

%% my_gauss - 
% Gaussian Quadrature Function; Quick Integration
function [wts, pts] = my_gauss(p)
    pts_need = ceil((p + 1) / 2);
    
    switch pts_need + 2
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
            wts = [0.5688888888888889; 0.4786286704993665; ...
                0.4786286704993665; 0.2369268850561891; ...
                0.2369268850561891];
            pts = [0; -0.5384693101056831;  0.5384693101056831; ...
                -0.9061798459386640; 0.9061798459386640];
    end
end


%% GLOSSARY:
% p            : Polynomial order of the shape function
% Ne           : Number of elements in the mesh
% x0           : Start of the domain
% h            : Element spacing (Ne x 1)
% x_globe      : Global node coordinates (1 x Nn)
% Nn           : Number of nodes in the domain
% conn         : Connectivity matrix (Ne x Nn)
% force        : Load function handle
% E_func       : Material stiffness function handle
% BC0          : Boundary condition at x = 0 (left boundary)
% BC_L         : Boundary condition at x = L (right boundary)
% BC_type      : Type of boundary condition (0 = Dirichlet, 1 = Neumann)
% du_true      : True displacement derivative (for error calculation)
% u_N          : Numerical solution at nodes (Nn x 1)
% error        : Error between true and numerical solution
% pts          : Gauss points (1 x N)
% shape_func   : Shape function values at Gauss points (N x (p+1))
% shape_der    : Derivatives of shape functions at Gauss points (N x (p+1))
% wts          : Gauss weights for integration (1 x pts_need)
