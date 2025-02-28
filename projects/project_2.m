clear;
clc;

k = 6; % Frequency constant
x0 = 0.1;   % Left endpoint location
L = 1.2;    % right endpoint location
Efunc = 0.2;% Youngs modulus function
BC0 = 1;    % u(x = x_0 = 0.1) = 1
BCL = -0.7; % Efunc*du/dx = -0.7
res = 10;   % Number of resample points

% Array to see if BC is Dirichlet (1) or Neumann (0)
% First entry (BCType(1)) is for left boundary
% Second entry (BCType(2)) is for right boundary
BCType = [1 0]; % Left boundary Dirichlet, right boundary Neumann

force = @(x) -k^2*sin(pi*k*x/L) - k*cos(2*pi*k*x/L);

c1 = -0.7 + (k .* L ./ pi) .* cos(pi .* k) - (L ./ (2 .* pi)) .* sin(2 .* pi .* k);
c2 = 0.2 + (L.^2 ./ pi.^2) .* sin((pi .* k .* x0) ./ L) + (L.^2 ./ (4 .* pi.^2 .* k)) .* cos((2 .* pi .* k .* x0) ./ L) - c1 .* x0;

uTrue = @(x) ((-L.^2 ./ pi.^2) .* sin((pi .* k .* x) ./ L) - (L.^2 ./ (4 .* pi.^2 .* k)) .* cos((2 .* pi .* k .* x) ./ L) + c1 .* x + c2) ./ Efunc;

duTrue = @(x) ((-k.*L./pi.*cos(pi.*k.*x./L)+L./(2.*pi).*sin(2.*pi.*k.*x./L)) + c1 ) ./ Efunc;

% Stuff to make plots look nice
markers = {'o','+','*','s','d','v','>','h'};
% List a bunch of colors; like the markers, they will be selected circularly. 
colors = {'b','g','r','k','c','m'};
% Same with line styles
linestyle = {'-','--','-.',':'};
% this function will do the circular selection
% Example:  getprop(colors, 7) = 'b'
getFirst = @(v)v{1}; 
getprop = @(options, idx)getFirst(circshift(options,-idx+1));
% End Plotting Stuff

%% Problem 3a Filling in Chart
tic;
disp("Problem (a)")
pVec = [1, 2, 3]; % Orders of the polynomial
minNe = zeros(numel(pVec),1);
errorVec = zeros(numel(pVec),1);

for i = 1:numel(pVec)
    p_i = pVec(i);

    errFlag = false;
    Ne = 4; % coarse initial mesh
    while ~errFlag
        % mesh and solve FEM problem
        h = (L-x0)/Ne * ones(Ne,1);
        
        % Discretize domain
        [xglobe, conn] = Mesh1D(p_i, Ne, x0, h);

        [~, ~, error] = myFEM1D(p_i, Ne, conn, xglobe, force, Efunc, BC0, BCL, BCType, duTrue, res);

        if error <= 0.05    % if error is above threshold, increment element count and try again
            minNe(i) = Ne;
            errorVec(i) = error;
            errFlag = true;
        else
            Ne = Ne + 1;
        end  
    end
end

% Store results in a table and print it
errorTable = table(pVec(:), minNe(:), errorVec(:), 'VariableNames', {'p', 'Ne', 'Error'})

%% Problem 3b Plotting for each p
disp("Problem (b)")
NeVec = [6, 12, 64];
for p = pVec  % Looping through each element polynomial order (1 = Linear, 2 = Quadratic, 3 = Cubic)
    figure(p) % Creating a plot for each p
    for i = 1:numel(NeVec) % Loop through each value of elements
        Ne = NeVec(i);

        h = (L-x0)/Ne * ones(Ne,1); % Spacing for each element

        % Discretize domain
        [xglobe, conn] = Mesh1D(p, Ne, x0, h);

        % Solve FEM problem
        [xh, uN, error] = myFEM1D(p, Ne, conn, xglobe, force, Efunc, BC0, BCL, BCType, duTrue, res);

        % Plot Solution
        hold on;
        plot(xh, uN,...
            'Marker',getprop(markers,i),...
            'MarkerSize',4,...
            'color',getprop(colors,i),...
            'linestyle',getprop(linestyle,i),...
            'DisplayName', ['Ne = ', num2str(Ne)]);
    end
    %Plotting the true solution
    xTrue = linspace(x0,L,1000);
    plot(xTrue, uTrue(xTrue), 'DisplayName', 'True Sol','LineWidth',3);
    title(['p = ', num2str(p)], FontSize=24)
    xlabel('x', FontSize=16);
    ylabel('uN', FontSize=16);
    xlim([0,1.3]);
    ylim([0,11.4]);
    legend('Location','best')
    grid on;
    hold off;

    set(gcf, 'Color', 'w');  % Set background to white
    saveas(gcf, ['p_' num2str(p) '.png'], 'png');  % Save as PNG
    print(gcf, ['eN_vs_1_N.png'], '-dpng', '-r300');  % Save with 300 dpi
    close;
end

%% Problem 3c Getting error plots
disp("Problem (c)")
NeVec = [16 32 64 128 256 512];
errorVec = zeros(1,numel(NeVec));
figure();

for p = 1:3
    for i = 1:numel(NeVec)
        Ne = NeVec(i);

        h = (L-x0)/Ne * ones(Ne,1); % Spacing for each element

        % Discretize domain
        [xglobe, conn] = Mesh1D(p, Ne, x0, h);

        % Solve FEM problem and get error
        [~, ~, errorVec(i)] = myFEM1D(p, Ne, conn, xglobe, force, Efunc, BC0, BCL, BCType, duTrue, res);
    end

    plot(1./NeVec, errorVec,...
        'Marker',getprop(markers,p),...
        'color',getprop(colors,p),...
        'linestyle',getprop(linestyle,p),...
        'DisplayName', ['p = ', num2str(p)],'LineWidth',1.5);
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    hold on;
    title('Error vs 1/N_e', FontSize=24)
    xlabel('Inverse Number of Elements, N_e^{-1}', FontSize=16);
    ylabel('Error, Energy Norm', FontSize=16);
    legend('Location','best')
    grid on;
end
set(gcf, 'Color', 'w');  % Set background to white
saveas(gcf, ['eN_vs_1_N.png'], 'png');  % Save as PNG
print(gcf, ['eN_vs_1_N.png'], '-dpng', '-r300');  % Save with 300 dpi
close;
toc;
%% myFEM1D - 
% Performs FEA on a 1D stick.
function [xh, uN, error] = myFEM1D(p, Ne, conn, xglobe, force, Efunc, BC0, BCL, BCType, duTrue, res)
    % Number of nodes
    Nn = Ne*p + 1;

    % Defining weights and Gauss points
    [wts, pts] = myGauss(p); 

    % Evaluating shape functions and their derivatives
    [ShapeFunc, ShapeDer] = evalShape(p,pts);

    % Initializing matrices
    K = zeros(Nn, Nn); % Stiffness Matrix
    uN = zeros(Nn, 1); % Solution Vector
    R = zeros(Nn, 1); % Initializing Forcing vector

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

    % Evaluating Error --
    errNum = 0;
    errDen = 0;

    for e = 1:Ne
        % Extract element ID
        id = conn(e, :);

        % Loop through Gauss points
        for q = 1:numel(pts)
            % Jacobian (see above)
            J = xglobe(id) * ShapeDer(q,:)';

            % Derivative of numerical solution
            duN = uN(id)' * ShapeDer(q,:)' * (1/J); % Given form DQ4

            % Mapping gauss point from local to true spatial coordinates
            x_zeta = xglobe(id) * ShapeFunc(q,:)'; % See above

            % Error numerator and denominator
            errNum = errNum + wts(q) * (duTrue(x_zeta)-duN) * Efunc * (duTrue(x_zeta)-duN) * J;
            errDen = errDen + wts(q) * (duTrue(x_zeta)) * Efunc * (duTrue(x_zeta)) * J;
        end
    end
    % Final error - Energy norm
    error = sqrt(errNum/errDen);
    
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

%% Mesh1D - 
% Generates a 1D mesh along a domain.
function [xglobe, conn] = Mesh1D(p, Ne, x0, h)

    % Assigning nodal spacings
    h_node = repelem(h, p) / p; 

    % Mesh vector
    xglobe = x0 + [0, cumsum(h_node)']; 

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
