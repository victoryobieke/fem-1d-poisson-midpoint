% Finite Element Method for -u'' = f on (0,1), u(0)=u(1)=0
% Uses piecewise linear elements and midpoint quadrature
% Computes L2, H1 seminorm, and max nodal errors over h = 2^{-j}

clear; clc;

% Define exact solution, its derivative, and source function f
Exact       = @(x) sin(pi*x);
ExactDeriv  = @(x) pi*cos(pi*x);
Loadf       = @(x) pi^2*sin(pi*x);

% Mesh sizes
hvals = 2.^(-(1:10));
L2_err = zeros(size(hvals));
H1_err = zeros(size(hvals));
Max_err = zeros(size(hvals));

for jj = 1:length(hvals)
    h = hvals(jj);
 x = (0:h:1)';
    N = length(x)-1;

    % === Assemble stiffness matrix A ===
    A = zeros(N+1, N+1);
    for i = 1:N
        h_elem = x(i+1) - x(i);
        idx = [i i+1];
        A(idx, idx) = A(idx, idx) + [1 -1; -1 1]/h_elem;
    end

    % === Assemble load vector F using midpoint rule ===
    F = zeros(N+1, 1);
    for i = 1:N
        h_elem = x(i+1) - x(i);
        mid = (x(i) + x(i+1))/2;
        fmid = Loadf(mid);
        F([i i+1]) = F([i i+1]) + fmid * h_elem * [0.5; 0.5];
    end

    % === Apply Dirichlet boundary conditions and solve ===
    uh = zeros(N+1, 1);          % FEM solution
    uh(2:N) = A(2:N,2:N) \ F(2:N);

    % === Compute exact solution ===
    u_exact = Exact(x);

    % === Error: L2 norm (trapezoidal integration) ===
    L2_err(jj) = sqrt(trapz(x, (uh - u_exact).^2));

   % === Error: Energy norm (H1 semi-norm using exact form) ===
H1sq = 0;
for i = 1:N
    % h_elem = x(i+1) - x(i);
    xa = x(i);
    xb = x(i+1);
    hi = xb-xa;
    du = (uh(i+1) - uh(i)) / hi;
    der_enda = pi*cos(pi*xa);
    der_endb = pi*cos(pi*xb);
    % FEM derivative (piecewise constant)
    % mid = (x(i+1) + x(i)) / 2;
    % dtrue = ExactDeriv(mid);                   % exact derivative at midpoint
    % H1sq = H1sq + h_elem * (du - dtrue)^2;
    % Midpoint rule for energy norm
    H1sq = H1sq + hi/2*((der_enda - du)^2 + (der_endb-du)^2);
end
H1_err(jj) = sqrt(H1sq);


    % === Error: Maximum nodal error ===
    Max_err(jj) = max(abs(uh - u_exact));

end

% === Compute convergence orders ===
ordL2  = log2(L2_err(1:end-1) ./ L2_err(2:end));
ordH1  = log2(H1_err(1:end-1) ./ H1_err(2:end));
ordMax = log2(Max_err(1:end-1) ./ Max_err(2:end));
% ordEnergy = log2(E_err(1:end-1) ./ E_err(2:end));

% === Display error table ===
fprintf('   h\t\tL2 Error\tOrder\tH1 Error\tOrder\tMax Error\tOrder\n');
for j = 1:length(hvals)-1
    fprintf('%1.4f\t%.2e\t%.2f\t%.2e\t%.2f\t%.2e\t%.2f\n', ...
        hvals(j), L2_err(j), ordL2(j), H1_err(j), ordH1(j), Max_err(j), ordMax(j));
end
% Last row (no order)
j = length(hvals);
fprintf('%1.4f\t%.2e\t ----\t%.2e\t ----\t%.2e\t ----\n', ...
    hvals(j), L2_err(j), H1_err(j), Max_err(j));

% === Plot errors on log-log scale ===
figure;
loglog(hvals, L2_err, '-o', ...
       hvals, H1_err, '-s', ...
       hvals, Max_err, '-d', ...
       'LineWidth', 1.5); hold on;

% === Reference lines for first- and second-order ===
C1 = L2_err(1) / hvals(1)^2; % scaling for 2nd order
C2 = H1_err(1) / hvals(1)^1; % scaling for 1st order

loglog(hvals, C1*hvals.^2, '--k', 'LineWidth', 1.2); % 2nd order ref
loglog(hvals, C2*hvals, '--b', 'LineWidth', 1.2);    % 1st order ref

legend('L^2 Error', 'H^1 Error', 'Max Error', ...
       '2nd Order (h^2)', '1st Order (h)', ...
       'Location', 'southeast');
xlabel('Mesh size h'); ylabel('Error');
title('Error Plots');
grid off;
saveas(gcf, 'error_convergence.png');


% === Plot FEM and exact solution on finest mesh ===
figure;
plot(x, uh, '-', 'LineWidth', 1.5); hold on;
plot(x, u_exact, '--', 'LineWidth', 1.5);
legend('FEM Solution', 'Exact Solution', 'Location', 'best');
xlabel('x'); ylabel('u(x)');
ylim([0, 1.2])
title('FEM Approximation vs Exact Solution');
grid off;
saveas(gcf, 'fem_vs_exact.png');  % saves the figure to a PNG file
