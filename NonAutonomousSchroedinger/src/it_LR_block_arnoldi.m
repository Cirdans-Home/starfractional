function [u_leg_left, u_leg_right, svec, res, flag, left, right] = it_LR_block_arnoldi( ...
    B0, M0, K0, K1, u0vec, H1aphi, ...
    maxit, tol, nit, trunc, ...
    FF, Tal, Usch, Hal, t0, tf)

% Initialize
[d, m] = size(u0vec * H1aphi.');  % Estimate dimensions
left = zeros(d,1);
right = zeros(m,1);
leftold = u0vec;
rightold = H1aphi;
svec = [];
res = zeros(1, maxit);

for k = 1:maxit
    % Step 1: Construct Aleft and Aright
    Aleft = u0vec;
    Aright = H1aphi;
    Aleft = [Aleft, - B0*(M0\((K1-K0)*B0'*left))];
    Aright = [Aright, FF{2}*right];

    % Step 2: Truncated SVD
    if verLessThan('matlab','24.2')
        [Q, R] = qr(Aright, 0);
    else
        [Q, R] = qr(Aright, 'econ');
    end
    [UX, SX, VX] = svd(full(Aleft*R.'), 0);
    s_index = find(diag(SX) < trunc * SX(1,1), 1);
    if isempty(s_index)
        s_index = min(size(SX)) + 1;
    end
    s_index = s_index - 1;
    Aleft = UX(:, 1:s_index) * SX(1:s_index, 1:s_index);
    Aright = Q * conj(VX(:, 1:s_index));
    svec = [svec, s_index];

    % Step 3: Block Arnoldi + DLYAP
    A0fun = @(x) -B0*(M0\(K0*B0'*x));
    [VK0, J0, nita, flag(k)] = block_arnoldi(A0fun, Aleft, nit);
    nita = nita-1;
    J = J0(1:(nita)*s_index,1:(nita)*s_index);
    Y = dlyap(J, Tal.', VK0(:,1:(nita)*s_index)' * Aleft * (Aright.' * conj(Usch)));

    % Step 4: Update left, right
    left = VK0(:,1:(nita)*s_index);
    right = Usch * Y.';

    % Step 5: Residual
    res(k) = abs(u0vec.' * left * (abs(H1aphi).' * right).' - ...
                 u0vec.' * leftold * (abs(H1aphi).' * rightold).');

    if res(k) < tol
        break
    end

    % Prepare for next iteration
    leftold = left;
    rightold = right;
end

% Final solution
u_leg_left = Hal * right;
u_leg_right = left * 2 / (tf - t0);

end
