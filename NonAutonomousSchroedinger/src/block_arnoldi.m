function [V, H, j, flag] = block_arnoldi(Afun, B, m)
% BLOCK_ARNOLDI: Block Arnoldi with reorthogonalization (double Gram-Schmidt)
%
% Inputs:
%   Afun - function handle for A*X
%   B    - initial block of vectors (n x s)
%   m    - number of block Arnoldi steps
%
% Outputs:
%   V    - Orthonormal basis of size (n x (m+1)*s)
%   H    - Block Hessenberg matrix of size ((m+1)*s x m*s)
%   j    - number of performed iterations
%   flag - flag = 0: no breakdown, flag = 1: breakdown

[n, s] = size(B);
V = zeros(n, (m+1)*s);
H = zeros((m+1)*s, m*s);
flag = 0;

% Orthonormalize the initial block B
[V(:,1:s), ~] = qr(B, 0);

for j = 1:m
    idx_j     = (j-1)*s + (1:s);
    idx_j1    = j*s + (1:s);

    % Apply A to the current block
    W = Afun(V(:, idx_j));

    % First orthogonalization pass
    for i = 1:j
        idx_i = (i-1)*s + (1:s);
        H_block = V(:, idx_i)' * W;
        H(idx_i, idx_j) = H_block;
        W = W - V(:, idx_i) * H_block;
    end

    % Reorthogonalization pass
    for i = 1:j
        idx_i = (i-1)*s + (1:s);
        H_block = V(:, idx_i)' * W;
        H(idx_i, idx_j) = H(idx_i, idx_j) + H_block;
        W = W - V(:, idx_i) * H_block;
    end

    % QR factorization of the new block
    [Q, R] = qr(W, 0);
    V(:, idx_j1) = Q;
    H(idx_j1, idx_j) = R;

    % Check for early termination
    if norm(R, 'fro') < 1e-12 %1e-14
        % warning('Block Arnoldi terminated early at step %d (rank deficiency)', j);
        V = V(:, 1:j*s);
        H = H(1:j*s, 1:(j-1)*s);
        flag = 1;
        return;
    end
end
end
