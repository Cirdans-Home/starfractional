function [Q,H] = arnoldi(A,q1,m)
%ARNOLDI    Arnoldi iteration
%   [Q,H] = ARNOLDI(A,q1,M) carries out M iterations of the
%   Arnoldi iteration with N-by-N matrix A and starting vector q1
%   (which need not have unit 2-norm).  For M < N it produces
%   an N-by-(M+1) matrix Q with orthonormal columns and an
%   (M+1)-by-M upper Hessenberg matrix H such that
%   A*Q(:,1:M) = Q(:,1:M)*H(1:M,1:M) + H(M+1,M)*Q(:,M+1)*E_M',
%   where E_M is the M'th column of the M-by-M identity matrix.

%  Modified by SP

n = length(A);
if nargin < 3, m = n; end
q1 = q1/norm(q1);
Q = zeros(n,m); Q(:,1) = q1;
H = zeros(min(m+1,m),m);

for k=1:m
    z = A*Q(:,k);
    for i=1:k
        H(i,k) = Q(:,i)'*z;
        z = z - H(i,k)*Q(:,i);
    end
    if k < n
       H(k+1,k) = norm(z);
       if abs(H(k+1,k)) < 1e-12,
	  H = H(1:k,1:k);
	  Q = Q(:,1:k);
	  return
       end
       Q(:,k+1) = z/H(k+1,k);
   end
end

H = H(1:m,1:m);
Q = Q(:,1:m);


