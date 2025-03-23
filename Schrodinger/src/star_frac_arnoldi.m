function [u_leg,X] = star_frac_arnoldi(A,u0vec,alpha,int,max_it,m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

t0 = int(1); tf = int(2);

name = ['schur_hside_',num2str(m),'.mat'];
load(name)
Usch = schur_hside.U;
Tsch = schur_hside.T*(tf-t0)/2;
phisch = schur_hside.Uphi;
Hleg = schur_hside.Hleg*(tf-t0)/2;

Tal = Tsch^alpha;

[Q,Hes] = arnoldi(conj(A),u0vec,max_it);    % It works only for B = 0
[~, nit] = size(Hes);
e1n = speye(nit,1);
Y = dlyap(Tal,Hes(1:end-1,:)', phisch*norm(u0vec)*e1n.');
X = Usch*Y*Q(:,1:end-1)';
u_leg = Hleg*X*2/(tf-t0);

end