function gen_schur_hside(m,string)
% Generate a file with the data necessary 
% for the stra_frac_arnoldi function

fV{1} = @(t) 1 + 0*t;

trunc = 0;
[Flegvec, ~] = genCoeffMatrix_SP_interval(fV,m,[-1,1],trunc);
Hleg = sparse(Flegvec{2});

phi = zeros(m,1);
for i=1:m
    phi(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2);
end

[U,T] = schur(full(Hleg),'complex');

schur_hside.Uphi = U'*phi;
schur_hside.T = T;
schur_hside.U = U;
schur_hside.Hleg = Hleg;

name = [string,'_',num2str(m),'.mat'];

save(name,'schur_hside')

end