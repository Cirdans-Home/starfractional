
gam = -10;
alpha = 0.5;
M  = 40;
TM = basisfunction(M,alpha);
IM = eye(M);

R = (IM - gam*TM);
for i=1:M
  w(i,1)=(-1)^(i-1)*sqrt((2*(i-1)+1)/2);
end
x = R\w;
cml = TM*x;

%%
c = TM*w;
c = c(1:16);

c_cheb = leg2cheb(c,'norm');
fun = chebfun(c_cheb,'coeffs');

x = linspace(-1,1,1000);
figure(1)
plot(x,fun(x),'r-',...
    x,(x+1).^(alpha-1)/gamma(alpha),'b--','LineWidth',2);

%% reconstruct ml
figure(18)
semilogy(abs(cml),'o')

x = linspace(-1,1,1000);
c_cheb_ml = leg2cheb(cml(1:30),'norm');
fun = chebfun(c_cheb_ml,'coeffs');
figure(2)
plot(x,fun(x),'r-',...
    x,ml((x+1).^alpha*gam,alpha),'b--','LineWidth',2);
