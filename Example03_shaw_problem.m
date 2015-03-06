clear


load shawexamp.mat
A=G
b=dspiken

m=size(A,1)
n=size(A,2) % M should larger than N

%get the singular values
[u,s,v]=svd(A)

%set regularization parameters
lambda=zeros(1,1000)
for i=1:1000
    lambda(1,i)=10^(-10+10*i/1000)
end

ll=size(lambda,2);
x_lambda = zeros(n,ll);
rho = zeros(1,ll); 
eta = zeros(1,ll);
k=zeros(1,ll);

%calculate residual norm and solution norm
for i=1:ll
    etaderi=0;
    for j=1:n
        beta(j,1)=u(:,j)'*b;
        f(j)=s(j,j)^2/(s(j,j)^2+lambda(1,i)^2);
        xi(j,1)=beta(j,1)/s(j,j);
       xx(j,1)=beta(j,1)*f(j)/s(j,j);
    end
     etaderi=etaderi*(-4)/lambda(1,i);
     x_lambda(:,i)=v*xx;
    eta1(1,i)=norm(x_lambda(:,i));
    rho1(1,i)=norm(A*x_lambda(:,i)-b);
    eta(1,i)=eta1(1,i)^2;
    rho(1,i)=rho1(1,i)^2;
end

%calculate  curvatue and maximum absolute curcvature
k = lcfun(lambda,diag(s),beta,xi)
[kmax,knum]=max(abs(k))

%plot L-curve and curvature
figure(1)
subplot(1,2,1)
loglog(rho1,eta1);
title('L-curve')
xlabel('Residual norm ||Ax_\lambda-b||_2')
ylabel('Solution norm ||x_\lambda||')
grid on
subplot(1,2,2)
semilogx(lambda,k)
title('           Curvature of the L-curve')
xlabel('Regularization parameter \lambda')
ylabel('Curvature of the L-curve  \kappa')
%plot solution of the maximum abs(curvature)
figure(2)
stairs(x_lambda(:,knum));
xlabel('Element')
ylabel('Intensity')
