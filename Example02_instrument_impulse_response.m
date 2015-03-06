clear


noise=0.05;
N=211;
M=211;
t=linspace(-5,100,N);
%instrument response is a critically-damped pulse
sigi=10;
for i=1:N-1
  if (t(i)<0)
    g(i) = 0;
  else
    g(i) = t(i)*exp(-t(i)/sigi);
  end
end
gmax=max(g);
g=g/gmax;

for i=2:M
  for j=1:N-1
    tp=t(j)-t(i);
    if (tp > 0)
      G(i-1,j)=0;
    else
      G(i-1,j)=-tp*exp(tp/sigi);
    end
  end
end
deltat=t(2)-t(1);
G=G/gmax*deltat;
%true signal is two pulses of sig deviation
sig=2;
mtrue = (exp(-((t(1:N-1)-8).^2/(sig^2*2)))+0.5*exp(-((t(1:N-1)-25).^2/(sig^2*2))))';
mtrue=mtrue/max(mtrue);
d=G*mtrue;
nn=noise*randn(M-1,1);
dn=G*mtrue+nn;
A=G
b=dn

m=size(A,1)
n=size(A,2) % M should larger than N

%get the singular values
[u,s,v]=svd(A)

%set regularization parameters
lambda=zeros(1,1000)
for i=1:1000
    lambda(1,i)=10^(-5+7*i/1000)
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

%calculate  curvatue
k = lcfun(lambda,diag(s),beta,xi)

%plot true model
figure(1)
plot(t(1:N-1),mtrue)
xlabel('Time (s)')
ylabel('V')
%plot the predicted data set without noise
figure(2)
plot(t(1:N-1),d)
xlabel('Time (s)')
ylabel('V')
%plot the predicted data set with noise
figure(3)
plot(t(1:N-1),dn)
xlabel('Time (s)')
ylabel('V')
%plot L-curve and curvature
figure(4)
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
%plot the picking out instrument impulse solution
figure(5)
plot(t(1:N-1),x_lambda(:,620));
xlabel('Time (s)')
ylabel('Acceleratiion (m/s)')
