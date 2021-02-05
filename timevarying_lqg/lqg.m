clc
clear

N =1e3;
A = 1.1;
B = 1;
Q = 1;
R = 1;
F = 1;

S = zeros(1,N);
S(N) = F;
L = zeros(1,N);
for k=N-1:-1:1
    S(k) = A'*(S(k+1) - S(k+1)*B*(B'*S(k+1)*B+R)^-1*B'*S(k+1))*A + Q;
    L(k) = (B'*S(k+1)*B+R)^-1*B'*S(k+1)*A;
end
% L = fliplr(L);
%%
C = 1;

x = zeros(N,1);
y = zeros(N,1);
u = zeros(N,1);
gain = zeros(N,1);
x_hat = zeros(N,1);
Kal = zeros(N,1);
P=zeros(N,1);
P(1) = 1;
x(1) = randn();
y(1) = randn();
u(1) = 0;
sigma2w = 1;
sigma2v = 1;
for k=1:N 
   y(k) = C  * x(k) + sqrt(sigma2v)*randn();
   Kal(k) = P(k)*C'*(C*P(k)*C' + sigma2w)^-1;
   P(k+1)= A*(P(k) - P(k)*C'*(C*P(k)*C' + sigma2w)^-1*C*P(k))*A' + sigma2v;
   gain(k) = L(k);

   u(k) = -gain(k)*x_hat(k);
   x_hat(k+1) = A * x_hat(k) + B*u(k) +Kal(k)*(y(k)- C*x_hat(k));;
   x(k+1) = A * x(k) + B*u(k)+ sqrt(sigma2w)*randn();
end

plot(1:N+1,x,'b',1:N+1,x_hat,'r');