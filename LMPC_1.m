% MATLAB program for Linear MPC: Single-input system
clear all;
close all
% System parameters and simulation parameters
A=[0.9 0.2;-0.4 0.8];
B=[0.1;0.01];
NT=50;N=5;n=2;m=1; 
Q=eye(n); QN=Q; R=eye(m);
xmin=-10;xmax=10;umin=-1;umax=1;

x0=[10;5]; 
x=zeros(n,NT+1); x(:,1)=x0;
Xk=zeros(n*(N+1),1); Xk(1:n,1)=x0;
u=zeros(m,NT);
Uk=zeros(m*N,1);
zk=[Xk;Uk];

% constructing AX,BU,QX,RU,H
for i=1:N+1
    AX((i-1)*n+1:i*n,:)=A^(i-1);
end
for i=1:N+1
  for j=1:N
      if i>j
          BU((i-1)*n+1:i*n,(j-1)*m+1:j*m)=A^(i-j-1)*B;
      else
          BU((i-1)*n+1:i*n,(j-1)*m+1:j*m)=zeros(n,m);
      end    
  end
end
QX=Q;RU=R;
for i=1:N-1
  QX=blkdiag(QX,Q); RU=blkdiag(RU,R);
end
QX=blkdiag(QX,QN);
H=blkdiag(QX,RU);

% simulating system with MPC
for k=1:NT
   xk=x(:,k);  
   fun = @(z)z'*H*z;
   F=[];g=[];Feq=[eye((N+1)*n) -BU];geq=AX*xk;
   lb=[xmin*ones(1,(N+1)*n),umin*ones(1,N*m)];
   ub=[xmax*ones(1,(N+1)*n),umax*ones(1,N*m)];
   z=fmincon(fun,zk,F,g,Feq,geq,lb,ub);
   u(:,k)=z((N+1)*n+1:(N+1)*n+m,1);
   x(:,k+1)=A*x(:,k)+B*u(:,k);
   zk=z;
end    

% plotting response
figure(1)
time = (0:NT);
subplot(2,1,1)
plot(time,x(1,:),'r.-','LineWidth',.7) 
hold on
plot(time,x(2,:),'k.-','LineWidth',.7) 
legend('$x_1$','$x_2$','Interpreter','latex');
%axis ([0 50-10 10])
xlabel('$k$','Interpreter','latex');ylabel('$\textbf{x}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:5:50])
set(gca,'ytick',[-10:5:10])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,1,2)
stairs(time(1:end-1),u,'r.-','LineWidth',.7)
%axis([0 50 -10 0])
xlabel('$k$','Interpreter','latex');ylabel('${u}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:5:50])
set(gca,'ytick',[-1:.5:1])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg lmpc1





