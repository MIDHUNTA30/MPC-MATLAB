% MATLAB program for Noninear MPC Set point tracking: Simple pendulum system 
function main;
clear all;
close all
global N Np n m Ts x0 xr ur
NT=50;N=6;n=2;m=1; Ts = 0.1;
Q=eye(n); QN=Q; R=eye(m);
xr=[0.5;0];M=1;c=9.8;l=1;ur=M*c*l*sin(xr(1))
Fx=[1 0;0 1;-1 0;0 -1];gx=[10;10;10;10]-Fx*xr;
Fu=[1;-1];gu=[5;0]-Fu*ur;
                       
x0=[2;1];
x=zeros(n,NT+1); x(:,1)=x0-xr;
Xk=zeros(n*(N+1),1); Xk(1:n,1)=x0;
u=zeros(m,NT);
Uk=zeros(m*N,1);
zk=[Xk;Uk];

% constructing QX,RU,FX,gX,FU,gU,H
QX=Q;RU=R;
FX=Fx;gX=gx;FU=Fu;gU=gu;
for i=1:N-1
  QX=blkdiag(QX,Q); RU=blkdiag(RU,R);
  FX=blkdiag(FX,Fx);gX=[gX;gx];
  FU=blkdiag(FU,Fu);gU=[gU;gu];
end
QX=blkdiag(QX,QN);
FX=blkdiag(FX,Fx);
gX=[gX;gx];
H=blkdiag(QX,RU);



% simulating system with MPC
for k=1:NT
   x0=x(:,k);  
   fun = @(z)z'*H*z;
   F=blkdiag(FX,FU);g=[gX;gU];Feq=[];geq=[];lb=[];ub=[];nonlcon=@nlcon;
   z=fmincon(fun,zk,F,g,Feq,geq,lb,ub,nonlcon);
   u(k)=z((N+1)*n+1,1);
   x(:,k+1)=f(x(:,k),u(k));
   zk=z;
end    


% plotting response
figure(1)
time = (0:NT);
subplot(2,1,1)
plot(time,x(1,:)+xr(1),'r.-','LineWidth',.7) 
hold on
plot(time,x(2,:)+xr(2),'k.-','LineWidth',.7) 
axis ([0 50 -4 4])
legend('$x_1$','$x_2$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');ylabel('$\textbf{x}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:5:50])
set(gca,'ytick',[-4:2:4])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,1,2)
stairs(time(1:end-1),u+ur,'r.-','LineWidth',.7)
axis ([0 50 4 5])
xlabel('$k$','Interpreter','latex');ylabel('${u}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:5:50])
set(gca,'ytick',[4:.25:5])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg nmpc2

end

function xnext=f(x,u)
global NT N n m Ts x0 xr ur    
c=9.8;M=1;B=3;l=1;dt=0;
xnext=x+Ts*[x(2);-(c/l)*sin(x(1))-(B/M*l^2)*(x(2))+(1/M*l^2)*u+dt];
end

function [c,ceq] =nlcon(z)
global NT N n m Ts x0 xr ur
c=[];
for i=1:N
    c1((i-1)*n+1:i*n,1)=z(i*n+1:(i+1)*n,1)-f(z((i-1)*n+1:i*n,1)+xr,z((N+1)*n+(i-1)*m+1:(N+1)*n+i*m)+ur)+xr;
end
ceq = [z(1:2)-x0+xr;c1];
end



