clc;
%clear all;

global Hu Qybar Smono Mmono M W S N ;

Ts=1/(50e3);
Cmpc=C;%eye(size(A,1));
Dmpc=D;%zeros(size(A,1),1);
[Ad,Bd,Cd,Dd]=c2dm(A,B,Cmpc,Dmpc,Ts);
m=size(Bd,2);
n=size(Ad,1);
r=size(Cd,1);
Amono=[Ad Bd;zeros(m,n) eye(m)];
Bmono=[Bd; eye(m)];
Cmono=[Cd zeros(r,m)];

N=10;
R=Ri;
Qy=1;%[.01 0 0 0; 0 .01 0 0 ;0 0 .01 0 ;0 0 0 .01 ];%eye(n);
%Primer caso
Qybar=[];
Mmono=[];
Tbar=[];
Smono=[];
Rbar=[];
for i=0:N-1
    row=zeros(size(Cmono,1),size(Bmono,2)*N);
    for j=0:i;
        row(:,j+1)=[Cmono*Amono^(i-j)*Bmono];
    end 
    Mmono=[Mmono;row];
    Qybar=blkdiag(Qybar,Qy);
    Rbar=blkdiag(Rbar,R);
    Smono=[Smono; Cmono*Amono^(i+1)];
end 

Umax=1-alpha;
Umin=0.5-alpha;
Mu=[eye(size(Bd,2)*N);-eye(size(Bd,2)*N)];
Ulim=[Umax*ones(size(Bd,2)*N,1);-Umin*ones(size(Bd,2)*N,1)];
M=[Mu];
W=[Ulim];
S=[zeros(size(Mu,1),n+1)];
size(Mmono)
size(Qybar)
size(Rbar)
Hu=Mmono'*Qybar*Mmono+Rbar;

u=[0.1;0.2;16;20;15;0.5];
X0=[u(1)];
X0=[X0;u(2)];
X0=[X0;u(3)];
X0=[X0;u(4)];
r=u(5)*ones(1,10);
ukm1=u(6);
X0=[X0;ukm1];
size(r*Qybar*Mmono)
size(X0'*Smono'*Qybar*Mmono)
Ubar= quadprog(Hu,X0'*Smono'*Qybar*Mmono-r*Qybar*Mmono,M,W+S*X0);
y=Ubar(1);

% DeltaU=-inv(Hu)*Gu1'*X0+inv(Hu)*Gu2'*r*ones(N,1);
% size(inv(Hu)*Gu2')
% size(r*ones(N,1))
% y=ukm1+DeltaU(1)
