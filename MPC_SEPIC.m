function y=MPCsistema_online(u)
global Hu Qybar Smono Mmono M W S N;
X0=[u(1)];
X0=[X0;u(2)];
X0=[X0;u(3)];
X0=[X0;u(4)];
r=u(5)*ones(1,N);
ukm1=u(6);
X0=[X0;ukm1];
DeltaU= quadprog(2*Hu,2*(X0'*Smono'*Qybar*Mmono-r*Qybar*Mmono),M,W+S*X0);
y=ukm1+DeltaU(1);
