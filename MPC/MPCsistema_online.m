function y=MPCsistema_online(u)
% global Hu Gu1 Gu2 N;
% X0=[u(1)];
% X0=[X0;u(2)];
% X0=[X0;u(3)];
% X0=[X0;u(4)];
% r=u(5);
% ukm1=u(5);
% X0=[X0;ukm1];
% DeltaU=-inv(Hu)*Gu1'*X0+inv(Hu)*Gu2'*r*ones(N,1);
% y=ukm1+DeltaU(1);
    global H G M W S;
    X0=[u(1);u(2);u(3);u(4)];
    Ubar= quadprog(H,G*X0,M,W+S*X0);
    y=Ubar(1);