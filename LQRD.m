function K=LQRD(A,B,Q,R,N)
Pk(:,:,N)=Q;
for i= N:-1:2
    Pk(:,:,i-1)=A'*(Pk(:,:,i)-Pk(:,:,i)*B*((B'*Pk(:,:,i)*B+R)^-1)*B'*Pk(:,:,i))*A+Q;
    K(:,i-1)=((B'*Pk(:,:,i)*B+R)^-1*B'*Pk(:,:,i)*A);
end
