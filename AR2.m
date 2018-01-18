
function [Q, B, V, E] = AR2( Y, L )

N=length(Y); y=Y(:,2); m=55;
b=zeros(3,1); R=eye(3); Q=0; Se=1;
B=zeros(N,3); V=B; E=B;

y0=y(1:m)-y(m)+y(1); %---------------------
for t=7:m
    z=[1 y0(t-1) y0(t-4)]';
    e=y0(t)-b'*z;
    R=L*R+z*z';
    b=b+inv(R)*z*e;
    u=y0(t)-b'*z;
    Se=L*Se+abs(1-L)*(e*u);
end
for t=7:m
    z=[1 y0(t-1) y0(t-4)]';
    e=y0(t)-b'*z;
    R=L*R+z*z';
    b=b+inv(R)*z*e;
    u=y0(t)-b'*z;
    Se=L*Se+abs(1-L)*(e*u);
end
for t=7:m
    z=[1 y0(t-1) y0(t-4)]';
    e=y0(t)-b'*z;
    R=L*R+z*z';
    b=b+inv(R)*z*e;
    u=y0(t)-b'*z;
    Se=L*Se+abs(1-L)*(e*u);
end

for t=7:N
    z=[1 y(t-1) y(t-4)]';
    e=y(t)-b'*z;
    Q = Q + e^2;
    R=L*R+z*z';
    b=b+inv(R)*z*e;
    B(t,:)= b';
    u=y(t)-b'*z;
    Se=L*Se+abs(1-L)*(e*u);
    E(t,:)=[Se, e, u];
    V(t,:)= diag(inv(R))'*Se/(1+L);
end 


