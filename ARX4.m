
function [Q, B, V, E, Vg] = ARX4( Y, L )

N=length(Y); y=Y(:,2); x=Y(:,3); m=55;
b=zeros(5,1); R=eye(5); Q=0; Se=1;
B=zeros(N,5); V=B; E=B(:,1:3);Vg=ones(N,1);

y0=y(1:m)-y(m)+y(1); %-----------------------
x0=x(1:m)-x(m)+x(1); % Inizializzazione ...
for t=7:m
    z=[1 y0(t-1) y0(t-4) x0(t-5) x0(t-6)]';
    e=y0(t)-b'*z;
    R=L*R+z*z';
    b=b+inv(R)*z*e;
    u=y0(t)-b'*z;
    Se=L*Se+abs(1-L)*(e*u);
end
for t=7:m
    z=[1 y0(t-1) y0(t-4) x0(t-5) x0(t-6)]';
    e=y0(t)-b'*z;
    R=L*R+z*z';
    b=b+inv(R)*z*e;
    u=y0(t)-b'*z;
    Se=L*Se+abs(1-L)*(e*u);
end
for t=7:m
    z=[1 y0(t-1) y0(t-4) x0(t-5) x0(t-6)]';
    e=y0(t)-b'*z;
    R=L*R+z*z';
    b=b+inv(R)*z*e;
    u=y0(t)-b'*z;
    Se=L*Se+abs(1-L)*(e*u);
end %----------------------------------------

for t=7:N
    z=[1 y(t-1) y(t-4) x(t-5) x(t-6)]';
    e=y(t)-b'*z;
    Q = Q + e^2;
    R=L*R+z*z';
    b=b+inv(R)*z*e;
    B(t,:)= b';
    u=y(t)-b'*z;
    Se=L*Se+abs(1-L)*(e*u);
    E(t,:)=[Se, e, u];
    P=inv(R)*Se/(1+L);
    V(t,:)= diag(P)';
    nu=sum(b(4:5));
    de=1-sum(b(2:3));
    gr = [nu/de nu/de 1 1]'/de;
    Vg(t) = gr'*(P(2:5,2:5))*gr;
end 


