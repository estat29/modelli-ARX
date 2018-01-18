
X=load('Clima.txt');

ll=linspace(.9,1,20);
for i=1:20
Ql(i)=ARX4(X,ll(i));
end

figure 
plot(ll,Ql,'.-')
ii=find(Ql==min(Ql));
L=ll(ii)  %  L=.96;

N=length(X); tt=X(7:N,1);
[qq, bb, vv, ee] = AR2( X, L );
[Qq, Bb, Vv, Ee, Vg] = ARX4( X, L );

B1=Bb*0+2*sqrt(Vv);
B2=Bb*0-2*sqrt(Vv);
G=sum(Bb(:,4:5),2)./(ones(N,1)-sum(Bb(:,2:3),2));
F=( (ee(:,1)-Ee(:,1))/(2*(1-L)) )./( Ee(:,1)/((1+L)-5*(1-L)) ); % eq 17a G96

figure
subplot(221); hold on
plot(tt,Bb(7:N,2))
plot(tt,B1(7:N,2),':k')
plot(tt,B2(7:N,2),':k')
title('\phi_1')
subplot(222); hold on
plot(tt,Bb(7:N,3))
plot(tt,B1(7:N,3),':k')
plot(tt,B2(7:N,3),':k')
title('\phi_4')
subplot(223); hold on
plot(tt,Bb(7:N,4))
plot(tt,B1(7:N,4),':k')
plot(tt,B2(7:N,4),':k')
title('\beta_5')
subplot(224); hold on
plot(tt,Bb(7:N,5))
plot(tt,B1(7:N,5),':k')
plot(tt,B2(7:N,5),':k')
title('\beta_6')

figure
subplot(221); hold on
plot(tt,Bb(7:N,1))
plot(tt,B1(7:N,1),':k')
plot(tt,B2(7:N,1),':k')
title('\alpha_0')
subplot(222); hold on
plot(tt,Ee(7:N,2),'r')
plot(tt,Ee(7:N,3),'b')
plot(tt,2*sqrt(Ee(7:N,1)),':k')
plot(tt,-2*sqrt(Ee(7:N,1)),':k')
title('Residuals')
subplot(223); hold on
plot(tt,G(7:N))
plot(tt,2*sqrt(Vg(7:N)),':k')
plot(tt,-2*sqrt(Vg(7:N)),':k')
title('G gain')
subplot(224); hold on 
plot(tt,F(7:N)) 
plot(tt,4.7*ones(N-7+1,1),':k')
title('F stat')

