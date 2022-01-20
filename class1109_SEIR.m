clear all
close all
clc
global dt initial_SEIRN

%%
dt = 0.01;
s_d = 1;
e_d = 14;
t = s_d : dt :e_d;
ttemp=s_d : 1 : e_d;
data=[3 6 25 73 222 294 258 237 191 125 69 27 11 4];


S(1)=762; 
E(1)=0;
I(1)=1;
R(1)=0;
N(1)=S(1)+E(1)+I(1)+R(1);
initial_SEIRN=[S(1) E(1) I(1) R(1) N(1)];

ini=[1 0.5 0.5];
ub=[10 1 1];
lb=[0 0 0];

y=lsqcurvefit(@ftSEIR,ini,t,data,lb,ub,[]);
beta=y(1);
kappa=y(2);
alpha=y(3);


for i=1:length(t)-1
    S(i+1)=S(i)+dt*(-beta*S(i)*I(i)/N(i));
    E(i+1)=E(i)+dt*(beta*S(i)*I(i)/N(i) - kappa*E(i));
    I(i+1)=I(i)+dt*(kappa*E(i)-alpha*I(i));
    R(i+1)=R(i)+dt*alpha*I(i);
    N(i+1)=S(i+1)+I(i+1)+R(i+1)+E(i+1);
end


figure(1)
set(gcf,'color','w');
plot(t,S,'b-','linewidth',2)
hold on
plot(t,E,'r-','linewidth',2)
hold on
plot(t,I,'k-','linewidth',2)
hold on
plot(t,R,'g-','linewidth',2)
hold on
plot(ttemp,data,'ko','markersize',5)
legend('S','E','I','R')
xlim([s_d,e_d])
title( sprintf('\\beta=%.4g ,\\kappa=%.4g,\\alpha=%.4g',beta,kappa,alpha) )
