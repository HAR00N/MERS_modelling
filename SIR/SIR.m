clear all
close all
clc

global dt initial_SIRN


%%
%time array
dt=0.01;
ini_data=1;
end_data=14;
time=ini_data:dt:end_data;
ttemp=ini_data:1:end_data;
data=[3 6 25 73 222 294 258 237 191 125 69 27 11 4];


S(1)=762; 
I(1)=1; 
R(1)=0; 
N(1)=S(1)+I(1)+R(1);
initial_SIRN=[S(1) I(1) R(1) N(1)];

%fitting
ini=[1 0.5];
ub=[10 1];
lb=[0 0];
y=lsqcurvefit(@ftSIR,ini,time,data,lb,ub);
beta=y(1); alpha=y(2);

for i=1:length(time)-1
    S(i+1)=S(i)+dt*(-beta*S(i)*I(i)/N(i));
    I(i+1)=I(i)+dt*(beta*S(i)*I(i)/N(i)-alpha*I(i));
    R(i+1)=R(i)+dt*alpha*I(i);
    N(i+1)=S(i+1)+I(i+1)+R(i+1);
end

%plot
figure(1)
set(gcf,'color','w');
plot(time,S,'b-','linewidth',2)
hold on
plot(time,I,'k-','linewidth',2)
plot(time,R,'r-','linewidth',2)
plot(ttemp,data,'ko','markersize',5)
hold off
legend('S','I','R')
xlim([ini_data,end_data])
title(sprintf('\\beta=%.4g ,\\alpha=%.4g',beta,alpha))
