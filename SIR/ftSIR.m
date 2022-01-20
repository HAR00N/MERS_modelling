function y=ftemp(x,time)
global dt initial_SIRN
beta=x(1);
alpha=x(2);

S(1)=initial_SIRN(1);
I(1)=initial_SIRN(2);
R(1)=initial_SIRN(3);
N(1)=initial_SIRN(4);

for i=1:length(time)-1
    S(i+1)=S(i)+dt*(-beta*S(i)*I(i)/N(i));
    I(i+1)=I(i)+dt*(beta*S(i)*I(i)/N(i)-alpha*I(i));
    R(i+1)=R(i)+dt*alpha*I(i);
    N(i+1)=S(i+1)+I(i+1)+R(i+1);
end
y=I(1:1/dt:end);
length(y)
end