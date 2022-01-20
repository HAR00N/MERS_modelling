function y=ftSEIR(x,t)
global dt initial_SEIRN

beta=x(1);
kappa=x(2);
alpha=x(3);

S(1)=initial_SEIRN(1);
E(1)=initial_SEIRN(2);
I(1)=initial_SEIRN(3);
R(1)=initial_SEIRN(4);
N(1)=initial_SEIRN(5);

for i=1:length(t)-1
    S(i+1)=S(i)+dt*(-beta*S(i)*I(i)/N(i));
    E(i+1)=E(i)+dt*(beta*S(i)*I(i)/N(i) - kappa*E(i));
    I(i+1)=I(i)+dt*(kappa*E(i)-alpha*I(i));
    R(i+1)=R(i)+dt*alpha*I(i);
    N(i+1)=S(i+1)+I(i+1)+R(i+1)+E(i+1);
end
y=I(1:1/dt:end);
end