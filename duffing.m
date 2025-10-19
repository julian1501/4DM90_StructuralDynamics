function f = duffing(t,x,alpha,beta,gamma,delta,f_e)
% Duffing system 
q = x(1);
v = x(2);
f = [v ; 
     -delta*v-alpha*q-beta*q^3+gamma*cos(2*pi*f_e*t)];


