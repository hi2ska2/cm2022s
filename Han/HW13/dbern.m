function outdb = dbern(a, b)

q = 1.602e-19;
Kb = 1.38e-23;
T = 300;
Vt = Kb*T/q;
x = (b-a)/Vt;

if abs(x) < 0.02502
    outdb = -0.5 + x/6*(1-x^2/30*(1-x^2/28));
elseif abs(x) < 0.15
    outdb = -0.5 + x/6*(1-x^2/30*(1-x^2/28*(1-x^2/30*(1-0.03156565656565656565657*x^2))));
else
    outdb = (exp(x)-1-x*exp(x)) / (exp(x)-1)^2;
end
