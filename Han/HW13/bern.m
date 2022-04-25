function outb = bern(a, b)

q = 1.602e-19;
Kb = 1.38e-23;
T = 300;
Vt = Kb*T/q;
x = (b-a)/Vt;

if abs(x) < 0.02502
    outb = 1 - x/2 + x^2/12*(1.0-x^2/60*(1.0-x^2/42)); 
elseif abs(x) < 0.15
    outb = 1 - x/2 + x^2/12*(1.0-x^2/60*(1.0-x^2/42*(1-x^2/40*(1-0.02525252525252525252525*x^2))));
else
    outb = x / (exp(x)-1);
end
