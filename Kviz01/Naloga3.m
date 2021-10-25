1;
format long;

f = @(x) tan(x) - x;

sedmaPozNicla = fsolve(f,23.5)