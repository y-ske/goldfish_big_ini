function dx = my_bibun(x,fs)
c = fs/2;
N = length(x);
dx = x*0;
for k=2:N-1
dx(k) = c*(x(k+1) - x(k-1));
end
dx(1) = (-x(1)+x(2))*fs;
dx(N) = (-x(N-1) + x(N))*fs;
end
