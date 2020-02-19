% method 3
function t = time_map3(N, read_out)
delta = 1/N;
x = -0.5+delta/2:delta:0.5-delta/2;
[X,Y] = meshgrid(x,x);
t = zeros(N,N);
for i = 1:N
    for j = 1:N
        r2 = X(i,j)^2 + Y(i,j)^2;
        if r2 > 0.5^2
            t(i,j) = 0; % Value outside the circle
        else
            t(i,j) = 4*read_out*r2;
        end
    end
end
t = t;
end