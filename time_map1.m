% method 1
function t = time_map1(N, k_points, interval)
[max_i, c] = size(k_points);
x = real(k_points);
x = reshape(x,[],1);
y = imag(k_points);
y  = reshape(y,[],1);
index = 1:max_i;
z = index';
z = repmat(z,c,1);
F = scatteredInterpolant(x,y,z);
delta = 1/N;
xx = -0.5+delta/2:delta:0.5-delta/2;
[X,Y] = meshgrid(xx,xx);
% F.Method = 'nearest';
F.Method = 'linear';
T = F(X,Y);
T(T>=max_i) = 0;
t = interval*T;
end