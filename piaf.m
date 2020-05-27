points = dlmread('points.mat', ' ');
tt = points(:, 1);
xx = points(:, 2);
yy = points(:, 3);

sx = spline(tt, xx);
sy = spline(tt, yy);

figure;
period = tt(length(tt)) - tt(1);

t = linspace(tt(1), tt(length(tt)), 1000);
x = ppval(sx, t);
y = ppval(sy, t);
plot(x, y);