function hd = ideal_lp(wc,M);

alpha = (M-1)/2;
n = [-alpha:1:alpha];
m = n + eps;
hd = sin(wc*m) ./ (pi*m);

