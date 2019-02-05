function plot_Z_ok(Z_ol, Y_lk)

Z_ok=Z_ol*Y_lk;

figure,
imagesc(Z_ok);
title('imagesc(Z_ok)');

end