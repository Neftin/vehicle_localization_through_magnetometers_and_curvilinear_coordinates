

t = 1; rx = 3; ry = 3; rz = 7; mx = -200; mz = -200; mx = 100;

mu0 = 10^-7;

cg5 = [mu0 * (rx(t) ^ 2 + ry(t) ^ 2 + rz(t) ^ 2) ^ (-0.5e1 / 0.2e1) * (3 * (rx(t) * mx(t) + ry(t) * my(t) + rz(t) * mz(t)) * mx(t) + (-rx(t) ^ 2 - ry(t) ^ 2 - rz(t) ^ 2) * mx(t)) / pi / 0.4e1 mu0 * (rx(t) ^ 2 + ry(t) ^ 2 + rz(t) ^ 2) ^ (-0.5e1 / 0.2e1) * (3 * (rx(t) * mx(t) + ry(t) * my(t) + rz(t) * mz(t)) * my(t) + (-rx(t) ^ 2 - ry(t) ^ 2 - rz(t) ^ 2) * my(t)) / pi / 0.4e1 mu0 * (rx(t) ^ 2 + ry(t) ^ 2 + rz(t) ^ 2) ^ (-0.5e1 / 0.2e1) * (3 * (rx(t) * mx(t) + ry(t) * my(t) + rz(t) * mz(t)) * mz(t) + (-rx(t) ^ 2 - ry(t) ^ 2 - rz(t) ^ 2) * mz(t)) / pi / 0.4e1];

cg5

3*(rx(t)*mx(t)+ry(t)*my(t)+rz(t)*mz(t)) *rx(t)

(-rx(t)^2-ry(t)^2-rz(t)^2)*mx(t)