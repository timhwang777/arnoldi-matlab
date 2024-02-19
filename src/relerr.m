relar = [1.2642e-10
1.4205e-15
7.7958e-16
3.1532e-16
1.2051e-15
2.6948e-16
6.0164e-16
6.6849e-17
5.102e-16
5.5125e-16];

relor = [1.2642e-10
8.0496e-16
5.7115e-16
4.7062e-16
7.9729e-16
3.0082e-16
1.4948e-16
1.3781e-16
4.4843e-16
2.9896e-16];

figure(7)
X = 10:10:100;
semilogy(X, relar, '-o',X, relor, '-*')
legend('Standard','re-Orthogonal','Location','northwest')
grid on