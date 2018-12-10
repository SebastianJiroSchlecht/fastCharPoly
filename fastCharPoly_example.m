% fast charpoly LaBudde
clear; clc;

N = 12;
A = randn(N) + 1i*randn(N);

disp('Comparison of Matlab charpoly and LaBudde charpoly')

%%-----------------
disp('Matlab')
tic
pMatlab = charpoly(A);
toc
disp(pMatlab)

%%-----------------
disp('LaBudde')
tic
pFast = fastCharPoly( A );
toc
disp(pFast)