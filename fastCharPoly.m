function p = fastCharPoly( A )
% Fast algorithm for characteristic polynomial of square matrices.
% Algortihm is described in 
% La Budde's Method For Computing Characteristic Polynomials
% by Rizwana Rehman and Ilse C.F. Ipsen
% 
% Input: square matrix A
% 
% Output: charactersitic polynomial p of A
%
% Author: Sebastian J. Schlecht
% Date: 11.04.2015


N = size(A, 1);
H = hess(A);
beta = diag(H(2:end,1:end-1));
alpha = diag(H);

pH = zeros(N+1,N+1);
pH(0+1, 1) = 1;
pH(1+1, 1:2) = [-alpha(1), 1];

for it = 2:N
    partB = flipud(cumprod(flipud(beta(1 : it-1))));
    partH = (H(1:it-1,it));
    partP = pH(1:it-1,:);
    
    rec = sum(bsxfun(@times, partB.*partH, partP),1);
    convp = conv(pH(it-1+1, :), [-alpha(it), 1]);
    
    pH(it+1, :) = convp(1:end-1) - rec;
end

p = fliplr(pH(end, :));