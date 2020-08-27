function [Tx] = Tprod2(m,Tvals,x)
%TPROD2 2-level Toeplitz product
%   TPROD2(N,TVALS,X) returns the product between the 2-level Toeplitz
%   matrix whose 'diagonal' values are given by TVALS, and the array X.
%
%   N = (N1,N2) is the dimensions of the 2-level Toeplitz matrix
%
%   TVALS is a matrix of size (2*N1-1,2*N2-1)
%   TODO: explain order
%
%   X is a matrix of size (N1*N2,r)
%
%   Internal use only

r = size(x,2);
x = reshape(x,[m,r]);

Px = padarray(x,[m-1, 0],'post');
Cx = ifft2( fft2(Tvals) .* fft2(Px) ); % broadcasting...

Tx = Cx(1:m(1),1:m(2),:);
Tx = reshape(Tx,prod(m),r);


end

