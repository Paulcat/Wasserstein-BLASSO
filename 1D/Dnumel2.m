function C = Dnumel2(m)
%DNUMEL(N) Number of 'diagonal' elements for 2-level Toeplitz matrix

O = ones(m);
O = padarray(O,m-1,'post');
C = ifftshift( ifft2(fft2(O).^2) );

end

