function [U,nit] = ffw_bfgs(U,f,g,options)
%FFW_BFGS BFGS step in FFW algorithm

deal2 = @(varargin) deal(varargin{1:nargout});

r = size(U,2);

% from complex to real and 'vice-versa'
reshc = @(u) reshape( u(1:length(u)/2)    , [length(u)/2/r, r] ) + ...
        1i * reshape( u(length(u)/2+1:end), [length(u)/2/r, r] );
%
flatc = @(Z) [real(Z(:)); imag(Z(:))];

fb = @(Z) f(reshc(Z));
gb = @(Z) flatc( g(reshc(Z)) );
fg = @(Z) deal2(fb(Z),gb(Z));

%checkgradient(fb,gb,flatc(U));
%drawnow;
%pause(1);

[U,val,exitflag,output] = lbfgs(fg,flatc(U),options);
U = reshc(U);

nit = output.iterations;
end

