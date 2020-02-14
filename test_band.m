function [err,res] = test_band( n, kl, ku )
% [err, res] = test_band( n, kl, ku )
% 
% simple test for bandfactor() and bandsolver()
% ---------------------------------------------
idebug = 1;

A = gen_banded(n,kl,ku);

use_sparse = 0;
if (use_sparse), 
        [ii,jj,aij] = find(A);
        A = sparse(ii,jj,aij, n,n);
end;

% -------------------------
% generate solution and rhs
% -------------------------
x = 2*rand(n,1)-1;
if (idebug >= 1),
   x = reshape(1:n,n,1);
end;
b = A * x;

% ---------------------
% perform factorization
%
% note new bandwidth may be larger due to pivoting
% kl2 ~ 2*(kl+ku), ku2 ~ 2*ku
% ---------------------
[L,U,old2new,kl2,ku2] = bandfactor(A);

[x2]  = bandsolve(n,kl2,ku2, L,U,old2new,  b);

res = norm( b - A*x2 );
err = norm( x - x2 );

if (idebug >= 1),
   disp(sprintf('test_band:n=%d, kl=%d ,ku=%d, kl2=%d, ku2=%d', ...
		 n,    kl,    ku,    kl2,    ku2 ));
   disp(sprintf('err=%g, res=%g', err, res));
end;


