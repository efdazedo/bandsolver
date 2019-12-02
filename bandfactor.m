function [L,U,old2new, kl,ku] = bandfactor( A )
% [L,U, old2new, kl,ku] = bandfactor( A )
%
% input matrix A is a conceptually banded matrix
% perform LU factorization
% to estimate the lower and upper bandwidth in L and U factors
% ------------------------------------------------------------

n = size(A,1);
is_square = (size(A,1) == size(A,2));
if (~is_square),
   error(disp('bandfactor: matrix is not square %d %d', ...
	       size(A,1), size(A,2) ));
   return;
end;

% ---------
% P*A = L*U
% ---------
[L,U,P] = lu(A);

is_lower = (norm(tril(L) - L,1) == 0);
is_upper = (norm(triu(U) - U,1) == 0);
isok = (is_lower && is_upper);
if (~isok),
   error(sprintf('is_lower=%g, is_upper=%g', ...
		 is_lower,    is_upper));
   return;
end;


kl = 0;
ku = 0;
for j=1:n,
for i=1:n,
    if (L(i,j) ~= 0),
       kl = max( kl, i-j);
    end;
    if (U(i,j) ~= 0),
       ku = max( ku, j-i);
    end;
end;
end;

% ------------------------------
% precompute the explicit inverse for
% triangular blocks in L and U
% ------------------------------

% ------------------------
% note L has unit diagonal
% ------------------------
for istart=1:kl:n,
    iend = min( n, istart+kl-1);
    isize = iend - istart + 1;

    % -----------------------------------------------------------------
    % Compute  using  DTRSM
    % L( istart:iend, istart:iend) = inv( L(istart:iend,istart:iend) );
    % -----------------------------------------------------------------
    L( istart:iend, istart:iend) =  L(istart:iend,istart:iend) \ eye( isize, isize);
end;


for istart=1:ku:n,
    iend = min(n, istart+ku-1);
    isize = iend - istart + 1;
    % -----------------------------------------------------------------
    % Compute using DTRSM
    % U( istart:iend, istart:iend ) = inv( U(istart:iend, istart:iend) );
    % -----------------------------------------------------------------
    U( istart:iend, istart:iend ) = U(istart:iend, istart:iend) \ eye(isize,isize);
end;


% ---------------------------
% generate permutation vector
% ---------------------------
ip = reshape( 1:n, n,1);
old2new(1:n) = P*ip(1:n);
