function [A] = gen_banded(n,kl,ku)
% [A] = gen_banded(n,kl,ku)
% ----------------------
% generate banded matrix of size n by n
% with lower bandwidth kl and upper bandwidth ku
% ----------------------
idebug = 2;
A = rand(n,n)*2-1;
if (idebug >= 2),
    A = -reshape( 1:(n*n), n,n);
end;

for j=1:n,
for i=1:n,
  if (i-j > kl),
    A(i,j) = 0;
  end;
  if (j-i > ku),
    A(i,j) = 0;
  end;
end;
end;

is_diag_dominant = 0;
if (is_diag_dominant),
   for j=1:n,
       i = j;
       A(i,i) = max( 10^4, 2*n*n);
   end;
end;
