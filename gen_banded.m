function [A] = gen_banded(n,kl,ku)
% [A] = gen_banded(n,kl,ku)
% ----------------------
% generate banded matrix of size n by n
% with lower bandwidth kl and upper bandwidth ku
% ----------------------
A = rand(n,n)*2-1;
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
