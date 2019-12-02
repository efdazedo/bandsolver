function [ x ] = bandsolve( n,kl,ku, L,U,old2new,  b )
%
% [ x ] = bandsolve( n,kl,ku, L,U,old2new,  b )
% -------------------------------
% L has lower bandwidth kl
% U has upper bandwith ku
% diagonal blocks of L has explicit inverse
% diagonal blocks of U has explicit inverse
% -------------------------------
idebug = 1;

n = size(L,1);
x = zeros(n,1);

if (idebug >= 1),
  [ilist,jlist,alist] = find( L );
  is_lower = (norm( L - tril(L),1) == 0);
  max_kl = max(ilist-jlist);
  is_kl_ok = (max_kl == kl);
  is_L_ok = (is_lower && is_kl_ok);
  if (~is_L_ok),
    error(sprintf('bandsolve: kl=%d, max_kl=%d, is_lower=%d', ...
		              kl,    max_kl,    is_lower));
    return;
  end;

  [ilist,jlist,alist] = find( U );
  is_upper = (norm( U - triu(U),1) == 0);
  max_ku = max( jlist - ilist);
  is_ku_ok = (max_ku == ku);
  is_U_ok = (is_upper && is_ku_ok);
  if (~is_U_ok),
    error(sprintf('bandsolve: ku=%d, max_ku=%d, is_upper=%d', ...
		              ku,    max_ku,    is_upper));
    return;
  end;

  clear ilist;
  clear jlist;
  clear alist;
end;


% --------------------
% solve  A*x = b,
% P*A = L * U
% P*A*x = P*b
% L*U*x = (P*b)
% L*y = (P*b), y = U*x
% --------------------

% -------------------
% perform permutation
% -------------------
x(1:n) = b(old2new(1:n));

% --------------
% solve L*y = P*b
% --------------
for istart=1:kl:n,
    iend = min(n,istart+kl-1);
    x(istart:iend) =  tril(L( istart:iend, istart:iend) )* x(istart:iend);
    i1 = iend + 1;
    i2 = min(n, i1 + kl-1);
    isize = i2-i1+1;
    if (isize >= 1),
      x( i1:i2) = x(i1:i2) - triu(L( i1:i2, istart:iend)) * x(istart:iend);
    end;
end;


% -------------
% solve y = U*x
% -------------

max_k = ceil(n/ku);

for k=max_k:-1:1,
    istart = 1 + (k-1)*ku;
    iend = min(n, istart+ku-1);
    x( istart:iend ) = triu(U( istart:iend, istart:iend) ) * x( istart:iend);

    i2 = istart -1 ;
    i1 = max(1, i2 - ku + 1);
    isize = i2-i1 + 1;
    if (isize >= 1),
       x(i1:i2) = x(i1:i2) - tril(U( i1:i2, istart:iend))*x(istart:iend);
    end;
end;
