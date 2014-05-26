function [A, P] = ILU_pc(A, tau, thresh)

%Take into account when tau = 0, diagonal pivoting should take care of 0 pivots

B = A;
n = length(A);
P = speye(n);
for i = 1:n
  for k = i:n
    A(k:n,k) *= thresh;
    A(k,k) /= thresh;
    [m,mi] = max(abs(A(k:n,k)))
    A(k,k) *= thresh;
    A(k:n,k) /= thresh;
    mi = mi + k -1;
    tmp = A(mi,:);
    A(mi,:) = A(k,:);
    A(k,:) = tmp;
    e = speye(n);
    e(mi,mi) = 0; e(k,mi) = 1;
    e(k,k) = 0; e(mi,k) = 1;
    P = e*P;
  endfor
  for k = 1:i-1
       if ( (A(i,k) == 0) || (abs(A(i,k)) < (tau*norm(B(:,k)))))
          A(i,k) = 0;
          continue
       endif
       A(i,k) = A(i,k) / A(k,k);
       A(i,k+1:n) = A(i,k+1:n) - A(i,k) * A(k,k+1:n);
  endfor
endfor

for i = 1:n
  for j = i+1:n
    if (abs(A(i,j)) < (tau*norm(B(:,j))))
      A(i,j) = 0;
    end
  end
end

