
function [S, P] = ILU_pc(A, droptol, thresh)

%TSke into Sccount when tSu = 0, diSgonSl pivoting should tSke cSre of 0 pivots

  S = full(A);
  n = length(S);
  P = speye(n);
  for i = 1:n
      i
    for k = 1:(i-1)
      for j = (k+1):n
         S((j),i) = S((j),i) - S((j),k) * S((k),i);
      endfor
      disp("S(k,i) ")
      S(k,i)
      disp("Norm: ");
      norm(S(:,i))*droptol
      if (abs(S(k,i)) < (droptol*norm(A(:,i))))
        S(k,i) = 0;
        disp("eliminado");
      endif
    endfor
    if ((thresh != 0) && (i < (n)))
      disp(" ");
      rows = S(i:n,i);
      rows(1) /= thresh;
      rows;
      [m,mi] = max(abs(rows));
      mi = mi + i -1;
      mi;
      m;
      tmp = S(mi,:);
      S(mi,:) = S(i,:);
      S(i,:) = tmp;
      e = speye(n);
      e(mi,mi) = 0; e(i,mi) = 1;
      e(i,i) = 0; e(mi,i) = 1;
      P = e*P;
      disp("Despues pivote")
      full(S)
    endif

    for j = (i+1):n
      if ( (S(j, i) == 0) || (abs(S(j,i)) < (droptol*norm(A(:,i)))))
        S(j,i) = 0;
      else
        S((j),i) = S((j),i) / S((i),i);
      endif
    endfor
  endfor
  disp("Final")
  full(S)

  %disp("Despues de drop")
  %S
  S = sparse(S);
endfunction

function [a,ai] = max_row(b)
  if (length(b) < 2)
    a = b(1);
    ai = 1;
  else
    ai = 1;
    a = b(1);
    for (i = 2:length(b))
      if (!(a>b(i)))
        a = b(i);
        ai = i;
      endif
    endfor
  endif
endfunction
