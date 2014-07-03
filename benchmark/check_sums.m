function check_sums (A, L, U, milu)
  b = L * U;
  dim = 1;
  if (milu == 'row')
    dim = 2;
  endif
  c = sum (b, dim);
  d = sum (A, dim);
  v = abs (c - d);
  num_zeros = length (find (v > sqrt (eps)));
  printf('Number of rows-columns not well compensated: %d\n', num_zeros);
  if (num_zeros > 0)
    v (find (v > sqrt (eps)))
  endif
endfunction
