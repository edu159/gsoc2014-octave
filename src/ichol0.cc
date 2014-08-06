/**
 * Copyright (C) 2014 Eduardo Ramos Fernández <eduradical951@gmail.com>
 *
 * This file is part of Octave.
 *
 * Octave is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * Octave is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Octave; see the file COPYING.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <octave/oct.h>
#include <octave/parse.h>
#include <time.h>

/**
 * Performs the multiplication needed by the Cholesky factorization in the
 * complex case.
 */
template < typename T > inline T
ichol_mult_complex (T a, T b)
{
  b.imag (-std::imag (b));
  return a * b;
}


template < typename T > inline bool
ichol_checkpivot_complex (T pivot)
{
  if (pivot.imag () != 0)
    {
      error ("Non-real pivot encountered. Input matrix must be hermitian.");
      return false;
    }
  else if (pivot.real () < 0)
    {
      error ("Non-positive pivot encountered.");
      return false;
    }
  return true;

}

template < typename T > inline bool
ichol_checkpivot_real (T pivot)
{
  if (pivot < T(0))
    {
      error ("Non-positive pivot encountered.");
      return false;
    }
  return true;

}

/**
 * Performs the multiplication needed by the Cholesky factorization in the
 * real case.
 */
template < typename T> inline T 
ichol_mult_real (T a, T b)
{
  return a * b;
}

template <typename octave_matrix_t, typename T, T (*ichol_mult) (T, T), bool (*ichol_checkpivot) (T)>
void ichol_0 (octave_matrix_t& sm, const std::string michol = "off") 
{

  const octave_idx_type n = sm.cols ();
  OCTAVE_LOCAL_BUFFER (octave_idx_type, iw, n);
  octave_idx_type j1, jend, j2, jrow, jjrow, jw, i, k, jj, Llist_len_aux, Llist_len, r;
  T tl;

  char opt;
  enum {OFF, ON};
  if (michol == "on")
    opt = ON;
  else
    opt = OFF;

clock_t t_i, t_f, t_ii, t_ff; 
double time_spent1, time_spent2;

  octave_idx_type* cidx = sm.cidx ();
  octave_idx_type* ridx = sm.ridx ();
  T* data = sm.data ();
  OCTAVE_LOCAL_BUFFER(octave_idx_type, Lfirst, n);
  OCTAVE_LOCAL_BUFFER(octave_idx_type, Llist, n);
  OCTAVE_LOCAL_BUFFER(T, dropsums, n);

  for (i = 0; i < n; i++)
    {
      iw[i] = -1;
      Llist[i] = -1;
      Lfirst[i] = -1;
      dropsums[i] = 0;
    }
  time_spent1 = 0;
  time_spent2 = 0;
  for (k = 0; k < n; k++)
    {
      //printf("k: %d\n", k);
      j1 = cidx[k];
      j2 = cidx[k+1];
      octave_idx_type j;
      for (j = j1; j < j2; j++)
        iw[ridx[j]] = j;
      jrow = Llist [k];
//      t_i = clock ();
      while (jrow != -1) 
        {
          jjrow = Lfirst[jrow];
          jend = cidx[jrow+1];
          for (jj = jjrow; jj < jend; jj++)
            {
              r = ridx[jj];
              jw = iw[r];
              tl = ichol_mult (data[jj], data[jjrow]);
              if (jw != -1)
                data[jw] -= tl;
              else
                if (opt == ON)
                {
                  dropsums[r] -= tl;
                  dropsums[k] -= tl;
                }
            }
          if ((jjrow + 1) < jend)
            {
              Lfirst[jrow]++;
              j = jrow;
              jrow = Llist[jrow];
              Llist[j] = Llist[ridx[Lfirst[j]]];
              Llist[ridx[Lfirst[j]]] = j;
            }
          else
            jrow = Llist[jrow];
        }
      /**
      printf("Llist_in: ");
      for (j = 0; j < n; j++)
        printf("%d ", Llist[j]);
      printf("\n");
      printf("Lfirst_in: ");
      for (j = 0; j < n; j++)
        printf("%d ", Lfirst[j]);
      printf("\n");
      **/
  //    t_f = clock ();
     // time_spent1 += (double) (t_f - t_i) /CLOCKS_PER_SEC;
    //  t_ii = clock ();
      if (opt == ON)
        data[j1] += dropsums[k];

      if (ridx[j1] != k)
        {
          error ("ilu0: There is a pivot equal to zero.");
          break;
        }
      if (!ichol_checkpivot (data[j1]))
        break;

      data[cidx[k]] = std::sqrt (data[j1]);

              
      if (k < (n - 1)) 
        {
          iw[ridx[j1]] = -1;
          for(i = j1 + 1; i < j2; i++)
            {
              iw[ridx[i]] = -1;
              data[i] /=  data[j1];
            }
          Lfirst[k] = j1;
          if ((Lfirst[k] + 1) < j2)
            {
              Lfirst[k]++;
              jjrow = ridx[Lfirst[k]];
              Llist[k] = Llist[jjrow];
              Llist[jjrow] = k;
            }
        }

      //printf("Llist_len: %d\n", Llist_len);
      /**
      printf("Llist_out: ");
      for (j = 0; j < n; j++)
        printf("%d ", Llist[j]);
      printf("\n");
      printf("Lfirst_out: ");
      for (j = 0; j < n; j++)
        printf("%d ", Lfirst[j]);
      printf("\n");
      **/
       // t_ff = clock ();
        //time_spent2 += (double) (t_ff - t_ii) /CLOCKS_PER_SEC;
    }
//  printf ("Time1: %f\n", time_spent1);
 //printf ("Time2: %f\n", time_spent2);
}

DEFUN_DLD (ichol0, args, nargout, "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{L} =} ichol0 (@var{A}, @var{michol})\n\
\n\
Computes the no fill Incomplete Cholesky factorization [IC(0)] of A \
which must be an square hermitian matrix in the complex case and a symmetric \
positive definite matrix in the real one. \
\n\
\n\
@code{[@var{L} = ichol0 (@var{A}, @var{michol})} \
computes the IC(0) of @var{A}, such that @code{@var{L} * @var{L}'} which \
is an approximation of the square hermitian matrix @var{A}. \
The parameter @var{michol} decides whether the Modified IC(0) should \
be performed. This compensates the main diagonal of \
@var{L}, such that @code{@var{A} * @var{e} = @var{L} * @var{L}' * @var{e}} \
with @code{@var{e} = ones (size (@var{A}, 2), 1))} holds. \n\
\n\
For more information about the algorithms themselves see:\n\
\n\
[1] Saad, Yousef. \"Preconditioning Techniques.\" Iterative Methods for Sparse Linear \
Systems. PWS Publishing Company, 1996. \
\n\
@seealso{ichol, icholt, chol, ilu}\n\
@end deftypefn")

{
  octave_value_list retval;

  int nargin = args.length ();
  std::string michol;
 

  if (nargout > 1 || nargin < 1 || nargin > 2)
    {
      print_usage ();
      return retval;
    }

  if (args (0).is_scalar_type () || !args (0).is_sparse_type ())
    error ("ichol0: 1. parameter must be a sparse square matrix.");

  if (args (0).is_empty ())
    {
      retval (0) = octave_value (SparseMatrix ());
      return retval;
    }


  if (nargin == 2)
    {
      michol = args (1).string_value ();
      if (error_state || !(michol == "on" || michol == "off"))
        error ("ichol0: 2. parameter must be 'on' or 'off' character string.");
      // maybe resolve michol to a numerical value / enum type already here!
    }


  if (!error_state)
    {
      // In ICHOL0 algorithm the zero-pattern of the input matrix is preserved so
      // it's structure does not change during the algorithm. The same input
      // matrix is used to build the output matrix due to that fact.
      octave_value_list param_list;
      if (!args (0).is_complex_type ())
        {
          SparseMatrix sm = args (0).sparse_matrix_value ();
          param_list.append (sm);
          sm = feval ("tril", param_list)(0).sparse_matrix_value (); 
          ichol_0 <SparseMatrix, double, ichol_mult_real, ichol_checkpivot_real> (sm, michol);
          if (!error_state)
            retval (0) = octave_value (sm);
        }
      else
        {
          SparseComplexMatrix sm = args (0).sparse_complex_matrix_value ();
          param_list.append (sm);
          sm = feval ("tril", param_list)(0).sparse_complex_matrix_value (); 
          ichol_0 <SparseComplexMatrix, Complex, ichol_mult_complex, ichol_checkpivot_complex> (sm, michol);
          if (!error_state)
            retval (0) = octave_value (sm);
        }

    }

  return retval;
}

/*
%!shared A_1, A_off_1, A_2, A_off_2, A_3, A_off_3, A_4, A_off_4, A_5, A_off_5
%! A_1 = [ 0.37, -0.05,  -0.05,  -0.07;
%!        -0.05,  0.116,  0.0,   -0.05;
%!        -0.05,  0.0,    0.116, -0.05;
%!        -0.07, -0.05,  -0.05,   0.202];
%! A_1 = sparse(A_1);
%! A_off_1 = (tril (A_1));
%!
%! A_2 = gallery ('poisson', 30);
%! A_off_2 = (tril (A_2));
%!
%! A_3 = gallery ('tridiag', 50);
%! A_off_3 = (tril (A_3));
%!
%! nx = 400; ny = 200;
%! hx = 1 / (nx + 1); hy = 1 / (ny + 1);
%! Dxx = spdiags ([ones(nx, 1), -2 * ones(nx, 1), ones(nx, 1)], [-1 0 1 ], nx, nx) / (hx ^ 2);
%! Dyy = spdiags ([ones(ny, 1), -2 * ones(ny, 1), ones(ny, 1)], [-1 0 1 ], ny, ny) / (hy ^ 2);
%! A_4 = -kron (Dxx, speye (ny)) - kron (speye (nx), Dyy);
%! A_off_4 = (tril (A_4));
%!
%! A_5 = [ 0.37, -0.05,          -0.05,  -0.07;
%!        -0.05,  0.116,          0.0,   -0.05 + 0.05i;
%!        -0.05,  0.0,            0.116, -0.05;
%!        -0.07, -0.05 - 0.05i,  -0.05,   0.202];
%! A_5 = sparse(A_5);
%! A_off_5 = (tril (A_5));
%!
%% Test input
%!test
%!error ichol0 ([]);
%!error ichol0 ([],[]);
%!error ichol0 ([],[],[]);
%!error [~,~] = ichol0 ([],[],[]);
%!error [L] = ichol0 ([], 'foo');
%!error [L] = ichol0 (A_off_1, [], false);
%!error [L, E] = ichol0 (A_off_1, false);
%!error ichol0 (sparse (0), false);
%!
%!test
%! L = ichol0 (sparse (1), 'off');
%! assert (L, sparse (1));
%! L = ichol0 (sparse (2), 'off');
%! assert (L, sparse (sqrt (2)));
%!
%!test
%! L = ichol0 (A_off_1, 'off');
%! assert (norm (A_1 - L*L', 'fro') / norm (A_1, 'fro'), 1e-2, 1e-2);
%! L = ichol0 (A_off_1, 'on');
%! assert (norm (A_1 - L*L', 'fro') / norm (A_1, 'fro'), 2e-2, 1e-2);
%!
%!test
%! L = ichol0 (A_off_2, 'off');
%! assert (norm (A_2 - L*L', 'fro') / norm (A_2, 'fro'), 1e-1, 1e-1)
%!
%!test
%! L = ichol0 (A_off_3, 'off');
%! assert (norm (A_3 - L*L', 'fro') / norm (A_3, 'fro'), eps, eps);
%! L = ichol0 (A_off_3, 'on');
%! assert (norm (A_3 - L*L', 'fro') / norm (A_3, 'fro'), eps, eps);
%!
%!test
%! L = ichol0 (A_off_4, 'off');
%! assert (norm (A_4 - L*L', 'fro') / norm (A_4, 'fro'), 1e-1, 1e-1);
%!
%!test
%! L = ichol0 (A_off_5, 'off');
%! assert (norm (A_5 - L*L', 'fro') / norm (A_5, 'fro'), 1e-2, 1e-2);
%! L = ichol0 (A_off_5, 'on');
%! assert (norm (A_5 - L*L', 'fro') / norm (A_5, 'fro'), 2e-2, 1e-2);
*/


