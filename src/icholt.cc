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
      error ("Non-real pivot encountered.");
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

/* 
 * That function implements the IKJ and JKI variants of gaussian elimination to
 * perform the ILUTP decomposition. The behaviour is controlled by michol
 * parameter. If michol = ['off'|'col'] the JKI version is performed taking
 * advantage of CCS format of the input matrix. If michol = 'row' the input matrix
 * has to be transposed to obtain the equivalent CRS structure so we can work
 * efficiently with rows. In this case IKJ version is used.
 */

template <typename octave_matrix_t, typename T,  T(*ichol_mult) (T, T), 
          bool(*ichol_checkpivot) (T)>
void ichol_t (const octave_matrix_t& sm, octave_matrix_t& L, const T* cols_norm,
              const T droptol, const std::string michol = "off")
              
{

  const octave_idx_type n = sm.cols ();
  const octave_idx_type nnz = sm.nnz ();
  OCTAVE_LOCAL_BUFFER (octave_idx_type, iw, n);
  octave_idx_type j, jrow,jjrow, jw, i, k, jj, Llist_len, total_len, w_len,
                  max_len;

  char opt;
  enum {OFF, ON};
  if (michol == "on")
    opt = ON;
  else
    opt = OFF;

//clock_t t_i, t_f, t_ii, t_ff; 
//double time_spent1, time_spent2;

  octave_idx_type* cidx_in = sm.cidx ();
  octave_idx_type* ridx_in = sm.ridx ();
  T* data_in = sm.data ();
  OCTAVE_LOCAL_BUFFER(T, w_data, n);
  OCTAVE_LOCAL_BUFFER(octave_idx_type, Lfirst, n);
  OCTAVE_LOCAL_BUFFER(octave_idx_type, Llist, n);
  OCTAVE_LOCAL_BUFFER(T, col_drops, n);
  OCTAVE_LOCAL_BUFFER(T, col_drops_future, n);
  // Arrays for L
  max_len = nnz;
  max_len += (0.1 * max_len) > n ? 0.1 * max_len : n;
  Array <octave_idx_type> cidx_out_l (dim_vector (n + 1,1));
  octave_idx_type* cidx_l = cidx_out_l.fortran_vec ();
  Array <octave_idx_type> ridx_out_l (dim_vector (max_len ,1));
  octave_idx_type* ridx_l = ridx_out_l.fortran_vec ();
  Array <T> data_out_l (dim_vector (max_len, 1));
  T* data_l = data_out_l.fortran_vec ();

  T zero = T(0);

  cidx_l[0] = cidx_in[0];
  for (i = 0; i < n; i++)
    {
      Llist[i] = 0;
      Lfirst[i] = 0;
      w_data[i] = 0;
      col_drops[i] = zero;
      col_drops_future[i] = zero;
    }
  Llist_len = 0;
  total_len = 0;
  //time_spent1 = 0;
  //time_spent2 = 0;
  for (k = 0; k < n; k++)
    {
    //  printf("k: %d\n", k);
      for (j = cidx_in[k]; j < cidx_in[k+1]; j++)
        w_data[ridx_in[j]] = data_in[j];
      j = 0;
      jrow = Llist[j];
      while (j < Llist_len) 
        {
          jjrow = Lfirst[jrow];
          for (jj = jjrow; jj < cidx_l[jrow+1]; jj++)
            w_data[ridx_l[jj]] -=  ichol_mult (data_l[jj], data_l[jjrow]);
          j++;
          jrow = Llist[j];
        }

      //t_i = clock ();
      if ((max_len - total_len) < n)
        {
          max_len += (0.1 * max_len) > n ? 0.1 * max_len : n;
          data_out_l.resize (dim_vector (max_len, 1));
          data_l = data_out_l.fortran_vec ();
          ridx_out_l.resize (dim_vector (max_len, 1));
          ridx_l = ridx_out_l.fortran_vec ();
        }
    //  t_f = clock ();
     // time_spent1 += (double) (t_f - t_i) /CLOCKS_PER_SEC;
      data_l[total_len] = w_data[k];
      ridx_l[total_len] = k;
      w_len = 1;
      for (i = k + 1; i < n; i++)
        {
          if (w_data[i] != zero)
            {
              if (std::abs (w_data[i]) < (droptol * cols_norm[k]))
                {
                  if (opt == ON)
                    {
                      col_drops[k] += w_data[i];
                      col_drops_future[i] += w_data[i];
                    }
                }
              else
                {
                  data_l[total_len + w_len] = w_data[i];
                  ridx_l[total_len + w_len] = i;
                  w_len++;
                }
            }
          w_data[i] = zero;
        }

      // Compensate column sums --> michol option
      if (opt == ON)
        data_l[total_len] += col_drops[k] + col_drops_future[k];

      if (data_l[total_len] == zero)
        {
          error ("ilu0: There is a pivot equal to zero.");
          break;
        }
      else if (!ichol_checkpivot (data_l[total_len]))
        break;
      data_l[total_len] = std::sqrt(data_l[total_len]);
      for (jj = total_len + 1; jj < (total_len + w_len); jj++)
        data_l[jj] /=  data_l[total_len];
      total_len += w_len;
      cidx_l[k+1] = cidx_l[k] - cidx_l[0] + w_len;
     // t_ii = clock ();
      if (k < (n - 1))
        {
          Llist_len = 0;
          Lfirst[k] = cidx_l[k];
          for (i = 0; i <= k; i++)
            {
              jj = ridx_l[Lfirst[i]];
              if (jj < (k + 1))
                if(Lfirst[i] < (cidx_l[i+1]))
                  {
                    if ((Lfirst[i] + 1) < cidx_l[i+1])
                      {
                        Lfirst[i]++;
                        jj = ridx_l[Lfirst[i]];
                      }
                  }
              if (jj == (k + 1)) 
                {
                  Llist[Llist_len] = i;
                  Llist_len++;
                }
            }
        }
     // t_ff = clock ();
      //time_spent2 += (double) (t_ff - t_ii) /CLOCKS_PER_SEC;
    }
     // printf ("Time1: %f\n", time_spent1);
     // printf ("Time2: %f\n", time_spent2);
  if (!error_state)
    {
      // Build the output matrices
      L = octave_matrix_t (n, n, total_len);
      for (i = 0; i <= n; i++)
        L.cidx (i) = cidx_l[i];
      for (i = 0; i < total_len; i++)
        {
          L.ridx (i) = ridx_l[i];
          L.data (i) = data_l[i];
        }
    }

}

DEFUN_DLD (icholt, args, nargout, "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {[@var{L}, @var{U}] =} ilu0 (@var{A})\n\
@deftypefnx  {Loadable Function} {[@var{L}, @var{U}] =} ilu0 (@var{A}, @var{michol})\n\
\n\
NOTE: No pivoting is performed.\n\
\n\
Computes the incomplete LU-factorization (ILU) with 0-order level of fill of \
@var{A}.\n\
\n\
@code{[@var{L}, @var{U}] = ilu0 (@var{A})} computes the zero fill-in ILU-\
factorization ILU(0) of @var{A}, such that @code{@var{L} * @var{U}} is an \
approximation of the square sparse matrix @var{A}. Parameter @var{michol} = \
['off'|'row'|'col'] set if no row nor column sums are preserved, row sums \
are preserved or column sums are preserved respectively.\n\
\n\
For a full description of ILU0 and its options see ilu documentation.\n\
\n\
For more information about the algorithms themselves see:\n\
\n\
[1] Saad, Yousef: Iterative Methods for Sparse Linear Systems. Second Edition. \
Minneapolis, Minnesota: Siam 2003.\n\
\n\
    @seealso{ilu, ilutp, iluc, ichol}\n\
    @end deftypefn")
{
  octave_value_list retval;

  int nargin = args.length ();
  std::string michol = "off";
  double droptol = 0;
 

  if (nargout > 1 || nargin < 1 || nargin > 3)
    {
      print_usage ();
      return retval;
    }

  if (args (0).is_scalar_type () || !args (0).is_sparse_type ())
    error ("icholt: 1. parameter must be a sparse square matrix.");

  if (args (0).is_empty ())
    {
      retval (0) = octave_value (SparseMatrix());
      return retval;
    }

  if (! error_state && (nargin >= 2))
    {
      droptol = args (1).double_value ();
      if (error_state || (droptol < 0) || ! args (1).is_real_scalar ())
        error ("icholt: 2. parameter must be a positive real scalar.");
    }

  if (! error_state && (nargin == 3))
    {
      michol = args (2).string_value ();
      if (error_state || !(michol == "on" || michol == "off"))
        error ("icholt: 3. parameter must be 'on' or 'off' character string.");
    }

  if (!error_state)
    {
      // In ICHOL0 algorithm the zero-pattern of the input matrix is preserved so
      // it's structure does not change during the algorithm. The same input
      // matrix is used to build the output matrix due to that fact.
      octave_value_list param_list;
      if (!args (0).is_complex_type ())
        {
          Array<double> cols_norm;
          SparseMatrix L;
          param_list.append (args (0).sparse_matrix_value ());
          SparseMatrix sm_l = feval ("tril", param_list)(0).sparse_matrix_value (); 
          param_list (1) = 1;
          param_list (2) = "cols";
          cols_norm = feval ("norm", param_list)(0).vector_value ();
          param_list.clear ();
          ichol_t <SparseMatrix, 
                   double, ichol_mult_real, ichol_checkpivot_real> 
                   (sm_l, L, cols_norm.fortran_vec (), droptol, michol);
          if (!error_state)
            retval (0) = octave_value (L);
        }
      else
        {
          Array<Complex> cols_norm;
          SparseComplexMatrix L;
          param_list.append (args (0).sparse_complex_matrix_value ());
          SparseComplexMatrix sm_l = feval ("tril", param_list)(0).sparse_complex_matrix_value (); 
          param_list (1) = "cols";
          cols_norm = feval ("norm", param_list)(0).complex_vector_value ();
          param_list.clear ();
          ichol_t <SparseComplexMatrix, 
                   Complex, ichol_mult_complex, ichol_checkpivot_complex> 
                   (sm_l, L, cols_norm.fortran_vec (), Complex (droptol), michol);
          if (!error_state)
            retval (0) = octave_value (L);
        }

    }

  return retval;
}
/*
%!shared A_1, A_1_in, A_2, A_2_in, A_3, A_3_in, A_4, A_4_in, A_5, A_5_in
%! A_1 = [ 0.37, -0.05,  -0.05,  -0.07;
%!        -0.05,  0.116,  0.0,   -0.05;
%!        -0.05,  0.0,    0.116, -0.05;
%!        -0.07, -0.05,  -0.05,   0.202];
%! A_1_in = sparse(tril (A_1));
%!
%! A_2 = gallery ('poisson', 30);
%! A_2_in = sparse(tril (A_2));
%!
%! A_3 = gallery ('tridiag', 50);
%! A_3_in = sparse(tril (A_3));
%!
%! nx = 400; ny = 200;
%! hx = 1 / (nx + 1); hy = 1 / (ny + 1);
%! Dxx = spdiags ([ones(nx, 1), -2 * ones(nx, 1), ones(nx, 1)], [-1 0 1 ], nx, nx) / (hx ^ 2);
%! Dyy = spdiags ([ones(ny, 1), -2 * ones(ny, 1), ones(ny, 1)], [-1 0 1 ], ny, ny) / (hy ^ 2);
%! A_4 = -kron (Dxx, speye (ny)) - kron (speye (nx), Dyy);
%! A_4_in = sparse(tril (A_4));
%!
 A_5 = [ 0.37, -0.05,          -0.05,  -0.07;
        -0.05,  0.116,          0.0,   -0.05 + 0.05i;
        -0.05,  0.0,            0.116, -0.05;
        -0.07, -0.05 - 0.05i,  -0.05,   0.202];
 A_5 = sparse(A_5);
 A_5_in = sparse(tril (A_5));
%!
%!test
%!error icholt ([]);
%!error icholt ([],[]);
%!error icholt ([],[],[]);
%!error [~] = icholt ([],[],[]);
%!error [L] = icholt ([],[],[]);
%!error [L] = icholt ([], 1e-4, 1)
%!error [L] = icholt (A_1_in, [], 'off')
%!error [L] = icholt (A_1_in, 1e-4, [])
%!error [L, E] = icholt (A_1_in, 1e-4, 'off')
%!error [L] = icholt (A_1_in, 1e-4, 'off', A_1)
%!error icholt (sparse (0), 1e-4, 'off');
%!error icholt (sparse (-0), 1e-4, 'off');
%!error icholt (sparse (-1), 1e-4, 'off');
%!error icholt (sparse (i), 1e-4, 'off');
%!error icholt (sparse (-i), 1e-4, 'off');
%!error icholt (sparse (1 + 1i), 1e-4, 'off');
%!error icholt (sparse (1 - 1i), 1e-4, 'off');
%!
%!test
%! [L] = icholt (sparse (1), 1e-4, 'off');
%! assert (L, sparse (1));
%! [L] = icholt (sparse (4), 1e-4, 'off');
%! assert (L, sparse (2));
%!
%!test
%! [L] = icholt (A_1_in, 1e-4, 'off');
%! assert (norm (A_1 - L*L', 'fro') / norm (A_1, 'fro'), eps, eps)
%! [L] = icholt (A_1_in, 1e-4, 'on');
%! assert (norm (A_1 - L*L', 'fro') / norm (A_1, 'fro'), eps, eps)
%!
%!test
%! [L] = icholt (A_2_in, 1e-4, 'off');
%! assert (norm (A_2 - L*L', 'fro') / norm (A_2, 'fro'), 1e-4, 1e-4)
%! [L] = icholt (A_2_in, 1e-4, 'on');
%! assert (norm (A_2 - L*L', 'fro') / norm (A_2, 'fro'), 3e-4, 1e-4)
%!
%!test
%! [L] = icholt (A_3_in, 1e-4, 'off');
%! assert (norm (A_3 - L*L', 'fro') / norm (A_3, 'fro'), eps, eps)
%! [L] = icholt (A_3_in, 1e-4, 'on');
%! assert (norm (A_3 - L*L', 'fro') / norm (A_3, 'fro'), eps, eps)
%!
%!test
%! [L] = icholt (A_4_in, 1e-4, 'off');
%! assert (norm (A_4 - L*L', 'fro') / norm (A_4, 'fro'), 2e-4, 1e-4)
%! [L] = icholt (A_4_in, 1e-4, 'on');
%! assert (norm (A_4 - L*L', 'fro') / norm (A_4, 'fro'), 7e-4, 1e-4)
%!
%%!test
%%!error [L] = icholt (A_5_in, 1e-4, 'off');
%%!error [L] = icholt (A_5_in, 1e-4, 'on');
*/


