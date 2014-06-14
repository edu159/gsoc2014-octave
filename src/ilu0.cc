/**
 * Copyright (C) 2013 Eduardo Ramos Fernández <eduradical951@gmail.com>
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


template <typename octave_matrix_t, typename T>
void ilu_0(octave_matrix_t& sm,const  std::string milu = "") {
  // sm is modified outside so big matrices are not copied twice in memmory 
  // sm is transposed before and after the algorithm because the algorithm is 
  // thought to work with CRS instead of CCS. That does the trick.
  sm = sm.transpose ();
  const octave_idx_type n = sm.cols();
  // Extract pointers to the arrays for faster access inside loops
  octave_idx_type* cidx = sm.cidx();
  octave_idx_type* ridx = sm.ridx();
  T* data = sm.data();
  OCTAVE_LOCAL_BUFFER (octave_idx_type, iw, n);
  OCTAVE_LOCAL_BUFFER (octave_idx_type, uptr, n);
  octave_idx_type j1, j2, jrow, jw, i, k, jj;
  T tl, r;
  char opt;
  // Map the strings into chars to faster comparation inside loops
  // That increase quite a lot the performance!
  if (milu == "row")
    opt = 1;
  else if (milu == "col")
    opt = 2;
  else
    opt = 0;

   
  for (i = 0; i < n; i++)
    iw[i] = -1;
  for (k = 0; k < n; k++)
    {
      j1 = cidx[k];
      j2 = cidx[k+1] - 1;
      octave_idx_type j;
      for (j = j1; j <= j2; j++)
        {
          iw[ridx[j]] = j;
        }
      for (i = 0; i < n; i++)
      r = 0;
      j = j1;
      jrow = ridx[j];
      while (jrow < k && j <= j2) 
        {
          tl = data[j] / data[uptr[jrow]];
          data[j] = tl;
          for (jj = uptr[jrow] + 1; jj <= cidx[jrow+1]-1; jj++)
            {
              jw = iw[ridx[jj]];
              if (jw != -1)
                data[jw] = data[jw] - tl * data[jj];
              else
                // That is for the milu='row'
                if (opt == 1)
                  r += tl*data[jj];
            }
          j++;
          jrow = ridx[j];
        }
      if (k != jrow)
        {
          error("ilu0: Your input matrix has a zero in the diagonal.");
          break;
        }
      uptr[k] = j;
      // That is for the milu='row'
      if(opt == 1)
        data[uptr[k]] -= r;
      if (data[j] == T(0))
        {
          error ("ilu0: There is a pivot equal to zero.");
          break;
        }
      for(i = j1; i <= j2 ; i++)
        iw[ridx[i]] = -1;
    }
    sm = sm.transpose ();
}

DEFUN_DLD (ilu0, args, nargout, "-*- texinfo -*-")
{
  octave_value_list retval;

  int nargin = args.length ();
  std::string milu = "";
 

  if (nargout > 3 || nargin < 1 || nargin > 2)
    {
      print_usage ();
      return retval;
    }

  if (nargin == 2) 
    if (args (1).is_string ())
      milu = args (1).string_value ();
    else
      error ("ilu0: Second parameter should be 'row' or 'col'.");
  if (! error_state)
    {
      
      if (!args (0).is_complex_type ())
      {
        SparseMatrix sm = args (0).sparse_matrix_value ();
        ilu_0<SparseMatrix, double> (sm, milu);
        if (! error_state)
          retval(0) = octave_value (sm);
      }
      else
      {
        SparseComplexMatrix sm = args (0).sparse_complex_matrix_value ();
        ilu_0 < SparseComplexMatrix, std::complex<double> > (sm, milu);
        if (! error_state)
          retval(0) = octave_value (sm);
      }

    }

  return retval;
}
