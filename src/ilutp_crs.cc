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


// Function to compute the 2nd-norm of each column of the matrix.
template <typename T>
void cols_2nd_norm (const octave_idx_type* cidx, const T* data, T* cols_norm, const octave_idx_type n) {
  octave_idx_type i, j;
  T nrm = T (0);
  for (i = 0; i < n ; i++)
    {
      j = cidx[i];
      while (j < cidx[i+1]) 
        {
          nrm += std::abs(data[j]) * std::abs(data[j]);
          j++;
        }
      cols_norm[i] = sqrt (nrm);
      nrm = T(0);
    }
}


template <typename T>
octave_idx_type max_row(const octave_idx_type* cidx, const octave_idx_type* ridx, const T* data, const T thresh, const octave_idx_type row)
  {
    octave_idx_type i, max_ind; 
    T maximun, max_candidate;
    maximun = std::abs(data[cidx[row]]);
    max_ind = ridx[cidx[row]];
    //FIXME:Better perform a binary search
    for (i = cidx[row] ; i < (cidx[row+1] -1); i++)
      {
        //Only search in subdiagonal elements
        if (ridx[i] >= row)
          {
            if (ridx[i] == row)
              max_candidate = std::abs(data[i])/thresh;
            else
              max_candidate = std::abs(data[i]);
            if (max_candidate > maximun)
              {
                maximun = max_candidate;
                max_ind = ridx[i];
              }
          }     
      }
    return max_ind;
  }


template <typename octave_matrix_t, typename T>
void ilu_tp (octave_matrix_t& sm, const T droptol, const T thresh, const  std::string milu, const bool udiag) {
  // sm is modified outside so big matrices are not copied twice in memmory 
  // sm is transposed before and after the algorithm because the algorithm is 
  // thought to work with CRS instead of CCS. That does the trick.
  
  const octave_idx_type n = sm.cols ();
  // Extract pointers to the arrays for faster access inside loops
  OCTAVE_LOCAL_BUFFER(T, cols_norm, n);
  cols_2nd_norm <T> (sm.cidx (), sm.data (), cols_norm , n);
  sm = sm.transpose ();
  octave_idx_type* cidx_in = sm.cidx ();
  octave_idx_type* ridx_in = sm.ridx ();
  T* data_in = sm.data();
  OCTAVE_LOCAL_BUFFER (octave_idx_type, uptr, n);
  octave_idx_type j1, j2, jrow, i, j, k, jj, p, diag, c, ndrop, total_len;
  T tl, r, aux;
  char opt;
  // Arrays to build later the output matrix and fortran vectors for fast access
  Array <octave_idx_type> cidx_out(dim_vector(n+1,1));
  octave_idx_type* cidx = cidx_out.fortran_vec();
  Array <octave_idx_type> ridx_out(dim_vector(n*n,1));
  octave_idx_type* ridx = ridx_out.fortran_vec();
  Array<T> data_out(dim_vector(n*n,1));
  T* data = data_out.fortran_vec();
  T zero = T(0);

  for (octave_idx_type i = 0; i <= n; i++)
  {
    cidx[i] = cidx_in[i];
   }
  for (octave_idx_type i = 0; i < (n*n); i++)
  {
    ridx[i] = 0;
    data[i] = 0;
   }
/**
  octave_idx_type res ;
  for (i = 0; i < n; i++)
  {
    res= max_row(cidx_in, ridx_in, data_in, thresh, i);
    printf("row %d, %d\n", i, res);
  }
**/
  // Map the strings into chars to faster comparation inside loops
  // That increase quite a lot the performance!
  OCTAVE_LOCAL_BUFFER (T, w_data, n);
  OCTAVE_LOCAL_BUFFER (octave_idx_type, iw, n);
  OCTAVE_LOCAL_BUFFER (octave_idx_type, jw, n);
  OCTAVE_LOCAL_BUFFER (octave_idx_type, perm, n);
  octave_idx_type w_ind, w_len;
  if (milu == "row")
    opt = 1;
  else if (milu == "col")
    opt = 2;
  else
    opt = 0;
   
  for (i = 0; i < n; i++)
    {
      w_data[i] = 0;
      iw[i] = -1;
      jw[i] = -1;
      perm[i] = i;
    }
  total_len = 0;
  for (k = 0; k < n; k++)
    {
      j1 = cidx_in[k];
      j2 = cidx_in[k+1] - 1;
      octave_idx_type j;
      //Initialize the working data vector with a new row from the input matrix
      w_ind = 0;
      for (j = j1; j <= j2; j++)
        {
          
          jw[ridx_in[j]] = w_ind;
          w_data[ridx_in[j]] = data_in[j];
          iw[w_ind] = ridx_in[j];
          w_ind++;
        }
      w_len = w_ind;
      r = 0;
      jrow = iw[0];
      ndrop = 0;
      //Find the best pivote
   
      //res= max_row(cidx_in, ridx_in, data_in, thresh, k);
      //perm[res]
      
      while (jrow < k) 
        {
          if (w_data[jrow] != zero)
          {
            if (std::abs (w_data[jrow]) < (droptol * cols_norm[jrow]))
              {
                w_data[jrow] = zero;
                ndrop++;
              }
            else
              {
                tl = w_data[jrow] / data[uptr[jrow]];
                w_data[jrow] = tl;
                aux;
                for (jj = uptr[jrow] + 1; jj <= cidx[jrow+1]-1; jj++)
                  {
                    aux = w_data[ridx[jj]];
                    w_data[ridx[jj]] = w_data[ridx[jj]] - tl * data[jj];
                    if ((aux == zero) && (w_data[ridx[jj]] != zero))
                        {
                          iw[w_len] = ridx[jj];
                          jw[ridx[jj]] = w_len;
                          w_len++;
                        }
                    else if ((aux != zero) && (w_data[ridx[jj]] == zero))
                      ndrop++;
                    // That is for the milu='row'
                    // if (opt == 1)
                    //  r += tl*w_data[ridx[jj]];
                  }
              }
          }
          jrow++;
        }
        if (w_data[k] == zero)
          {
            error ("ilu0: There is a pivot equal to zero.");
            break;
          }
        c = w_len - ndrop;
        /**
        octave_idx_type l = ridx_out.length();
        ridx_out.resize1(l + c); 
        ridx = ridx_out.fortran_vec();
        data_out.resize1(l + c); 
        data = data_out.fortran_vec();
        **/
        p = 0;
        diag = 0;

        for (octave_idx_type m = 0; m < n; m++)
        {
          if (w_data[m] != zero)
          {
            if (p == 0) {
              cidx[k+1] = cidx[k] - cidx[0] + c;
            }
            data[total_len + p] = w_data[m];
            ridx[total_len + p] = iw[jw[m]];
            if (iw[jw[m]] == k)
              diag = total_len + p;
            p++;
          }
      }
        total_len += c;
      uptr[k] = diag;
      // That is for the milu='row'
      //if(opt == 1)
      //  data[uptr[k]] -= r;
      p = 0;
      for(i = cidx[k]; i < cidx[k+1]; i++)
        {

          w_data[ridx[i]] = 0;
          iw[p] = -1;
          jw[ridx[i]] = -1;
          p++;
        }
    }
  //Drop the elements in U according to droptol 
  for (i = 0; i < n ; i++)
  {
    for (j = (uptr[i] + 1); j < (cidx[i+1]); j++)
    {
      if (std::abs (data[j]) < (droptol * cols_norm[ridx[j]]))
        data[j] = zero;
    }
   }
  //Expand cidx to fit with the constructor
  if (!error_state) {
    printf("total_len%d", total_len);
        ridx_out.resize1(total_len); 
        data_out.resize1(total_len); 
    Array <octave_idx_type> cidx_out_expand(dim_vector(ridx_out.length(),1));
    octave_idx_type* cidx_out_pt = cidx_out_expand.fortran_vec();
    octave_idx_type j;
    for (octave_idx_type i=0 ; i<n; i++) 
    {
      j = cidx[i];
      while (j < cidx[i+1])
      {
        cidx_out_pt[j] = i;
        j++;
      }
    }
    sm = octave_matrix_t(data_out, ridx_out, cidx_out_expand);

    
  }
  sm = sm.transpose();
}

DEFUN_DLD (ilutp, args, nargout, "-*- texinfo -*-")
{
  octave_value_list retval;

  int nargin = args.length ();
  std::string milu = "";
  double droptol, thresh;
  bool udiag;

 

  if (nargout > 3 || nargin < 1 || nargin > 5)
    {
      print_usage ();
      return retval;
    }
    //That is really nasty, valid for input checked on ilu.m and called ilutp with 5 params.
    //FIXME: Have to implement the apropiate checking.
    droptol = args(1).double_value();
    thresh = args(2).double_value();
    milu = args (3).string_value ();
    udiag = args(4).bool_value();

  if (! error_state)
    {
      
      if (!args (0).is_complex_type ())
      {
        SparseMatrix sm = args (0).sparse_matrix_value ();
        ilu_tp<SparseMatrix, double> (sm, droptol, thresh, milu, udiag);
        if (! error_state)
          retval(0) = octave_value (sm);
      }
      else
      {
        SparseComplexMatrix sm = args (0).sparse_complex_matrix_value ();
        ilu_tp < SparseComplexMatrix, std::complex<double> > (sm, std::complex<double>(droptol), 
                                                              std::complex<double>(thresh), milu, udiag);
        if (! error_state)
          retval(0) = octave_value (sm);
      }

    }

  return retval;
}
