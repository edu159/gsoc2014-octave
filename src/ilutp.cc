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


template <typename octave_matrix_t, typename T>
void ilu_tp (octave_matrix_t& sm, octave_matrix_t& L, octave_matrix_t& U, octave_idx_type* perm, const T droptol, const T thresh, const  std::string milu, const bool udiag)
  {
  // sm is modified outside so big matrices are not copied twice in memmory 
  // sm is transposed before and after the algorithm because the algorithm is 
  // thought to work with CRS instead of CCS. That does the trick.
  
  const octave_idx_type n = sm.cols ();
  // Extract pointers to the arrays for faster access inside loops
  OCTAVE_LOCAL_BUFFER(T, cols_norm, n);
  octave_idx_type* cidx_in = sm.cidx ();
  octave_idx_type* ridx_in = sm.ridx ();
  T* data_in = sm.data();
//  OCTAVE_LOCAL_BUFFER (octave_idx_type, uptr, n);
  octave_idx_type j1, j2, jrow, i, j, k, jj, p, diag, c, total_len_l, total_len_u,p_perm, res, max_ind;
  T tl, r, aux, maximun;
  char opt;

  // Data for L
  // ridx and data are fixed to the maximun possible ((n^2+n)/2)
  Array <octave_idx_type> cidx_out_l(dim_vector(n+1,1));
  octave_idx_type* cidx_l = cidx_out_l.fortran_vec();
  Array <octave_idx_type> ridx_out_l(dim_vector((n*n+n)/2,1));
  octave_idx_type* ridx_l = ridx_out_l.fortran_vec();
  Array<T> data_out_l(dim_vector((n*n+n)/2,1));
  T* data_l = data_out_l.fortran_vec();
  // Data for U
  Array <octave_idx_type> cidx_out_u(dim_vector(n+1,1));
  octave_idx_type* cidx_u = cidx_out_u.fortran_vec();
  Array <octave_idx_type> ridx_out_u(dim_vector((n*n+n)/2,1));
  octave_idx_type* ridx_u = ridx_out_u.fortran_vec();
  Array<T> data_out_u(dim_vector((n*n+n)/2,1));
  T* data_u = data_out_u.fortran_vec();

  T zero = T(0);
  if (droptol > zero)
    cols_2nd_norm <T> (sm.cidx (), sm.data (), cols_norm , n);
  

  cidx_l[0] = cidx_in[0];
  cidx_u[0] = cidx_in[0];
  for (i = 0; i < ((n*n)+n)/2; i++)
    {
      ridx_u[i] = 0;
      data_u[i] = 0;
      ridx_l[i] = 0;
      data_l[i] = 0;
    }

  // Working arrays and permutation arrays
  octave_idx_type w_ind_u, w_len_u, w_len_l, w_ind_l, ndrop_u, ndrop_l;
  T col_sum, partial_col_sum;
  std::set<octave_idx_type> iw_l;
  std::set<octave_idx_type> iw_u;
  std::set<octave_idx_type>::iterator  it, it2;
  //std::set<octave_idx_type>::reverse_iterator  rit;

  OCTAVE_LOCAL_BUFFER (T, w_data, n);
  OCTAVE_LOCAL_BUFFER (octave_idx_type, iperm, n);

  // Map the strings into chars to faster comparation inside loops
  // That increase quite a lot the performance!
  #define ROW  1
  #define COL  2
  #define OFF  0
  if (milu == "row")
    opt = ROW;
  else if (milu == "col")
    opt = COL;
  else
    opt = OFF;
  for (i = 0; i < n; i++)
    {
      w_data[i] = 0;
      perm[i] = i;
      iperm[i] = i;
    }
  total_len_u = 0;
  total_len_l = 0;

  for (k = 0; k < n; k++)
    {
      w_len_u = 0;
      w_len_l = 0;
      // Initialize the working data vector with a new row from the input matrix
      // Also the working binary trees for indexes of the U and L part are initialized
      // This block of code is ready to handle the milu="row" implementation.
      for (j = cidx_in[k]; j < cidx_in[k+1]; j++)
        {
          p_perm = iperm[ridx_in[j]];
          w_data[iperm[ridx_in[j]]] = data_in[j];
          if (p_perm > k)
            if (opt == ROW) 
              {
                iw_u.insert(ridx_in[j]);
                w_len_u++;
              }
            else
              {
                iw_l.insert(ridx_in[j]);
                w_len_l++;
              }
          else
            if (opt == ROW) 
              {
                iw_l.insert(p_perm);
                w_len_l++;
              }
            else
              {
                iw_u.insert(p_perm);
                w_len_u++;
              }
        }

      
      r = 0;
      ndrop_u = 0;
      ndrop_l = 0;
      //TODO: Implementation for milu="row". It is nice to place it here or maybe separated?
      if (opt == ROW)
        {
          it = iw_l.begin(); 
          while ((jrow < k) && (it != iw_l.end())) 
            {
              if (w_data[jrow] != zero)
                {
                  if (std::abs (w_data[jrow]) < (droptol * cols_norm[jrow]))
                    {
                      w_data[jrow] = zero;
                      ndrop_l++;
                    }
                  else
                    {
                      //TODO: Finish milu="row" algorithm
                      /**
                      tl = w_data[jrow] / data[jrow];
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
                        **/
                    }
                }
              *(++it);
            }
        }
      else
        {
          it = iw_u.begin ();
          jrow = *it;
          col_sum = zero;
          while ((jrow < k) && (it != iw_u.end ())) 
            {
              if (opt == COL)
                partial_col_sum = w_data[jrow];
              if (w_data[jrow] != zero)
                {
                  for (jj = cidx_l[jrow]; jj < cidx_l[jrow+1]; jj++)
                    {
                      p_perm = iperm[ridx_l[jj]];
                      aux = w_data[p_perm];
                      tl = data_l[jj] * w_data[jrow];
                      if (opt == COL)
                        partial_col_sum += tl;
                      w_data[p_perm] -= tl;

                      // In this case a new element has appeared in a position with 0
                      // TODO: Maybe nodes from the binary tree that correspond to indexes of elements
                      // that are no longer non-zero should be removed. The strategy adopted is to
                      // dealy the removal when expanding the w_data to the U and L sparse structure.
                      if ((aux == zero) && (w_data[p_perm] != zero))
                        {
                          if (p_perm > k)
                            {
                              iw_l.insert (ridx_l[jj]);
                              w_len_l++;
                            }
                          else
                            {
                              iw_u.insert (p_perm);
                              w_len_u++;
                            }
                        }
                      // The case when a zero is obtained due to a gaussian elimination step
                      else if ((aux != zero) && (w_data[p_perm] == zero))
                        {
                          if (p_perm > k)
                            ndrop_l++;
                          else
                            ndrop_u++;
                        }
                    }
                  // Drop elements from the U part
                  if (std::abs (w_data[jrow]) < (droptol * cols_norm[k]))
                    {
                      if (opt == COL)
                        col_sum += partial_col_sum;
                      w_data[jrow] = zero;
                      ndrop_u++;
                    }

                }
              jrow = *(++it);
            }


          // Search for the pivot and update iw_l and iw_u if the pivote is not the diagonal element
          if ((thresh) > zero && (k <(n-1)))
            {
              maximun =std::abs (w_data[k]) / thresh;
              max_ind = perm[k];
              for (it = iw_l.begin (); it != iw_l.end (); ++it) 
                {
                  p_perm = iperm[*it];
                  if (std::abs (w_data[p_perm]) > maximun)
                    {
                      maximun = w_data[p_perm];
                      max_ind = *it;
                      it2 = it; 
                    }
                }
              // If the pivot is not the diagonal element update all.
              p_perm = iperm[max_ind];
              if (max_ind != perm[k])
                {
                  iw_l.erase (it2);
                  if (w_data[k] != zero)
                    iw_l.insert (perm[k]);
                  else
                    {
                      ndrop_u--;
                      ndrop_l++;
                      iw_u.insert (k);
                     }
                  // Swap data and update permutation vectors
                  aux = w_data[k];
                  iperm[perm[p_perm]] = k;
                  iperm[perm[k]] = p_perm;
                  c = perm[k];
                  perm[k] = perm[p_perm];
                  perm[p_perm] = c;
                  w_data[k] = w_data[p_perm];
                  w_data[p_perm] = aux;
                }
              
          }              
          // Add the L terms that  are going to be dropped to the diagonal
          if (opt == COL)
            for (it = iw_l.begin (); it != iw_l.end (); ++it) 
              {
                  p_perm = iperm[*it];
                  if (droptol > zero)
                    if (std::abs (w_data[p_perm]) < (droptol * cols_norm[k]))
                      if (opt == COL)
                        col_sum += w_data[p_perm];
              }

          if (opt == COL)
            w_data[k] += col_sum;

          // Scale al the L terms by the pivote and drop the elements  
          for (it = iw_l.begin (); it != iw_l.end (); ++it) 
            {
                p_perm = iperm[*it];
                if (droptol > zero)
                  if (std::abs (w_data[p_perm]) < (droptol * cols_norm[k]))
                    {
                      w_data[p_perm] = zero;
                      ndrop_l++;
                      continue;
                    }
                w_data[p_perm] = w_data[p_perm] / w_data[k];
            }
          

        }

        if (w_data[k] == zero)
          {
            if (udiag == 1)
              {
                w_data[k] = droptol;
                iw_u.insert(k);
                w_len_u++;
                printf("Diagonal sustituida\n");
              }
            else
              {
                error ("ilutp: There is a pivot equal to zero.");
                break;
              }
          }

        c = w_len_u - ndrop_u;
        p = 0;

        // Extract the U part
        cidx_u[k+1] = cidx_u[k] - cidx_u[0] + c;
        for (it = iw_u.begin (); it != iw_u.end (); ++it)
        {
          if (w_data[*it] != zero)
          {
            data_u[total_len_u + p] = w_data[*it];
            ridx_u[total_len_u + p] = *it;
            p++;
          }
        }
        total_len_u += c;

        // Extract the L part
        c = w_len_l - ndrop_l;
        p = 0;
        cidx_l[k+1] = cidx_l[k] - cidx_l[0] + c;
        for (it = iw_l.begin (); it != iw_l.end (); ++it)
        {
          p_perm = iperm[*it];
          if (w_data[p_perm] != zero)
          {
            data_l[total_len_l + p] = w_data[p_perm];
            ridx_l[total_len_l + p] = *it;
            p++;
          }
        }
        total_len_l += c;

        // Clear the auxiliar data structures
        for (i=0; i < n; i++)
          {
            w_data[i] = 0;
          }
        iw_l.clear ();
        iw_u.clear ();
    }
  //Expand cidx to fit with the constructor
  if (!error_state) {
    ridx_out_l.resize1 (total_len_l); 
    data_out_l.resize1 (total_len_l); 
    ridx_out_u.resize1 (total_len_u); 
    data_out_u.resize1 (total_len_u); 
    Array <octave_idx_type> cidx_out_expand_l (dim_vector (ridx_out_l.length (),1));
    Array <octave_idx_type> cidx_out_expand_u (dim_vector (ridx_out_u.length (),1));
    octave_idx_type* cidx_out_pt_l = cidx_out_expand_l.fortran_vec ();
    octave_idx_type* cidx_out_pt_u = cidx_out_expand_u.fortran_vec ();
    octave_idx_type j;
    // Expand L and U
    for (octave_idx_type i=0 ; i<n; i++) 
    {
      j = cidx_l[i];
      while (j < cidx_l[i+1])
      {
        cidx_out_pt_l[j] = i;
        j++;
      }
      j = cidx_u[i];
      while (j < cidx_u[i+1])
      {
        cidx_out_pt_u[j] = i;
        j++;
      }
    }
    U = octave_matrix_t (data_out_u, ridx_out_u, cidx_out_expand_u, n, n);
    L = octave_matrix_t (data_out_l, ridx_out_l, cidx_out_expand_l, n, n);
  }
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
    droptol = args(1).double_value ();
    thresh = args(2).double_value ();
    milu = args (3).string_value ();
    udiag = args(4).bool_value();

  if (! error_state)
    {
      if (!args (0).is_complex_type ())
      {
        SparseMatrix sm = args (0).sparse_matrix_value ();
        Array <octave_idx_type> iperm (dim_vector (sm.cols (), 1)); 
        SparseMatrix U;
        SparseMatrix L;
        ilu_tp<SparseMatrix, double> (sm, L, U, iperm.fortran_vec (), droptol, thresh, milu, udiag);
        if (! error_state)
          {
            retval(0) = octave_value (L);
            retval(1) = octave_value (U);
            retval(2) = iperm.as_row ();
          }
      }
      else
      {
        SparseComplexMatrix sm = args (0).sparse_complex_matrix_value ();
        Array <octave_idx_type> iperm (dim_vector (sm.cols (), 1)); 
        SparseComplexMatrix U;
        SparseComplexMatrix L;
        ilu_tp < SparseComplexMatrix, std::complex<double> > (sm, L, U, iperm.fortran_vec (), std::complex<double> (droptol), 
                                                              std::complex<double> (thresh), milu, udiag);
        if (! error_state)
          {
            retval(0) = octave_value (L);
            retval(1) = octave_value (U);
            retval(2) = iperm.as_row ();
          }
      }

    }

  return retval;
}
