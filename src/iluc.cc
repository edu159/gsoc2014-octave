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
struct list_elem {
    octave_idx_type ridx;
    T data;
};

template <typename octave_matrix_t, typename T>
void ilu_crout (octave_matrix_t& sm_l, octave_matrix_t& sm_u, octave_matrix_t& L, octave_matrix_t& U, const T droptol, const  std::string milu)
  {
  // sm is modified outside so big matrices are not copied twice in memmory 
  // sm is transposed before and after the algorithm because the algorithm is 
  // thought to work with CRS instead of CCS. That does the trick.
  
  // Map the strings into chars to faster comparation inside loops
  // That increase quite a lot the performance!
  #define ROW  1
  #define COL  2
  #define OFF  0
  char opt;
  if (milu == "row")
    opt = ROW;
  else if (milu == "col")
    opt = COL;
  else
    opt = OFF;
  
  const octave_idx_type n = sm_u.cols ();
  sm_u = sm_u.transpose();
  // Extract pointers to the arrays for faster access inside loops
  OCTAVE_LOCAL_BUFFER(T, cols_norm_u, n);
  OCTAVE_LOCAL_BUFFER(T, cols_norm_l, n);
  octave_idx_type* cidx_in_u = sm_u.cidx ();
  octave_idx_type* ridx_in_u = sm_u.ridx ();
  T* data_in_u = sm_u.data();
  octave_idx_type* cidx_in_l = sm_l.cidx ();
  octave_idx_type* ridx_in_l = sm_l.ridx ();
  T* data_in_l = sm_l.data();
  octave_idx_type j1, j2, jrow, i, j, k, jj, p, diag, c, total_len_l, total_len_u,p_perm, res, max_ind;
  T tl, r, aux, maximun;

  T zero = T(0);
  if (droptol > zero)
    {
      cols_2nd_norm <T> (sm_u.cidx (), sm_u.data (), cols_norm_u , n);
      cols_2nd_norm <T> (sm_l.cidx (), sm_l.data (), cols_norm_l , n);
    }

  // Working arrays and permutation arrays
  OCTAVE_LOCAL_BUFFER(typename std::list <list_elem <T> >, cidx_l, n);
  OCTAVE_LOCAL_BUFFER(typename std::list <list_elem <T> >, cidx_u, n);
  OCTAVE_LOCAL_BUFFER(typename std::list <octave_idx_type>, cols_list, n);
  OCTAVE_LOCAL_BUFFER(typename std::list <octave_idx_type>, rows_list, n);
  octave_idx_type w_ind_u, w_len_u, w_len_l, w_ind_l, ndrop_u, ndrop_l;
  T total_sum, partial_col_sum, partial_row_sum;
  std::set<octave_idx_type> iw_l;
  std::set<octave_idx_type> iw_u;
  std::set<octave_idx_type>::iterator  it, it2;
  typename std::list <list_elem <T> >::iterator it_list, it_row_u, it_col_l;
  std::list <octave_idx_type>::iterator it_col_u, it_row_l;
  OCTAVE_LOCAL_BUFFER (T, w_data_l, n);
  OCTAVE_LOCAL_BUFFER (T, w_data_u, n);
  list_elem<T> elm;
  std::list <list_elem <T> > list_u;
  std::list <list_elem <T> > list_l;
  OCTAVE_LOCAL_BUFFER(typename std::list <list_elem <T> >::iterator, Ufirst, n);
  OCTAVE_LOCAL_BUFFER(typename std::list <list_elem <T> >::iterator, Lfirst, n);

  for (i = 0; i < n; i++)
    {
      w_data_u[i] = zero;
      w_data_l[i] = zero;
    }
  total_len_u = 0;
  total_len_l = 0;

  for (k = 0; k < n; k++)
    {
      ndrop_u = 0;
      ndrop_l = 0;
      w_len_l = 0;
      w_len_u = 0;
      
      for (i = cidx_in_l[k]; i < cidx_in_l[k+1]; i++)
        {
          if (ridx_in_l[i] != k)
            {
              iw_l.insert (ridx_in_l[i]);
              w_data_l[ridx_in_l[i]] = data_in_l[i];
              w_len_l++;
            }
        }
      for (i = cidx_in_u[k]; i < cidx_in_u[k+1]; i++)
        {
          iw_u.insert (ridx_in_u[i]);
          w_data_u[ridx_in_u[i]] = data_in_u[i];
          w_len_u++;
        }
      it_row_l = rows_list[k].begin ();
      /**
      printf("k: %d\n", k);
        printf("w_data_l antes: \n");
        for (i = 0; i<n; i++)
          printf("%f ", w_data_l[i]);
        printf("\n");
        printf("w_data_u antes: \n");
        for (i = 0; i<n; i++)
          printf("%f ", w_data_u[i]);
        printf("\n");
      printf("it_row: %d\n", *it_row_l);
      **/
      if (!rows_list[k].empty ())
        while (it_row_l != rows_list[k].end ())
          {
            for (it_row_u = Ufirst[*it_row_l]; it_row_u != cidx_u[*it_row_l].end (); it_row_u++)
              {
                aux = w_data_u[it_row_u->ridx];
                //printf("row_l: %f, row_u: %f , col: %d\n", Lfirst[*it_row_l]->data, it_row_u->data, it_row_u->ridx);
                w_data_u[it_row_u->ridx] -= Lfirst[*it_row_l]->data * it_row_u->data;
                if ((aux == zero) && (w_data_u[it_row_u->ridx] != zero))
                  {
                    iw_u.insert (it_row_u->ridx);
                    w_len_u++;
                  }
                else if ((aux != zero) && (w_data_u[it_row_u->ridx] == zero))
                  ndrop_u++;
              }
            it_row_l++;
          }
            
      it_col_u = cols_list[k].begin ();
      if (!cols_list[k].empty ())
        while (it_col_u != cols_list[k].end ())
          {
            for (it_col_l = Lfirst[*it_col_u]; it_col_l != cidx_l[*it_col_u].end (); it_col_l++)
              {
                if (it_col_l->ridx > k)
                  {
                    aux = w_data_l[it_col_l->ridx];
                    //printf("col_u: %f, col_l: %f\n", Ufirst[*it_col_u]->data, it_col_l->data);
                    w_data_l[it_col_l->ridx] -= Ufirst[*it_col_u]->data * it_col_l->data;
                    if ((aux == zero) && (w_data_l[it_col_l->ridx] != zero))
                      {
                        iw_l.insert (it_col_l->ridx);
                        w_len_l++;
                      }
                    else if ((aux != zero) && (w_data_l[it_col_l->ridx] == zero))
                      ndrop_l++;
                  }
              }
            it_col_u++;
          }
      /**
        printf("w_data_l: \n");
        for (i = 0; i<n; i++)
          printf("%f ", w_data_l[i]);
        printf("\n");
        printf("w_data_u: \n");
        for (i = 0; i<n; i++)
          printf("%f ", w_data_u[i]);
        printf("\n");
        **/
        if (w_data_u[k] == zero)
          {
                error ("ilutp: There is a pivot equal to zero.");
                break;
          }
   //   printf("diagonal: %f\n", w_data_u[k]);
        // Extract the U part
        for (it = iw_u.begin (); it != iw_u.end (); ++it)
        {
          if (w_data_u[*it] != zero)
            {
               
              elm.data = w_data_u[*it];
              elm.ridx = *it;
              cidx_u[k].push_back (elm);
            }
        }
        // Extract the L part
        c = w_len_u - ndrop_u;
        total_len_u += c;
        for (it = iw_l.begin (); it != iw_l.end (); ++it)
        {
          if (w_data_l[*it] != zero)
            {
              elm.data = w_data_l[*it]/w_data_u[k];
              elm.ridx = *it;
              cidx_l[k].push_back (elm);
            }
        }
        /**
        printf("cidx_u: ");
        it_list = cidx_u[k].begin();
        while (it_list != cidx_u[k].end())
        {
          printf ("%d ", (*it_list).ridx);
          it_list++;
        }
        printf("\n");
        printf("cidx_l: ");
        it_list = cidx_l[k].begin();
        while (it_list != cidx_l[k].end())
        {
          printf ("%d ", (*it_list).ridx);
          it_list++;
        }
        printf("\n");
        **/
        c = w_len_l - ndrop_l;
        total_len_l += c;
        // Clear the auxiliar data structures
        for (i=0; i < n; i++)
          {
            w_data_u[i] = zero;
            w_data_l[i] = zero;
          }
        iw_l.clear ();
        iw_u.clear ();
        /**
        printf("rows_list: \n");
        for (i=0;i<n;i++)
        {
          printf("row %d: ", i);
          if (!rows_list[i].empty())
          {
            it_row_l = rows_list[i].begin();
            while (it_row_l != rows_list[i].end())
              {
              printf( "%d ", *it_row_l);
              it_row_l++;
              }
          }
        }
        printf("\n");
        printf("cols_list: \n");
        for (i=0;i<n;i++)
        {
          printf("col %d: ", i);
          if (!cols_list[i].empty())
          {
            it_row_l = cols_list[i].begin();
            while (it_row_l != cols_list[i].end())
              {
              printf( "%d ", *it_row_l);
              it_row_l++;
              }
          }
        }
        printf("\n");
        **/
      if (k < (n-1))
        {
          Ufirst[k] = cidx_u[k].begin ();
          Lfirst[k] = cidx_l[k].begin ();
          for (i = 0; i <= k; i++)
            {
              if (Ufirst[i]->ridx <= k)
                {
                  ++(Ufirst[i]);
                  if (Ufirst[i] != cidx_u[i].end ())
                    {
                      cols_list[Ufirst[i]->ridx].push_back (i);
                    }
                }
              else if (i == k)
                  cols_list[Ufirst[i]->ridx].push_back (i);

              if (Lfirst[i]->ridx <= k)
                {
                   ++(Lfirst[i]);
                   if (Lfirst[i] != cidx_l[i].end ())
                     {
                       rows_list[Lfirst[i]->ridx].push_back (i);
                     }
                }
              else if (i == k)
                  rows_list[Lfirst[i]->ridx].push_back (i);
            }
        }
       // printf("%d ", *(rows_list[Lfirst[k]->ridx].begin()));
       // printf("\n");
       /**
        printf("Lfirst: ");
        for (i=0;i<=k;i++)
          printf("%d ",Lfirst[k]->ridx);
        printf("\n");
        printf("Ufirst: ");
        for (i=0;i<=k;i++)
          printf("%d ",Ufirst[k]->ridx);
        printf("\n");
        **/

    }
  //Expand cidx to fit with the constructor
  if (!error_state) {
    // Expand LIL format into CCS. U matrix.
    // First do with U and then with L. This way less memory is used at the same time.
    // TODO: Maybe inline function?
    Array <octave_idx_type> cidx_out (dim_vector (total_len_u, 1));
    Array <octave_idx_type> ridx_out (dim_vector (total_len_u, 1));
    Array <T> data_out (dim_vector (total_len_u, 1));
    octave_idx_type* cidx_out_pt = cidx_out.fortran_vec ();
    octave_idx_type* ridx_out_pt = ridx_out.fortran_vec ();
    T* data_out_pt = data_out.fortran_vec ();
    p = 0;
    for (i=0 ; i < n; i++) 
    {
      if (!cidx_u[i].empty ())
        {
          it_list = cidx_u[i].begin ();
          while (it_list != cidx_u[i].end ())
            {
              ridx_out_pt[p] = (*it_list).ridx; ;
              data_out_pt[p] = (*it_list).data;
              cidx_out_pt[p] = i;
              it_list++;
              p++;
            }
        }
    }
    /**
        printf("cidx: \n");
        for (i = 0; i<total_len_u; i++)
          printf("%d ", cidx_out_pt[i]);
        printf("\n");
        printf("ridx: \n");
        for (i = 0; i<total_len_u; i++)
          printf("%d ", ridx_out_pt[i]);
        printf("\n");
        printf("data: \n");
        for (i = 0; i<total_len_u; i++)
          printf("%f ", data_out_pt[i]);
        printf("\n");
        **/
    U = octave_matrix_t (data_out, ridx_out, cidx_out, n, n);
    U = U.transpose();
    // Expand LIL format into CCS. L matrix.
    cidx_out = Array<octave_idx_type> (dim_vector (total_len_l, 1));
    ridx_out = Array<octave_idx_type> (dim_vector (total_len_l, 1));
    data_out = Array<T> (dim_vector (total_len_l, 1));
    cidx_out_pt = cidx_out.fortran_vec ();
    ridx_out_pt = ridx_out.fortran_vec ();
    data_out_pt = data_out.fortran_vec ();
    p = 0;
    for (i=0 ; i < n; i++) 
    {
      it_list = cidx_l[i].begin ();
      while (it_list != cidx_l[i].end ())
        {
          ridx_out_pt[p] = (*it_list).ridx; ;
          data_out_pt[p] = (*it_list).data;
          cidx_out_pt[p] = i;
          it_list++;
          p++;
        }
    }
    L = octave_matrix_t (data_out, ridx_out, cidx_out, n, n);
  }
}

DEFUN_DLD (iluc, args, nargout, "-*- texinfo -*-")
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
    droptol = args(2).double_value ();
    milu = args (3).string_value ();

  if (! error_state)
    {
      if (!args (0).is_complex_type ())
      {
        SparseMatrix sm_l = args (0).sparse_matrix_value ();
        SparseMatrix sm_u = args (1).sparse_matrix_value ();
        SparseMatrix U;
        SparseMatrix L;
        ilu_crout <SparseMatrix, double> (sm_l, sm_u, L, U, droptol, milu);
        if (! error_state)
          {
            retval(0) = octave_value (L);
            retval(1) = octave_value (U);
          }
      }
      else
      {
        SparseComplexMatrix sm_l = args (0).sparse_complex_matrix_value ();
        SparseComplexMatrix sm_u = args (1).sparse_complex_matrix_value ();
        SparseComplexMatrix U;
        SparseComplexMatrix L;
        ilu_crout < SparseComplexMatrix, std::complex<double> > (sm_l, sm_u, L, U, std::complex<double> (droptol), milu);
        if (! error_state)
          {
            retval(0) = octave_value (L);
            retval(1) = octave_value (U);
          }
      }

    }

  return retval;
}
