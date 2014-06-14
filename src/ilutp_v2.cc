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
#include <math.h>

template <typename T>
void cols_2nd_norm (const octave_idx_type* cidx, const T* data, T* cols_norm, const octave_idx_type n) {
  octave_idx_type i, j;
  T nrm = T (0);
  for (i = 0; i < n ; i++)
    {
      j = cidx[i];
      while (j < cidx[i+1]) 
        {
          nrm += data[j] * data[j];
          j++;
        }
      cols_norm[i] = sqrt (nrm);
      nrm = T (0);
    }
}


//That can be enhanced using a binary search tree, but for a first approach its fine.
inline void bin_insert(octave_idx_type* iw, octave_idx_type iw_len, octave_idx_type e) 
  {
    octave_idx_type imid, imin, imax, aux;
    imin = 0; imax = iw_len;
    while (true) {
      imid = ceil(imax-imin/2);
      if (e < iw[mid]) 
        {
          if (iw[mid-1] < e) 
            {
              pos = mid;
              break;
            }
          else
              imax = imid - 1; 
        }
      else if (iw[mid] < e)
        {
          if (e < iw[mid+1])
            {
              pos = mid + 1;
              break;
            }
          else
            imin = imid + 1;
        }
    }
    aux = iw[pos];
    for (i = pos; i < iw_len; i++)


  }

template <typename octave_matrix_t, typename T>
void ilu_0 (octave_matrix_t& sm,const  std::string milu) {
  // sm is modified outside so big matrices are not copied twice in memmory 
  // sm is transposed before and after the algorithm because the algorithm is 
  // thought to work with CRS instead of CCS. That does the trick.
  
  const octave_idx_type n = sm.cols ();
  // Extract pointers to the arrays for faster access inside loops
  OCTAVE_LOCAL_BUFFER(T, cols_norm, n);
  cols_2nd_norm <T> (sm.cidx (), sm.data (), cols_norm , n);
  sm = sm.transpose ();
  //TODO Mirar a ver si se puede quitar los punteros a cidx y ridx
  octave_idx_type* cidx_in = sm.cidx ();
  octave_idx_type* ridx_in = sm.ridx ();
  T* data_in = sm.data();
  OCTAVE_LOCAL_BUFFER (octave_idx_type, uptr, n);
  octave_idx_type j1, j2, jrow, i, k, jj, p, diag, c;
  T tl, r;
  char opt;
  T droptol = T(0);
  // Arrays to build later the output matrix and fortran vectors for fast access
  Array <octave_idx_type> cidx_out(dim_vector(n+1,1));
  octave_idx_type* cidx = cidx_out.fortran_vec();
//  Array <octave_idx_type> ridx_out(dim_vector(cidx_in[1] - cidx_in[0],1));
 // Array<T> data_out(dim_vector(cidx_in[1] - cidx_in[0],1));
  Array <octave_idx_type> ridx_out;
  octave_idx_type* ridx ;
  Array<T> data_out;
  T* data;;
  //Ininialise the first row of the output matrix
  for (octave_idx_type i = 0; i <= n; i++)
  {
    cidx[i] = cidx_in[i];
   }
  /**
  for (octave_idx_type i = 0; i< (cidx_in[1] - cidx_in[0]); i++) 
    {
      ridx[i] = ridx_in[i];
      data[i] = data_in[i];
    }
    **/

  // Map the strings into chars to faster comparation inside loops
  // That increase quite a lot the performance!
  OCTAVE_LOCAL_BUFFER (T, w_data, n);
  OCTAVE_LOCAL_BUFFER (octave_idx_type, iw, n);
  OCTAVE_LOCAL_BUFFER (octave_idx_type, jw, n);
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
    }
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
      /**
        printf("k = %d\n",k);
        printf("w_data_1: ");
        for (int ii = 0; ii < n; ii++)
          printf("%f ", w_data[ii]);
        printf("\n");
      printf("\n");
      **/
      r = 0;
      jrow = iw[0];
      printf("K =%d\n", k);
      printf("w_len: %d\n", w_len);
      while (jrow < k) 
        {
          if (w_data[jrow] != T(0))
          {
            if (abs (w_data[jrow]) < (droptol /** *cols_norm[k]**/))
              {
                w_data[jrow] = T(0);
              }
            else
              {
                tl = w_data[jrow] / data[uptr[jrow]];
                w_data[jrow] = tl;
               // printf("cidx = %d ", cidx[jrow+1] -1);
                T aux;
                printf(" jrow: %d \n", jrow);
                for (jj = uptr[jrow] + 1; jj <= cidx[jrow+1]-1; jj++)
                  {
                    aux = w_data[ridx[jj]];
                    //printf("jj = %d ", jj);
                    w_data[ridx[jj]] = w_data[ridx[jj]] - tl * data[jj];
                    if ((aux == T(0)) && (w_data[ridx[jj]] != T(0)))
                        {
                          iw_in[w_len] = ridx[jj];
                          w_len++;
                        }
                    //printf("caca");
                    // That is for the milu='row'
                    // if (opt == 1)
                    //  r += tl*w_data[ridx[jj]];
                  printf("w_data: ");
                  for (int ii = 0; ii < n; ii++)
                    printf("%f ", w_data[ii]);
                  printf("\n");
                  printf("iw: ");
                  for (int ii = 0; ii < w_len; ii++)
                    printf("%d ", iw[ii]);
                  printf("\n");
                  }
                  //printf("\n");
              }
                //printf("j= %d jrow = %d \n",j,  jrow);
          }
          jrow++;
        }
      /**
      if (k != jrow)
        {
          error("ilu0: Your input matrix has a zero in the diagonal.");
          break;
        }
        **/
    //  if (k >0){
        c = 0;
        for (octave_idx_type m = 0; m < n; m++)
          if (w_data[m] != T(0))
            c++;
        printf("Before: ridx length: %d , data length: %d , cidx length: %d\n", ridx_out.length(), data_out.length(), cidx_out.length());
        octave_idx_type l = ridx_out.length();
        ridx_out.resize1(l + c); 
        ridx = ridx_out.fortran_vec();
        data_out.resize1(l + c); 
        data = data_out.fortran_vec();
        printf("After: ridx length: %d , data length: %d , cidx length: %d\n", ridx_out.length(), data_out.length(), cidx_out.length());
        p = 0;
        diag = 0;

        for (octave_idx_type m = 0; m < w_len; m++)
          if (w_data[iw[m]] != T(0))
          {
            if (p != 0) {
              cidx[k+1] = cidx[k] - cidx[0] + c;
            }
            data[l + p] = w_data[iw[m]];
            ridx[l + p] = iw[m];
            if (iw[m] == k)
              diag = l + p;
            p++;
      //  }
        /**
        printf("cidx: ");
        for (int ii = 0; ii < cidx_out.length(); ii++)
          printf("%d ", cidx_out(ii));
        printf("\n");
        printf("ridx: ");
        for (int ii = 0; ii < ridx_out.length(); ii++)
          printf("%d ", ridx_out(ii));
        printf("\n");
        printf("data: ");
        for (int ii = 0; ii < data_out.length(); ii++)
          printf("%f ", data_out(ii));
        printf("\n");
        **/
      }
      //TODO Usar uptr para luego filtrar U ya que sabemos donde están las normales
      uptr[k] = diag;
      printf("uptr[k] %d, jrow %d\n", uptr[k],jrow);
      // That is for the milu='row'
      if(opt == 1)
        data[uptr[k]] -= r;
      /**
      if (data[diag] == T(0))
        {
          error ("ilu0: There is a pivot equal to zero.");
          break;
        }
      **/
      //TODO Hacer 0 usando cidx ya que tenemos la información de los que son != 0
      for(i = 0; i < n ; i++)
        {
          w_data[i] = 0;
          iw[i] = -1;
          jw[i] = -1;
        }
    }

  //Expand cidx to fit with the constructor
  if (!error_state) {
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
    printf("holaa mama");
    sm = octave_matrix_t(data_out, ridx_out, cidx_out_expand);
/**
    for (int i = 0; i < n; i++)
      printf("%d:%d ", i,cidx_out_expand(i));
    printf("holaa\n");
    **/
  }
  sm = sm.transpose();
}

DEFUN_DLD (ilutp, args, nargout, "-*- texinfo -*-")
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
      error ("ilu0: Second parameter should be 'row' or 'col'.");
    else
      milu = args (1).string_value ();

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
