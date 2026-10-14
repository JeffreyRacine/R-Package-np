void np_p_ukernelv(const int KERNEL, 
                   const int P_KERNEL,
                   const int P_IDX,
                   const int P_NIDX,
                   const double * const xt, const int num_xt, 
                   const int do_xw,
                   const double x, const double lambda, const int ncat,
                   const double cat,
                   double * const result,
                   double * const p_result,
                   const XL * const xl,
                   const XL * const p_xl,
                   const int swap_xxt,
                   const int do_ocg,
                   const int do_perm,
                   double * const scratch_kbuf){

  /* 
     this should be read as:
     an array of constant pointers to functions that take a double
     and return a double
  */

  int i,j,r,l; 
  const int bin_do_xw = do_xw > 0;
  double unit_weight = 1.0;
  double * const xw = (bin_do_xw ? result : &unit_weight);

  double * const pxw = (bin_do_xw ? p_result : &unit_weight);

  const double ex = do_ocg ? cat : x;

  double (* const k[])(int, double, int) = { np_uaa, np_unli_racine,
                                             np_econvol_uaa, np_econvol_unli_racine,
                                             np_score_uaa, np_score_unli_racine,
                                             np_score_convol_uaa, np_score_convol_unli };
  const int kernel = KERNEL;
  const int p_kernel = P_KERNEL;

  /* Unordered kernels depend only on same/different category; cache both values once. */
  const double kn_same = k[kernel](1, lambda, ncat);
  const double kn_diff = k[kernel](0, lambda, ncat);
  const double pkn_same = k[p_kernel](1, lambda, ncat);
  const double pkn_diff = k[p_kernel](0, lambda, ncat);
  const int use_const_k = NP_CF316_DISC_UPPER(kernel, lambda, ncat) &&
    NP_CF316_DISC_CONSTANT(kn_same, kn_diff);
  const int use_const_pk = NP_CF316_DISC_UPPER(p_kernel, lambda, ncat) &&
    NP_CF316_DISC_CONSTANT(pkn_same, pkn_diff);
  const double kn_const = 0.5*(kn_same + kn_diff);
  const double pkn_const = 0.5*(pkn_same + pkn_diff);
  const int p_iscat_const = (swap_xxt && do_ocg);
  const int p_iscat_const_val = (cat == x);

  double *kbuf = scratch_kbuf;
  const int own_kbuf = (kbuf == NULL);
  if(own_kbuf){
    kbuf = (double *)malloc(num_xt*sizeof(double));
    if(kbuf == NULL) error("memory allocation failed");
  }

  if(xl == NULL){
    for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
      const double kn = use_const_k ?
        kn_const :
        (((xt[i] == x) ? kn_same : kn_diff));

      result[i] = xw[j]*kn;
      kbuf[i] = kn;

      if(do_perm){
        const double pkn = use_const_pk ?
          pkn_const :
          (((p_iscat_const ? p_iscat_const_val : (xt[i] == ex)) ? pkn_same : pkn_diff));
        p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*pkn;
      }
    }

    for(l = 0, r = 0; l < P_NIDX; l++, r += bin_do_xw){
      if((l == P_IDX) && do_perm) continue;
      for (i = 0, j = 0; i < num_xt; i++, j += bin_do_xw){
        p_result[l*num_xt + i] = pxw[r*num_xt + j]*kbuf[i];
      }
    }

  }
  else{
    for (int m = 0; m < xl->n; m++){
      const int istart = xl->istart[m];
      const int nlev = xl->nlev[m];
      for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
        const double kn = use_const_k ?
          kn_const :
          (((xt[i] == x) ? kn_same : kn_diff));

        result[i] = xw[j]*kn;
        kbuf[i] = kn;
      }
    }

    if(do_perm){
      for (int m = 0; m < p_xl->n; m++){
        const int istart = p_xl->istart[m];
        const int nlev = p_xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          const double pkn = use_const_pk ?
            pkn_const :
            (((p_iscat_const ? p_iscat_const_val : (xt[i] == ex)) ? pkn_same : pkn_diff));
          p_result[P_IDX*num_xt + i] = pxw[bin_do_xw*P_IDX*num_xt + j]*pkn;
        }
      }
    }


    for(l = 0, r = 0; l < P_NIDX; l++, r+=bin_do_xw){
      if((l == P_IDX) && do_perm) continue;
      for (int m = 0; m < xl->n; m++){
        const int istart = xl->istart[m];
        const int nlev = xl->nlev[m];
        for (i = istart, j = bin_do_xw*istart; i < istart+nlev; i++, j += bin_do_xw){
          p_result[l*num_xt + i] = pxw[r*num_xt + j]*kbuf[i];
        }
      }
    }
  }

  if(own_kbuf)
    free(kbuf);
}
