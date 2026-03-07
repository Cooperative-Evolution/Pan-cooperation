// this model is an adaptation of the Pan model for a human group
// there are four items/behavior

functions {
  vector chol2corvec(matrix inputobj, int n_axes, int n_cor) {
    matrix[n_axes, n_axes] cor_mat = multiply_lower_tri_self_transpose(inputobj);
    vector[n_cor] out;
    int idx = 1;
    for (k in 1:rows(cor_mat)) {
      for (j in 1:(k - 1)) {
        out[idx] = cor_mat[j, k];
        idx = idx + 1;
      }
    }
  return out;
  }

}
data {
  int<lower=0> n_ids;
  int<lower=0> n_dyads;
  int<lower=0> n_obs;
  int<lower=0> n_behs;
  int<lower=0> n_cors;

  // indices per observation
  array[n_obs] int<lower=0> id1_for_obs;
  array[n_obs] int<lower=0> id2_for_obs;
  array[n_obs] int<lower=0> dyad_for_obs;

  array[n_dyads] int<lower=1,upper=3> combi_for_dyad;
  array[n_ids] int<lower=1,upper=2> sex_for_id;

  // dyadic covariates
  vector<lower=0,upper=1>[n_dyads] relatedness_deg;
  vector<lower=0,upper=1>[n_dyads] dyad_is_spouse;

  // interaction data
  array[n_obs] int<lower=0,upper=4> food;
  array[n_obs] int<lower=0,upper=4> assoc;
  array[n_obs] int<lower=0,upper=4> help;
  array[n_obs] int<lower=0,upper=4> friendship;
  array[n_obs] int<lower=0,upper=4> obseff;

  // indices
  int<lower=1> n_females;
  array[n_females] int<lower=1> index_female_ids;
  int<lower=1> n_males;
  array[n_males] int<lower=1> index_male_ids;

  int<lower=1> n_ff;
  array[n_ff] int<lower=1> index_ff_dyads;
  int<lower=1> n_mix;
  array[n_mix] int<lower=1> index_mix_dyads;
  int<lower=1> n_mm;
  array[n_mm] int<lower=1> index_mm_dyads;
}

transformed data {
  // centered version of spouse/kin predictor
  vector[n_dyads] relatedness_deg_std = relatedness_deg - mean(relatedness_deg);
  vector[n_dyads] dyad_is_spouse_std = dyad_is_spouse - mean(dyad_is_spouse);
}

parameters {
  matrix<lower=0>[n_behs, 2] sigma_sexspec; // [X behaviours, 2 sex]
  matrix[n_behs, 2] icpts_sexspec;
  matrix[n_behs, n_ids] indivals_z;

  matrix<lower=0>[n_behs, 3] sigma_combispec; // [X behaviours, 3 sexcombi]
  matrix[n_behs, 3] icpts_combispec;
  matrix[n_behs, n_dyads] dyadvals_z;

  vector[n_behs] b; // 1=food,2=assoc,3=help,4=fried
  vector[n_behs] b_relatedness_deg; // estimate kin effect for dyads
  vector[n_behs] b_spouse; // estimate spouse effect for dyads


  cholesky_factor_corr[n_behs] chol_indi_dims_sexspec1; // cholesky factor of correlation matrix for individual-level dims FEMALES
  cholesky_factor_corr[n_behs] chol_indi_dims_sexspec2; // cholesky factor of correlation matrix for individual-level dims MALES

  cholesky_factor_corr[n_behs] chol_dyad_dims_combispec1; // cholesky factor of correlation matrix for dyad-level dims FF dyads
  cholesky_factor_corr[n_behs] chol_dyad_dims_combispec2; // cholesky factor of correlation matrix for dyad-level dims mix dyads
  cholesky_factor_corr[n_behs] chol_dyad_dims_combispec3; // cholesky factor of correlation matrix for dyad-level dims MM dyads

}

transformed parameters {
  matrix[n_ids, n_behs] indivals;
  matrix[n_dyads, n_behs] dyadvals;

  for (i in 1:n_behs) {
    indivals[, i] = to_vector(indivals_z[i, ]) .* to_vector(sigma_sexspec[i, sex_for_id]);
    dyadvals[, i] = to_vector(dyadvals_z[i, ] .* sigma_combispec[i, combi_for_dyad]);
  }

  indivals[index_female_ids, ] = transpose(diag_pre_multiply(sigma_sexspec[, 1], chol_indi_dims_sexspec1) *
                                       indivals_z[, index_female_ids]);
  indivals[index_male_ids, ] = transpose(diag_pre_multiply(sigma_sexspec[, 2], chol_indi_dims_sexspec2) *
                                       indivals_z[, index_male_ids]);

  dyadvals[index_ff_dyads, ] = transpose(diag_pre_multiply(sigma_combispec[, 1], chol_dyad_dims_combispec1) *
                                           dyadvals_z[, index_ff_dyads]);
  dyadvals[index_mix_dyads, ] = transpose(diag_pre_multiply(sigma_combispec[, 2], chol_dyad_dims_combispec2) *
                                           dyadvals_z[, index_mix_dyads]);
  dyadvals[index_mm_dyads, ] = transpose(diag_pre_multiply(sigma_combispec[, 3], chol_dyad_dims_combispec3) *
                                           dyadvals_z[, index_mm_dyads]);

  for (i in 1:n_behs) {
    indivals[, i] = indivals[, i] + to_vector(icpts_sexspec[i, sex_for_id]) ;
    dyadvals[, i] = dyadvals[, i] + to_vector(icpts_combispec[i, combi_for_dyad]);
    // kin effect
    dyadvals[, i] = dyadvals[, i] + (relatedness_deg_std .* b_relatedness_deg[i]) + (dyad_is_spouse_std .* b_spouse[i]);
  }

}
model {
  vector[n_obs] lp1 = rep_vector(0.0, n_obs); // food
  vector[n_obs] lp2 = rep_vector(0.0, n_obs); // assoc
  vector[n_obs] lp3 = rep_vector(0.0, n_obs); // help
  vector[n_obs] lp4 = rep_vector(0.0, n_obs); // friend

  lp1 = lp1 + b[1];
  lp1 = lp1 + sqrt(0.5) * (indivals[id1_for_obs, 1] + indivals[id2_for_obs, 1]) + dyadvals[dyad_for_obs, 1];
  food ~ binomial_logit(obseff, lp1);

  lp2 = lp2 + b[2];
  lp2 = lp2 + sqrt(0.5) * (indivals[id1_for_obs, 2] + indivals[id2_for_obs, 2]) + dyadvals[dyad_for_obs, 2];
  assoc ~ binomial_logit(obseff, lp2);

  lp3 = lp3 + b[3];
  lp3 = lp3 + sqrt(0.5) * (indivals[id1_for_obs, 3] + indivals[id2_for_obs, 3]) + dyadvals[dyad_for_obs, 3];
  help ~ binomial_logit(obseff, lp3);

  lp4 = lp4 + b[4];
  lp4 = lp4 + sqrt(0.5) * (indivals[id1_for_obs, 4] + indivals[id2_for_obs, 4]) + dyadvals[dyad_for_obs, 4];
  friendship ~ binomial_logit(obseff, lp4);

  // priors
  for (i in 1:n_behs) {
    sigma_sexspec[i, ] ~ student_t(3, 0, 1);
    icpts_sexspec[i, ] ~ student_t(3, 0, 1);
    indivals_z[i, ] ~ normal(0, 1);

    sigma_combispec[i, ] ~ student_t(3, 0, 1);
    icpts_combispec[i, ] ~ student_t(3, 0, 1);
    dyadvals_z[i, ] ~ normal(0, 1);
  }

  b ~ student_t(3, -3, 2.5);
  b_relatedness_deg ~ student_t(3, 0, 1);
  b_spouse ~ student_t(3, 0, 1);

  chol_indi_dims_sexspec1 ~ lkj_corr_cholesky(4);
  chol_indi_dims_sexspec2 ~ lkj_corr_cholesky(4);
  chol_dyad_dims_combispec1 ~ lkj_corr_cholesky(4);
  chol_dyad_dims_combispec2 ~ lkj_corr_cholesky(4);
  chol_dyad_dims_combispec3 ~ lkj_corr_cholesky(4);
}

generated quantities {
  array[n_obs] int<lower=0> food_pred = rep_array(0, n_obs);
  array[n_obs] int<lower=0> assoc_pred = rep_array(0, n_obs);
  array[n_obs] int<lower=0> help_pred = rep_array(0, n_obs);
  array[n_obs] int<lower=0> friendship_pred = rep_array(0, n_obs);

  // extract correlations
  vector<lower=-1,upper=1>[n_cors] cors_indi_female = chol2corvec(chol_indi_dims_sexspec1, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_indi_male = chol2corvec(chol_indi_dims_sexspec2, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_dyad_femfem = chol2corvec(chol_dyad_dims_combispec1, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_dyad_mix = chol2corvec(chol_dyad_dims_combispec2, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_dyad_malemale = chol2corvec(chol_dyad_dims_combispec3, n_behs, n_cors);

  // predicted values
  {
    for (i in 1:4) {
      vector[n_obs] lp_pred = rep_vector(0.0, n_obs); // b_agg[1]

      lp_pred = lp_pred + b[i];
      lp_pred = lp_pred + sqrt(0.5) * (to_vector(indivals[id1_for_obs, i]) + to_vector(indivals[id2_for_obs, i])) + dyadvals[dyad_for_obs, i];

      if (i == 1) {
        food_pred = binomial_rng(obseff, inv_logit(lp_pred));
      }
      if (i == 2) {
        assoc_pred = binomial_rng(obseff, inv_logit(lp_pred));
      }
      if (i == 3) {
        help_pred = binomial_rng(obseff, inv_logit(lp_pred));
      }
      if (i == 4) {
        friendship_pred = binomial_rng(obseff, inv_logit(lp_pred));
      }
    }
  }
}

