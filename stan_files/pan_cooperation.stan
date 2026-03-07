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
  int<lower=0> n_grps;
  int<lower=0> n_years;
  int<lower=0> n_pops;

  // indices per observation
  array[n_obs] int<lower=0> id1_for_obs;
  array[n_obs] int<lower=0> id2_for_obs;
  array[n_obs] int<lower=0> dyad_for_obs;
  array[n_obs] int<lower=0> group_for_obs;
  array[n_obs] int<lower=0> year_for_obs;
  array[n_obs] int<lower=0> pop_for_obs;

  // other indices
  array[n_ids] int<lower=0> specsex_for_id;
  array[n_dyads] int<lower=0> combi_for_dyad;
  vector<lower=0,upper=1>[n_dyads] is_kin;

  // indices to map sexspec to individuals
  int<lower=0> n_sexspec1_ids;
  array[n_sexspec1_ids] int<lower=0> sexspec1_ids;
  int<lower=0> n_sexspec2_ids;
  array[n_sexspec2_ids] int<lower=0> sexspec2_ids;
  int<lower=0> n_sexspec3_ids;
  array[n_sexspec3_ids] int<lower=0> sexspec3_ids;
  int<lower=0> n_sexspec4_ids;
  array[n_sexspec4_ids] int<lower=0> sexspec4_ids;
  // indices to map combispec to dyads
  int<lower=0> n_combispec1_dyads;
  array[n_combispec1_dyads] int<lower=0> combispec1_dyads;
  int<lower=0> n_combispec2_dyads;
  array[n_combispec2_dyads] int<lower=0> combispec2_dyads;
  int<lower=0> n_combispec3_dyads;
  array[n_combispec3_dyads] int<lower=0> combispec3_dyads;
  int<lower=0> n_combispec4_dyads;
  array[n_combispec4_dyads] int<lower=0> combispec4_dyads;
  int<lower=0> n_combispec5_dyads;
  array[n_combispec5_dyads] int<lower=0> combispec5_dyads;
  int<lower=0> n_combispec6_dyads;
  array[n_combispec6_dyads] int<lower=0> combispec6_dyads;

  // interaction data
  vector<lower=0,upper=1>[n_obs] grooming;
  array[n_obs] int<lower=0> agg;
  vector<lower=0>[n_obs] obseff;
  array[n_obs] int<lower=0> coal;
  array[n_obs] int<lower=0> foodshare;
  vector<lower=0>[n_obs] obseff_adlib;
}
transformed data {
  // centered version of kin predictor
  vector[n_dyads] is_kin_std = is_kin - mean(is_kin);
}
parameters {
  matrix<lower=0>[n_behs, 4] sigma_sexspec; // [X behaviours, 4 sex-spec combis]
  matrix[n_behs, 4] icpts_sexspec;
  matrix[n_behs, n_ids] indivals_z;

  matrix<lower=0>[n_behs, 6] sigma_combispec; // [X behaviours, 6 sexcombi-spec combis]
  matrix[n_behs, 6] icpts_combispec;
  matrix[n_behs, n_dyads] dyadvals_z;

  vector<lower=0>[n_behs] group_sd;
  matrix[n_grps, n_behs] groupvals_z;

  vector<lower=0>[n_behs] pop_sd;
  matrix[n_pops, n_behs] popvals_z;

  real<lower=0> year_sd;
  vector[n_years] yearvals_z;

  vector[n_behs] b; // 1=groom,2=aggr,3=coal,4=foodshare
  matrix[n_behs, 6] b_kin; // estimate kin effect for dyads per sex combi and species

  // shape parameter for beta (grooming)
  // unconstrained
  vector[n_years] shape_beta_groom_log;

  cholesky_factor_corr[n_behs] chol_indi_dims_sexspec1; // cholesky factor of correlation matrix for individual-level dims
  cholesky_factor_corr[n_behs] chol_indi_dims_sexspec2; // cholesky factor of correlation matrix for individual-level dims
  cholesky_factor_corr[n_behs] chol_indi_dims_sexspec3; // cholesky factor of correlation matrix for individual-level dims
  cholesky_factor_corr[n_behs] chol_indi_dims_sexspec4; // cholesky factor of correlation matrix for individual-level dims

  cholesky_factor_corr[n_behs] chol_dyad_dims_combispec1; // cholesky factor of correlation matrix for dyad-level dims
  cholesky_factor_corr[n_behs] chol_dyad_dims_combispec2; // cholesky factor of correlation matrix for dyad-level dims
  cholesky_factor_corr[n_behs] chol_dyad_dims_combispec3; // cholesky factor of correlation matrix for dyad-level dims
  cholesky_factor_corr[n_behs] chol_dyad_dims_combispec4; // cholesky factor of correlation matrix for dyad-level dims
  cholesky_factor_corr[n_behs] chol_dyad_dims_combispec5; // cholesky factor of correlation matrix for dyad-level dims
  cholesky_factor_corr[n_behs] chol_dyad_dims_combispec6; // cholesky factor of correlation matrix for dyad-level dims

}
transformed parameters {
  matrix[n_ids, n_behs] indivals;
  matrix[n_dyads, n_behs] dyadvals;
  matrix[n_grps, n_behs] groupvals;
  matrix[n_pops, n_behs] popvals;
  vector[n_years] yearvals = rep_vector(0.0, n_years);

  for (i in 1:n_behs) {
    // bring SDs to actual scale
    indivals[, i] = to_vector(indivals_z[i, ]) .* to_vector(sigma_sexspec[i, specsex_for_id]);
    dyadvals[, i] = to_vector(dyadvals_z[i, ] .* sigma_combispec[i, combi_for_dyad]);
  }


  indivals[sexspec1_ids, ] = transpose(diag_pre_multiply(sigma_sexspec[, 1], chol_indi_dims_sexspec1) *
                                       indivals_z[, sexspec1_ids]);
  indivals[sexspec2_ids, ] = transpose(diag_pre_multiply(sigma_sexspec[, 2], chol_indi_dims_sexspec2) *
                                       indivals_z[, sexspec2_ids]);
  indivals[sexspec3_ids, ] = transpose(diag_pre_multiply(sigma_sexspec[, 3], chol_indi_dims_sexspec3) *
                                       indivals_z[, sexspec3_ids]);
  indivals[sexspec4_ids, ] = transpose(diag_pre_multiply(sigma_sexspec[, 4], chol_indi_dims_sexspec4) *
                                       indivals_z[, sexspec4_ids]);

  dyadvals[combispec1_dyads, ] = transpose(diag_pre_multiply(sigma_combispec[, 1], chol_dyad_dims_combispec1) *
                                           dyadvals_z[, combispec1_dyads]);
  dyadvals[combispec2_dyads, ] = transpose(diag_pre_multiply(sigma_combispec[, 2], chol_dyad_dims_combispec2) *
                                           dyadvals_z[, combispec2_dyads]);
  dyadvals[combispec3_dyads, ] = transpose(diag_pre_multiply(sigma_combispec[, 3], chol_dyad_dims_combispec3) *
                                           dyadvals_z[, combispec3_dyads]);
  dyadvals[combispec4_dyads, ] = transpose(diag_pre_multiply(sigma_combispec[, 4], chol_dyad_dims_combispec4) *
                                           dyadvals_z[, combispec4_dyads]);
  dyadvals[combispec5_dyads, ] = transpose(diag_pre_multiply(sigma_combispec[, 5], chol_dyad_dims_combispec5) *
                                           dyadvals_z[, combispec5_dyads]);
  dyadvals[combispec6_dyads, ] = transpose(diag_pre_multiply(sigma_combispec[, 6], chol_dyad_dims_combispec6) *
                                           dyadvals_z[, combispec6_dyads]);

  for (i in 1:n_behs) {
    indivals[, i] = indivals[, i] + to_vector(icpts_sexspec[i, specsex_for_id]) ;
    dyadvals[, i] = dyadvals[, i] + to_vector(icpts_combispec[i, combi_for_dyad]);
    // kin effect
    dyadvals[, i] = dyadvals[, i] + (is_kin_std .* to_vector(b_kin[i, combi_for_dyad]));
    // group and pop effs
    groupvals[, i] = groupvals_z[, i] .* group_sd[i];
    popvals[, i] = popvals_z[, i] .* pop_sd[i];
  }

  yearvals = yearvals_z .* year_sd;
}
model {
  vector[n_obs] lp_gro = rep_vector(0.0, n_obs); // grooming
  vector[n_obs] lp_agg = rep_vector(0.0, n_obs); // aggression
  vector[n_obs] lp_coa = rep_vector(0.0, n_obs); // coalitions
  vector[n_obs] lp_foo = rep_vector(0.0, n_obs); // food sharing
  vector[n_obs] shapetemp = rep_vector(0.0, n_obs);
  int bindex = 0;

  // grooming
  lp_gro = lp_gro + b[1];
  lp_gro = lp_gro + groupvals[group_for_obs, 1] + popvals[pop_for_obs, 1] + yearvals[year_for_obs];
  lp_gro = lp_gro + sqrt(0.5) * (indivals[id1_for_obs, 1] + indivals[id2_for_obs, 1]) + dyadvals[dyad_for_obs, 1];
  shapetemp = exp(shape_beta_groom_log[year_for_obs]); //
  grooming ~ beta(inv_logit(lp_gro) .* shapetemp, (1 - inv_logit(lp_gro)) .* shapetemp);

  // aggression
  bindex = 2;
  lp_agg = lp_agg + log(obseff);
  lp_agg = lp_agg + b[bindex];
  lp_agg = lp_agg + groupvals[group_for_obs, bindex] + popvals[pop_for_obs, bindex] + yearvals[year_for_obs];
  lp_agg = lp_agg + sqrt(0.5) * (indivals[id1_for_obs, bindex] + indivals[id2_for_obs, bindex]) + dyadvals[dyad_for_obs, bindex];
  agg ~ poisson(exp(lp_agg));

  bindex = 3;
  lp_coa = lp_coa + log(obseff_adlib);
  lp_coa = lp_coa + b[bindex];
  lp_coa = lp_coa + groupvals[group_for_obs, bindex] + popvals[pop_for_obs, bindex] + yearvals[year_for_obs];
  lp_coa = lp_coa + sqrt(0.5) * (indivals[id1_for_obs, bindex] + indivals[id2_for_obs, bindex]) + dyadvals[dyad_for_obs, bindex];
  coal ~ poisson(exp(lp_coa));

  bindex = 4;
  lp_foo = lp_foo + log(obseff_adlib);
  lp_foo = lp_foo + b[bindex];
  lp_foo = lp_foo + groupvals[group_for_obs, bindex] + popvals[pop_for_obs, bindex] + yearvals[year_for_obs];
  lp_foo = lp_foo + sqrt(0.5) * (indivals[id1_for_obs, bindex] + indivals[id2_for_obs, bindex]) + dyadvals[dyad_for_obs, bindex];
  foodshare ~ poisson(exp(lp_foo));

  // priors
  for (i in 1:n_behs) {
    sigma_sexspec[i, ] ~ student_t(3, 0, 1);
    icpts_sexspec[i, ] ~ student_t(3, 0, 1);
    indivals_z[i, ] ~ normal(0, 1);

    sigma_combispec[i, ] ~ student_t(3, 0, 1);
    icpts_combispec[i, ] ~ student_t(3, 0, 1);
    dyadvals_z[i, ] ~ normal(0, 1);
    groupvals_z[, i] ~ normal(0, 1);
    popvals_z[, i] ~ normal(0, 1);
    b_kin[i, ] ~ student_t(3, 0, 1);
  }

  yearvals_z ~ normal(0, 1);
  group_sd ~ student_t(3, 0, 1);
  pop_sd ~ student_t(3, 0, 1);
  year_sd ~ student_t(3, 0, 1);
  b ~ student_t(3, -2, 2.5);
  shape_beta_groom_log ~ normal(0, 1);

  chol_indi_dims_sexspec1 ~ lkj_corr_cholesky(4);
  chol_indi_dims_sexspec2 ~ lkj_corr_cholesky(4);
  chol_indi_dims_sexspec3 ~ lkj_corr_cholesky(4);
  chol_indi_dims_sexspec4 ~ lkj_corr_cholesky(4);
  chol_dyad_dims_combispec1 ~ lkj_corr_cholesky(4);
  chol_dyad_dims_combispec2 ~ lkj_corr_cholesky(4);
  chol_dyad_dims_combispec3 ~ lkj_corr_cholesky(4);
  chol_dyad_dims_combispec4 ~ lkj_corr_cholesky(4);
  chol_dyad_dims_combispec5 ~ lkj_corr_cholesky(4);
  chol_dyad_dims_combispec6 ~ lkj_corr_cholesky(4);

}
generated quantities {
  array[n_obs] int<lower=0> agg_pred = rep_array(0, n_obs);
  vector<lower=0,upper=1>[n_obs] groom_pred = rep_vector(0, n_obs);
  vector<lower=0,upper=1>[n_obs] prox_pred = rep_vector(0, n_obs);
  array[n_obs] int<lower=0> coal_pred = rep_array(0, n_obs);
  array[n_obs] int<lower=0> foodshare_pred = rep_array(0, n_obs);

  // extract correlations
  vector<lower=-1,upper=1>[n_cors] cors_indi_sexspec1 = chol2corvec(chol_indi_dims_sexspec1, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_indi_sexspec2 = chol2corvec(chol_indi_dims_sexspec2, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_indi_sexspec3 = chol2corvec(chol_indi_dims_sexspec3, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_indi_sexspec4 = chol2corvec(chol_indi_dims_sexspec4, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_dyad_combispec1 = chol2corvec(chol_dyad_dims_combispec1, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_dyad_combispec2 = chol2corvec(chol_dyad_dims_combispec2, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_dyad_combispec3 = chol2corvec(chol_dyad_dims_combispec3, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_dyad_combispec4 = chol2corvec(chol_dyad_dims_combispec4, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_dyad_combispec5 = chol2corvec(chol_dyad_dims_combispec5, n_behs, n_cors);
  vector<lower=-1,upper=1>[n_cors] cors_dyad_combispec6 = chol2corvec(chol_dyad_dims_combispec6, n_behs, n_cors);

  // attempt to be explicit in naming output vectors:
  real cors_indi_bon_fem_agg_gro;
  real cors_indi_bon_fem_coa_gro;
  real cors_indi_bon_fem_agg_coa;
  real cors_indi_bon_fem_food_gro;
  real cors_indi_bon_fem_agg_food;
  real cors_indi_bon_fem_coa_food;

  {
    vector[n_cors] cors_indi_sexspec1_temp = chol2corvec(chol_indi_dims_sexspec1, n_behs, n_cors);
    cors_indi_bon_fem_agg_gro = cors_indi_sexspec1_temp[1];
    cors_indi_bon_fem_coa_gro = cors_indi_sexspec1_temp[2];
    cors_indi_bon_fem_agg_coa = cors_indi_sexspec1_temp[3];
    cors_indi_bon_fem_food_gro = cors_indi_sexspec1_temp[4];
    cors_indi_bon_fem_agg_food = cors_indi_sexspec1_temp[5];
    cors_indi_bon_fem_coa_food = cors_indi_sexspec1_temp[6];
  }

  // predicted values
  {
    for (i in 1:4) {
      vector[n_obs] lp_pred = rep_vector(0.0, n_obs); // b_agg[1]
      if (i == 2) { // aggression
        lp_pred = lp_pred + log(obseff);
      }
      if (i == 3) {
        lp_pred = lp_pred + log(obseff_adlib);
      }
      if (i == 4) {
        lp_pred = lp_pred + log(obseff_adlib);
      }

      lp_pred = lp_pred + b[i];
      lp_pred = lp_pred + groupvals[group_for_obs, i] + popvals[pop_for_obs, i] + yearvals[year_for_obs];
      lp_pred = lp_pred + sqrt(0.5) * (to_vector(indivals[id1_for_obs, i]) + to_vector(indivals[id2_for_obs, i])) + dyadvals[dyad_for_obs, i];

      if (i == 1) {
        vector[n_obs] shapetemp = exp(shape_beta_groom_log[year_for_obs]);
        groom_pred = to_vector(beta_rng(shapetemp .* (inv_logit(lp_pred)),
                                        shapetemp .* (1 - inv_logit(lp_pred))));
      }
      if (i == 2) {
         agg_pred = poisson_rng(exp(lp_pred));
      }
      if (i == 3) {
         coal_pred = poisson_rng(exp(lp_pred));
      }
      if (i == 4) {
         foodshare_pred = poisson_rng(exp(lp_pred));
      }
    }
  }

}
