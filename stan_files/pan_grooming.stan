data {
  int n_obs;
  int n_indi;
  int n_groups;
  array[n_obs] int<lower=0,upper=n_indi> index_indi;
  array[n_obs] int<lower=0,upper=n_groups> index_group;
  array[n_obs] int<lower=1,upper=4> sexspec_indicator;
  array[n_obs] int<lower=0> partners_obs;
  array[n_obs] int<lower=0> ntrials;
  vector[n_obs] groupsize;
  vector<lower=0>[n_obs] obseff;

  // for predictions
  int n_preds;
  array[n_preds] int<lower=0,upper=n_groups> pred_group_index;
  array[n_preds] int<lower=0,upper=4> pred_sexspec_index;
  array[n_preds] int<lower=0> pred_size;
  real pred_obseff;
}

parameters {
  vector[4] b_raw;
  real<lower=0> indi_sd;
  vector[n_indi] indivals_z;
  real<lower=0> group_sd;
  vector[n_groups] groupvals_z;
  vector[4] gs_slope;
}

transformed parameters {
  vector[4] b = exp(b_raw);
  vector<lower=0,upper=1>[n_obs] probvec;
  vector[n_obs] lp;
  vector[n_indi] indivals = indivals_z * indi_sd;
  vector[n_groups] groupvals = groupvals_z * group_sd;
  lp = b_raw[sexspec_indicator] +
    indivals[index_indi] +
    groupvals[index_group] +
    gs_slope[sexspec_indicator] .* groupsize
    ;
  for (i in 1:n_obs) {
    probvec[i] = 1 - exp(-(exp(lp[i]) * obseff[i])/ntrials[i]);
  }
}

model {
  partners_obs ~ binomial(ntrials, probvec);
  // priors
  b_raw ~ normal(0, 1);
  indi_sd ~ exponential(1);
  indivals_z ~ normal(0, 1);
  group_sd ~ exponential(1);
  groupvals_z ~ normal(0, 1);
  gs_slope ~ normal(0, 0.2);
}

generated quantities {
  vector<lower=0,upper=1>[n_preds] pred_prob;
  vector<lower=0,upper=1>[n_obs] indi_preds;
  {
    real bereal = 0.0;
    for (i in 1:n_preds) {
      bereal = exp(b_raw[pred_sexspec_index[i]] +
        groupvals[pred_group_index[i]] +
        gs_slope[pred_sexspec_index[i]] * pred_size[i]
        );
      pred_prob[i] = 1 - exp(-(bereal * pred_obseff)/(pred_size[i] - 1));
    }

    vector[n_obs] lp_pred = b_raw[sexspec_indicator] + indivals[index_indi] +
      groupvals[index_group] +
      gs_slope[sexspec_indicator] .* groupsize
      ;
    for (i in 1:n_obs) {
      indi_preds[i] = binomial_rng(ntrials[i], 1 - exp(-(exp(lp_pred[i]) * obseff[i])/ntrials[i])) / (1.0 * ntrials[i]);
    }
  }
}
