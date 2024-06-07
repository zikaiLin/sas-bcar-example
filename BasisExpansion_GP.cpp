#include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

class BE_MSE
{
private:
  struct BEData
  {
    int num_subjects;   // n
    int num_predictors; // p
    vec y;              // n*1
    mat X;              // n*p
    cube B;             // n*L*j
    umat partition;     // image partition matrix
    uvec radius_candidates; // Radius candidates
    uvec excluded_vox; // excluded voxel, contains no information
  } dat;
  
  struct BEParas
  {
    // bool spatial_independent;
    vec delta;    // p*1 <=> delta_j 1*1
    mat theta;    // L*p <=>theta_j L*1
    double alpha; // 1*1
    vec z;        // n*1
    mat d_B_th;
    mat B_th;
    double inv_sigma_sq;
    double inv_a;
    double loglik;
    vec loglik_vec;
    double logpost;
    mat prob_mat;
    mat delta_mat;
    uvec radius_post;
    double prior_min;
    double prior_max;
    // mat BTzi;
    // mat z_outj;
    // mat theta_test;
  } paras;
  
  struct MCMCSample
  {
    mat delta;
    vec alpha;
    vec sigma_sq;
    vec a;
    mat z;
    // mat loglik;
    // cube d_B_th;
    // cube B_th;
    // cube BTzi;
    // cube z_outj;
    cube theta;
    // cube theta_test;
  } paras_sample;
  
  struct GibbsSamplerProfile
  {
    vec loglik;
    vec logpost;
    mat loglik_mat;
  } gibbs_profile;
  
  struct GibbsSamplerControl
  {
    int total_iter;
    int burnin;
    int mcmc_sample;
    int thinning;
    int radius_update_thinning;
    int verbose;
    int save_profile;
    int total_profile;
  } gibbs_control;
  
  int iter;
  
public:
  // truncated normal reject sampling for y_i=1
  double slice_truncnorm_positive(double mu, double sigma, double a)
  {
    Environment pkg = Environment::namespace_env("RcppTN");
    Function f = pkg["rtn"];
    double x0 = as<double>(f(mu, sigma, a, INFINITY));
    return (x0);
  }
  
  // truncated normal reject sampling for y_i=0
  double slice_truncnorm_negative(double mu, double sigma, double b)
  {
    Environment pkg = Environment::namespace_env("RcppTN");
    Function f = pkg["rtn"];
    double x0 = as<double>(f(mu, sigma, -INFINITY, b));
    return (x0);
  }
  
  void Hermite_polynomial(const mat &X,
                          const int &poly_degree,
                          const double &poly_a,
                          const double &poly_b)
  {
    dat.B.zeros(dat.num_subjects, poly_degree + 1, dat.num_predictors);
    Environment pkg = Environment::namespace_env("BayesGPfit");
    Function f = pkg["GP.eigen.funcs.fast"];
    for (int j = 0; j < dat.num_predictors; j++)
    {
      mat v_list = X.col(j);
      mat A = as<mat>(f(v_list, poly_degree, poly_a, poly_b));
      dat.B.slice(j) = A;
    }
  }
  
  void load_data(const vec &in_y,
                 const mat &in_X,
                 const int &poly_degree,
                 const double &poly_a,
                 double &poly_b,
                 umat &partition,
                 uvec &radius_candidates,
                 uvec &excluded_vox,
                 double prior_max = 1,
                 double prior_min = 0)
  {
    dat.y = in_y;
    dat.X = in_X;
    dat.num_predictors = dat.X.n_cols;
    dat.num_subjects = dat.X.n_rows;
    dat.partition= partition;
    dat.radius_candidates = radius_candidates;
    dat.excluded_vox = excluded_vox;
    paras.prior_max = prior_max;
    paras.prior_min = prior_min;
    // set initial B
    Hermite_polynomial(in_X, poly_degree, poly_a, poly_b);
  };
  
  void set_paras_initial_values(const double &initial_sigma_sq,
                                const double &initial_a,
                                double initial_alpha, 
                                const vec &initial_delta,
                                vec initial_z,
                                mat initial_theta,
                                mat initial_prob_mat)
  {
    // if (dat.radius_candidates.n_elem == 1 && dat.radius_candidates(0) == 0){
    //   paras.spatial_independent = true;
    // }else{
    //   paras.spatial_independent = false;
    // }
    
    paras.prob_mat = initial_prob_mat;
    // set initial value for sigma_sq
    
    paras.inv_sigma_sq = 1.0 / initial_sigma_sq;
    
    // set initial value for a
    paras.inv_a = 1/initial_a;
    
    // set initial value for theta_j
    paras.theta = initial_theta;
    
    // set initial value for alpha
    paras.alpha = initial_alpha;
    
    // set initial value for delta_j
    paras.delta = initial_delta;
    // set initial value for middle
    
    paras.d_B_th = zeros<mat>(dat.num_subjects,dat.num_predictors);
    
    paras.B_th = zeros<mat>(dat.num_subjects,dat.num_predictors);
    
    uvec unique_partition = unique(dat.partition);
    
    paras.radius_post = ones<uvec>(unique_partition.n_elem);
    
    // Rcout << dat.num_predictors << "\n";
    // paras.theta.brief_print();
    // paras.B_th.brief_print();
    for (int j = 0; j < dat.num_predictors; j++){
      paras.B_th.col(j) = dat.B.slice(j)*paras.theta.col(j);
    }
    //delta_j*B_j.row(i)*theta_j
    for (int j = 0; j < dat.num_predictors; j++){
      paras.d_B_th.col(j) = paras.delta(j)*paras.B_th.col(j);
    }
    // set initial value for z
    paras.z = initial_z;
    // paras.z = randn<vec>(dat.num_subjects) + paras.alpha + sum(paras.d_B_th,1);
    // paras.z = zeros<vec>(dat.num_subjects);
    // update_z();
    // paras.BTzi.zeros(dat.B.n_cols,dat.num_predictors);
    // paras.z_outj.zeros(dat.num_subjects,dat.num_predictors);
    // paras.theta_test.zeros(dat.B.n_cols,dat.num_predictors);


    
  };
  
  
  vec update_pj_mat_partial(int part,int radius){
    
    paras.delta_mat = reshape(paras.delta, paras.prob_mat.n_rows,  paras.prob_mat.n_cols);
    uvec part_indices = find(dat.partition == part);
    umat part_sub = ind2sub(size(dat.partition), part_indices);
    
    int nr = paras.delta_mat.n_rows;
    int nc = paras.delta_mat.n_cols;
    
    vec pj_mat_part(part_indices.n_elem);
    // pj_mat_part.brief_print();
    Rcpp::Function expGrid("expand.grid");
    
    for(int p = 0; p<part_indices.n_elem;p++){
      int i = part_sub(0,p);
      int j = part_sub(1,p);
      
      List ij_circle_list = expGrid(seq(i-radius, i+radius), seq(j-radius, j+radius));
      IntegerVector col1 = ij_circle_list[0];
      IntegerVector col2 = ij_circle_list[1];
      IntegerMatrix ij_circle= cbind(col1, col2);
      // Rcout << ij_circle << "\n";
      int count = 0;
      double sum = 0;
      for(int s=0; s<ij_circle.nrow();s++){
        IntegerVector s_idx = ij_circle.row(s);
        if(s_idx(0) < 0 | s_idx(1) < 0 | s_idx(0) >= nr |  s_idx(1) >= nc){
          continue;
        }else{
          sum += paras.delta_mat(s_idx(0), s_idx(1));
          count += 1;
        }
        
        // Rcout << "Update P_ij where i = " << i << " j = " << j <<"\n";
        // paras.prob_mat(i,j) = sum/count;
        if(sum/count > 0.1){
          pj_mat_part(p) = 1;
        }else{
          pj_mat_part(p) = sum/count;
        }
        pj_mat_part(p) = sum/count;
      }
      
      
    }
    
    return pj_mat_part;
  }
  
  
  
  
  void update_prob_mat(){
    int nr = paras.delta_mat.n_rows;
    int nc = paras.delta_mat.n_cols;
    
    
    Rcpp::Function expGrid("expand.grid");
    
    for(int i=0; i<nr; i++){
      for(int j=0;j<nc; j++){
        int radius = paras.radius_post(dat.partition(i,j)-1);
        
        if (radius == 0){
          paras.prob_mat(i,j) = 0.5;
          continue;
        }
        // Rcout << "Expanding Grids near i=" << i << " j = " << j << "\n";
        List ij_circle_list = expGrid(seq(i-radius, i+radius), seq(j-radius, j+radius));
        IntegerVector col1 = ij_circle_list[0];
        IntegerVector col2 = ij_circle_list[1];
        umat ij_circle= as<umat>(cbind(col1, col2));
        // Rcout << ij_circle << "\n";
        int count = 0;
        double sum = 0;
        
        // Rcout << "Computing Spatial Map Probability"<<"\n";
        for(int s=0; s<ij_circle.n_rows;s++){
          // Rcout << "Extracting row" << s<<"\n";
          urowvec s_idx = ij_circle.row(s);
          
          // Rcout << "Counting" << s<<"\n";
          if(s_idx(0) < 0 | s_idx(1) < 0 | s_idx(0) >= nr |  s_idx(1) >= nc){
            continue;
          }else{
            uword s_idx_uword = sub2ind(size(paras.delta_mat), s_idx(0), s_idx(1));
            if(dat.excluded_vox(s_idx_uword) == 1) continue;
            sum += paras.delta_mat(s_idx(0), s_idx(1));
            count += 1;
          }
        }
        paras.prob_mat(i,j) = sum/count;
      }
      
    }
    
    paras.prob_mat = paras.prob_mat * (paras.prior_max - paras.prior_min) + paras.prior_min;
    
  }
  
  
  void update_radius(){
    IntegerMatrix delta_mat_cpp = Rcpp::as<Rcpp::IntegerMatrix>(wrap(paras.delta_mat));
    uvec partition_unique = unique(dat.partition);
    // std::cout << "Partition_unique = " << partition_unique;  
    mat expectation_mat = paras.prob_mat;
    int n_candidates = dat.radius_candidates.n_elem;
    vec log_likelihood(n_candidates);
    
    
    for(int p=0; p<partition_unique.n_elem; p++){
      
      // Rcout << "Partition: "<<p+1<<"\n";
      uvec voxels_p = find(dat.partition == p+1);
      // voxels_p.brief_print();
      
      for(int i=0; i < n_candidates; i++){
        int radius_i = dat.radius_candidates(i);
        mat expecation_pj_R = expectation_mat;
        
        
        // Rcout << "Updatintg partitial expectation matrix" <<"\n";
        
        vec sub_expectation_pj_R = update_pj_mat_partial(p+1, radius_i);
        
        // sub_expectation_pj_R.brief_print();
        // Rcout << "Copy submat to expectation R\n";
        expecation_pj_R.elem(voxels_p) = sub_expectation_pj_R;
        rowvec expecation_pj_R_vec = expecation_pj_R.as_row();
        mat basis_x = paras.B_th.each_row() % expecation_pj_R_vec;
        
        vec temp = paras.z - paras.alpha - sum(basis_x,1);
        log_likelihood(i) = -0.5*accu(temp%temp);
      }
      
      log_likelihood -= max(log_likelihood);
      log_likelihood = exp(log_likelihood);
      log_likelihood /= accu(log_likelihood);
      
      // Rcout << "Select Radius\n";
      // Rcout << "p = " << p;
      // paras.radius_post.brief_print();
      paras.radius_post(p) = dat.radius_candidates(log_likelihood.index_max());
    }
    
    // Rcout << "Finished Update Radius\n";
    
    
    
  }
  
  
  void update_delta()
  {
    // double log_p_loglik;
    double log_p_loglik0;
    double log_p_loglik1;
    double posterior_p;
    vec z_outj;
    vec x_beta_outj0;
    vec x_beta_outj1;
    colvec prob_vec = paras.prob_mat.as_col();
    // bool printed = false;
    for (int j = 0; j < dat.num_predictors; j++)
    {
      // paras.d_B_th.col(j).brief_print();
      if(dat.excluded_vox(j) == 1){
        paras.delta(j) = 0;
        continue;
      }else{
        
        x_beta_outj0 = paras.alpha + sum(paras.d_B_th,1) - paras.d_B_th.col(j);
        x_beta_outj1 = x_beta_outj0 + paras.B_th.col(j);
        
      
        // std::cout << "x_beta_outj0 = " << x_beta_outj0 << "\n";
        // std::cout << "x_beta_outj1 = " << x_beta_outj1 << "\n";
        vec xbeta_normcdf0 = normcdf(x_beta_outj0, 0, 1);
        vec xbeta_normcdf1 = normcdf(x_beta_outj1, 0, 1);

        log_p_loglik0 = accu(dat.y % log(xbeta_normcdf0 + arma::datum::eps) + (1.0 - dat.y) % log(1-xbeta_normcdf0 + arma::datum::eps));
        log_p_loglik1 = accu(dat.y % log(xbeta_normcdf1 + arma::datum::eps) + (1.0 - dat.y) % log(1-xbeta_normcdf1 + arma::datum::eps));
        
        // std::cout << "log_p_loglik0 = " << log_p_loglik0 << "\n";
        // std::cout << "log_p_loglik1 = " << log_p_loglik1 << "\n";
        // z_outj = paras.z - paras.alpha - sum(paras.d_B_th,1) + paras.d_B_th.col(j);
        // log_p_loglik = (0.5*accu(paras.B_th.col(j)%paras.B_th.col(j))-accu(z_outj%paras.B_th.col(j)));

        log_p_loglik0 += log(1-prob_vec(j) + datum::eps);
        log_p_loglik1 += log(prob_vec(j) + datum::eps);
        // std::cout << "log_p_loglik0 = " << log_p_loglik0 << "\n";
        // std::cout << "log_p_loglik1 = " << log_p_loglik1 << "\n";
        
        if(log_p_loglik0 < log_p_loglik1){
          log_p_loglik0 -= log_p_loglik1;
          log_p_loglik1 = 0;
        }else{
          log_p_loglik1 -= log_p_loglik0;
          log_p_loglik0 = 0;
        }
        
        posterior_p = exp(log_p_loglik1)/(exp(log_p_loglik0)+exp(log_p_loglik1));
        // std::cout << "posterior_p = " << posterior_p << "\n";
        // posterior_p = 1.0/(1+exp(log_p_loglik));
        paras.delta(j) = R::rbinom(1, posterior_p);
        // if(posterior_p>=0.5) paras.delta(j) = 0;
        // else paras.delta(j) = 1;
        
        
        // if(!printed){
        //   x_beta_outj0.brief_print();
        //   x_beta_outj1.brief_print();
        //   std::cout << paras.alpha << "\n";
        //   paras.B_th.col(j).brief_print();
        //   std::cout << "log_p_loglik0 = " << log_p_loglik0 << "\n";
        //   std::cout << "log_p_loglik1 = " << log_p_loglik1 << "\n";
        //   std::cout << "posterior_p = " << posterior_p << "\n";
        //   printed = true;
        // }
        paras.d_B_th.col(j) = paras.delta(j)*paras.B_th.col(j);
      }

    }
  }
  
  // 
  // void update_delta(){
  //   double prob_1;
  //   double prob_0;
  //   
  //   int n = gp_cur.n_rows;
  //   int p = gp_cur.n_cols;
  //   arma::mat probs(2,p);
  //   
  //   
  //   for(int j = 0; j<p; j++){
  //     arma::uvec delta_j = delta;
  //     // delta_j.brief_print();
  //     delta_j(j) = 0;
  //     
  //     arma::mat gp_cur_delta = gp_cur.cols(arma::find(delta_j == 1));
  //     arma::vec sum_mean = y_prime - alpha - arma::sum(gp_cur_delta, 1);
  //     // sum_mean.brief_print();
  //     // prob_0 = pow(1-log(pj_vec(j)), 1/n) + arma::mean(dnormLog_vec(sum_mean,0,1.0));
  //     vec dloglik = dnormLog_vec(sum_mean,0,1.0);
  //     prob_0 = log(1-pj_vec(j)+datum::eps)+ arma::sum(dloglik);
  //     // prob_0 = arma::sum(dnormLog_vec(sum_mean,0,1.0));
  //     
  //     sum_mean -= gp_cur.col(j);
  //     // gp_cur.col(j).brief_print();
  //     // sum_mean.brief_print();
  //     // prob_1 = pow( log(pj_vec(j)), 1/n)+ arma::mean(dnormLog_vec(sum_mean,0,1.0));
  //     prob_1 = log(pj_vec(j) + datum::eps) + arma::sum(dnormLog_vec(sum_mean,0,1.0));
  //     // prob_1 = arma::sum(dnormLog_vec(sum_mean,0,1.0));
  //     
  //     
  //     probs(0,j) = prob_0;
  //     probs(1,j) = prob_1;
  //   }
  //   
  //   paras.prob_mat = probs;  
  // }
  // 
  
  void set_gibbs_control(int in_mcmc_sample, int in_burnin, int in_thinning, 
                         int in_verbose, int in_save_profile,int in_radius_update_thinning)
  {
    gibbs_control.mcmc_sample = in_mcmc_sample;
    gibbs_control.burnin = in_burnin;
    gibbs_control.thinning = in_thinning;
    gibbs_control.total_iter = gibbs_control.burnin;
    gibbs_control.radius_update_thinning = in_radius_update_thinning,
    gibbs_control.total_iter += gibbs_control.mcmc_sample * gibbs_control.thinning;
    gibbs_control.verbose = in_verbose;
    gibbs_control.save_profile = in_save_profile;
    if (gibbs_control.save_profile > 0)
    {
      gibbs_control.total_profile = gibbs_control.total_iter / gibbs_control.save_profile;
    }
    else
    {
      gibbs_control.total_profile = 0;
    }
  };
  
  vec dnormLog_vec(vec x, double mu, double precision) {
    int n = x.n_elem;
    vec res(n);
    for(int i = 0; i < n; i++) {
      res(i) = R::dnorm4(x(i), mu, std::sqrt(1/precision), 1);
    }
    return res;
  }
  
  
  void update_inv_sigma_sq()
  {
    paras.inv_sigma_sq =
      randg(distr_param(0.5 * dat.num_predictors * dat.B.n_cols + 1,
                        1.0 / (paras.inv_a +
                          0.5 * accu(paras.theta % paras.theta) + 0.5 * paras.alpha * paras.alpha)));
  };
  
  void update_inv_a()
  {
    paras.inv_a = randg(distr_param(1.0, 1.0 / (paras.inv_sigma_sq + 1)));
  };
  
  void update_alpha()
  {
    double temp = randn<double>();
    double mu = sum(paras.z - sum(paras.d_B_th,1));
    paras.alpha = (temp + mu*sqrt(1 / (dat.num_subjects + paras.inv_sigma_sq)))
      * sqrt(1 / (dat.num_subjects + paras.inv_sigma_sq));
  }
  
  void update_z()
  {
    vec temp =  sum(paras.d_B_th,1);
    for (int i = 0; i < dat.num_subjects; i++)
    {
      if (dat.y(i) < 1)
      {
        paras.z(i) = slice_truncnorm_negative(paras.alpha + temp(i),
                1, 0);
      }
      else
      {
        paras.z(i) = slice_truncnorm_positive(paras.alpha + temp(i),
                1, 0);
      }
    }
  }
  
  void update_theta()
  {
    for (int j = 0; j < dat.num_predictors; j++)
    {
      if (dat.excluded_vox(j) == 1) continue;
      if (paras.delta(j) == 1)
      {
        mat BiTBi = zeros<mat>(dat.B.n_cols, dat.B.n_cols);
        for (int i = 0; i < dat.num_subjects ; i++)
        {
          BiTBi += dat.B.slice(j).row(i).t()
          * dat.B.slice(j).row(i);
        }
        BiTBi.diag() += paras.inv_sigma_sq;
        mat V_kj;
        vec eigval;
        eig_sym(eigval, V_kj, BiTBi);
        mat D_kj = zeros<mat>(eigval.size(), eigval.size());
        D_kj.diag() += eigval;
        mat half_matrix = D_kj;
        half_matrix.diag() = sqrt(1 / half_matrix.diag());
        vec BTzi = zeros<vec>(dat.B.n_cols);
        vec z_outj = paras.z - paras.alpha - sum(paras.d_B_th,1) + paras.d_B_th.col(j);
        for(int i=0;i<dat.num_subjects;i++){
          BTzi += dat.B.slice(j).row(i).t()*z_outj(i);
        }
        vec temp = randn<vec>(dat.B.n_cols);
        paras.theta.col(j) =
          V_kj * half_matrix * (temp + half_matrix * V_kj.t() * BTzi);
        // Environment pkg = Environment::namespace_env("MASS");
        // Function f = pkg["mvrnorm"];
        // paras.theta.col(j) = as<vec>(f(1,BiTBi.i()*BTzi , BiTBi.i()));
        // paras.BTzi.col(j) = BTzi;
        // paras.z_outj.col(j) = z_outj;
      }
      else{
        paras.theta.col(j) =
          randn<vec>(dat.B.n_cols)*sqrt(1/paras.inv_sigma_sq);
      }
      paras.B_th.col(j) = dat.B.slice(j)*paras.theta.col(j);
      paras.d_B_th.col(j) = paras.delta(j)*paras.B_th.col(j);
    }
    // paras.theta_test =  paras.theta;
  }
  
  
  
  
  void update_loglik()
  {
    vec xbeta =  paras.alpha + sum(paras.d_B_th,1);
    // Rcpp::NumericVector xbeta_rcpp = NumericVector(xbeta.begin(), xbeta.end());
    
    // vec loglik_vec = zeros(dat.y.n_elem);
    vec xbeta_normcdf = normcdf(xbeta, 0, 1);
    // xbeta_normcdf.brief_print();
    // loglik_vec = dat.y % log(normcdf(xbeta, 0, 1)) + (1 - dat.y) % log(1-normcdf(xbeta,0,1));
    paras.loglik_vec = dat.y % log(xbeta_normcdf + arma::datum::eps) + (1.0 - dat.y) % log(1-xbeta_normcdf + arma::datum::eps);
    paras.loglik = accu(paras.loglik_vec);
  }
  
  void update_logpost()
  {
    paras.logpost = 2 * log(paras.inv_a) - paras.inv_a
    + 2 * log(paras.inv_sigma_sq)
    - paras.inv_a * paras.inv_sigma_sq - 0.5*paras.alpha*paras.alpha * paras.inv_sigma_sq ;
    double inv_sigma_2L = pow(paras.inv_sigma_sq,dat.B.n_cols);
    paras.logpost += 0.5*log(inv_sigma_2L)*dat.num_predictors
      - 0.5 * accu(paras.theta % paras.theta)*paras.inv_sigma_sq;
    paras.logpost += paras.loglik;
  }
  
  void initialize_paras_sample()
  {
    paras_sample.z.zeros(dat.num_subjects, gibbs_control.mcmc_sample);
    // paras_sample.d_B_th.zeros(dat.num_subjects, dat.num_predictors, gibbs_control.mcmc_sample);
    // paras_sample.B_th.zeros(dat.num_subjects, dat.num_predictors, gibbs_control.mcmc_sample);
    // paras_sample.BTzi.zeros(dat.B.n_cols, dat.num_predictors, gibbs_control.mcmc_sample);
    // paras_sample.z_outj.zeros(dat.num_subjects, dat.num_predictors, gibbs_control.mcmc_sample);
    // paras_sample.theta_test.zeros(dat.B.n_cols,dat.num_predictors, gibbs_control.mcmc_sample);
    paras_sample.delta.zeros(dat.num_predictors, gibbs_control.mcmc_sample);
    paras_sample.sigma_sq.zeros(gibbs_control.mcmc_sample);
    paras_sample.a.zeros(gibbs_control.mcmc_sample);
    paras_sample.alpha.zeros(gibbs_control.mcmc_sample);
    paras_sample.theta.zeros(dat.B.n_cols,gibbs_control.mcmc_sample,dat.num_predictors);
    // paras_sample.loglik.zeros(gibbs_control.mcmc_sample, dat.num_subjects);
    
  }
  
  void save_paras_sample()
  {
    if (iter > gibbs_control.burnin)
    {
      if ((iter - gibbs_control.burnin) % gibbs_control.thinning == 0)
      {
        int mcmc_iter = (iter - gibbs_control.burnin) / gibbs_control.thinning;
        paras_sample.z.col(mcmc_iter) = paras.z;
        // paras_sample.d_B_th.slice(mcmc_iter) = paras.d_B_th;
        // paras_sample.B_th.slice(mcmc_iter) = paras.B_th;
        // paras_sample.BTzi.slice(mcmc_iter) = paras.BTzi;
        // paras_sample.theta_test.slice(mcmc_iter) = paras.theta_test;
        // paras_sample.z_outj.slice(mcmc_iter) = paras.z_outj;
        paras_sample.delta.col(mcmc_iter) = paras.delta;
        paras_sample.sigma_sq(mcmc_iter) = 1.0 / paras.inv_sigma_sq;
        paras_sample.a(mcmc_iter) = 1.0 / paras.inv_a;
        paras_sample.alpha(mcmc_iter) = paras.alpha;
        // paras_sample.loglik.row(mcmc_iter) = paras.loglik_vec.as_row();
        
        for(int j=0; j < dat.num_predictors; j++){
          paras_sample.theta.slice(j).col(mcmc_iter) = paras.theta.col(j);
        }
      }
    }
  };
  
  void initialize_gibbs_profile()
  {
    std::cout << "total_profile: " << gibbs_control.total_profile << std::endl;
    if (gibbs_control.save_profile > 0)
    {
      gibbs_profile.loglik.zeros(gibbs_control.total_profile);
      gibbs_profile.loglik_mat.zeros(gibbs_control.total_profile, dat.num_subjects);
      gibbs_profile.logpost.zeros(gibbs_control.total_profile);
    }
  }
  
  void save_gibbs_profile()
  {
    if (gibbs_control.save_profile > 0)
    {
      if (iter % gibbs_control.save_profile == 0)
      {
        int profile_iter = iter / gibbs_control.save_profile;
        update_loglik();
        update_logpost();
        gibbs_profile.loglik(profile_iter) = paras.loglik;
        gibbs_profile.loglik_mat.row(profile_iter) = paras.loglik_vec.as_row();
        gibbs_profile.logpost(profile_iter) = paras.logpost;
      }
    }
  }
  
  void monitor_gibbs()
  {
    if (gibbs_control.verbose > 0)
    {
      if (iter % gibbs_control.verbose == 0)
      {
        std::cout << "iter: " << iter
                  << " inv_sigma_sq: " << paras.inv_sigma_sq
                  << " loglik: " << paras.loglik 
                  << " logpost: " << paras.logpost << std::endl;
      }
    }
  }
  
  void run_gibbs(bool &oc_sigma,
                 bool &oc_a,
                 bool &oc_alpha,
                 bool &oc_delta,
                 bool &oc_theta,
                 bool &oc_z,
                 int poly_degree,
                 double poly_a,
                 double poly_b)
  {
    
    std::cout << "Initializing Parameters" << std::endl;
    initialize_paras_sample();
    
    std::cout << "Initializing Gibbs Sampler" << std::endl;
    initialize_gibbs_profile();
    std::cout << "total iter:" << gibbs_control.total_iter << std::endl;
    for (iter = 0; iter < gibbs_control.total_iter; iter++)
    {
      // std::cout << "Iter: "<< iter << " --------------- " << std::endl;
        // std::cout << "Updating Basis Coefficients... \n";
        update_theta();
        // std::cout << "Updating Intercept... \n";
        update_alpha();
  
        
        // std::cout << "Updating Radius... \n";
          update_radius();
        
        
        // std::cout << "Updating Selection Parameters... \n";
        update_delta();
        

        // std::cout << "Updating Spatial inclusion Probability... \n";
        update_prob_mat();
      
        // std::cout << "Updating Latent Outcome Z... \n";
        update_z();
        // std::cout << "Updating Sigma... \n";
        update_inv_sigma_sq();
      // std::cout << "Updating Cauchy Scale... \n";
      update_inv_a();
      save_paras_sample();
      save_gibbs_profile();
      monitor_gibbs();
    }
  };
  
  // output for R
  List get_gibbs_post_mean(double threshold)
  {
    vec z = mean(paras_sample.z, 1);
    vec delta = mean(paras_sample.delta, 1);
    for(int i=0;i<delta.size();i++){
      if(delta(i)>= threshold){
        delta(i) = 1;
      }else{
        delta(i) = 0;
      }
    }
    return List::create(Named("z") = z,
                        Named("delta") = delta,
                        Named("sigma_sq") = mean(paras_sample.sigma_sq),
                        Named("a") = mean(paras_sample.a),
                        Named("alpha") = mean(paras_sample.alpha),
                        Named("theta") = mean(paras_sample.theta,1),
                        Named("radius") = paras.radius_post
    );
  };
  
  List get_gibbs_sample()
  {
    return List::create(Named("z") = paras_sample.z,
                        Named("delta") = paras_sample.delta,
                        Named("sigma_sq") = paras_sample.sigma_sq,
                        // Named("d_B_th") = paras_sample.d_B_th,
                        // Named("B_th") = paras_sample.B_th,
                        // Named("BTzi") = paras_sample.BTzi,
                        // Named("z_outj") = paras_sample.z_outj,
                        // Named("theta_test") = paras_sample.theta_test,
                        Named("a") = paras_sample.a,
                        Named("alpha") = paras_sample.alpha
                        // Named("theta") = paras_sample.theta
    );
  };
  
  List get_gibbs_trace()
  {
    uvec iters = linspace<uvec>(1, gibbs_control.total_iter, gibbs_control.total_profile);
    return List::create(Named("iters") = iters,
                        Named("loglik") = gibbs_profile.loglik,
                        Named("loglik_mat") = gibbs_profile.loglik_mat,
                        Named("logpost") = gibbs_profile.logpost);
  }
  
  List get_gibbs_control()
  {
    return List::create(Named("total_iter") = gibbs_control.total_iter,
                        Named("burnin") = gibbs_control.burnin,
                        Named("mcmc_sample") = gibbs_control.mcmc_sample,
                        Named("thinning") = gibbs_control.thinning,
                        Named("verbose") = gibbs_control.verbose,
                        Named("save_profile") = gibbs_control.save_profile,
                        Named("total_profile") = gibbs_control.total_profile);
  }
  
  int get_iter()
  {
    return iter;
  };
};



//[[Rcpp::export]]
List simul_dat_gam(int n, double a, vec& delta, int poly_degree, 
                   double poly_a, double poly_b) {
  double sigma_sq = 1/(randg(distr_param(0.5, a)));
  mat X = randn<mat>(n, delta.n_elem);
  Environment pkg = Environment::namespace_env("BayesGPfit");
  Function function = pkg["GP.eigen.funcs.fast"];
  cube B = zeros<cube>(n,poly_degree+1,delta.n_elem);
  for (int j = 0; j < delta.n_elem; j++)
  {
    mat v_list = X.col(j);
    mat A = as<mat>(function(v_list, poly_degree, poly_a, poly_b));
    B.slice(j) = A;
  }
  mat theta = zeros<mat>(poly_degree+1, delta.n_elem);
  for (int j = 0; j < delta.n_elem; j++){
    theta.col(j)=randn<vec>(poly_degree+1)*sqrt(sigma_sq) ;
  }
  vec temp;
  vec y = zeros<vec>(n);
  mat f = zeros<mat>(n,delta.n_elem);
  double alpha = randn<double>()*sqrt(sigma_sq);
  vec z = zeros<vec>(n);
  // vec z_1 = zeros<vec>(n);
  mat d_f = zeros<mat>(n,delta.n_elem);
  // for (int j = 0; j < delta.n_elem; j++){
  //   d_f.col(j) =  delta(j)*B.slice(j)*theta.col(j);
  // }
  // z_1 = sum(d_f,1);
  for (int j = 0; j < delta.n_elem; j++){
    f.col(j) = B.slice(j)*theta.col(j);
  }
  for (int i = 0; i < n; i++)
  {
    temp = zeros<vec>(delta.n_elem);
    for (int j = 0; j < delta.n_elem; j++){
      temp(j) = delta(j)*f(i,j);
    }
    z(i) = randn<double>() + alpha + accu(temp);
    if(z(i)>0){
      y(i) = 1;
    }
    else{
      y(i) = 0;
    }
  }
  
  return List::create(Named("y") = y,
                      Named("X") = X,
                      Named("alpha") = alpha,
                      Named("delta") = delta,
                      Named("poly_degree") = poly_degree,
                      Named("poly_a") = poly_a,
                      Named("poly_b") = poly_b,
                      Named("f") = f,
                      Named("theta") = theta,
                      Named("B") = B,
                      Named("a") = a,
                      Named("z") = z,
                      // Named("z_1") = z_1,
                      // Named("d_f") = d_f,
                      Named("sigma_sq") = sigma_sq
  );
}





// [[Rcpp::export]]
List Basis_Expansion(vec &y,
                     mat &X,
                     double &initial_alpha,
                     double &initial_a,
                     vec &initial_delta,
                     vec &initial_z,
                     mat &initial_theta,
                     mat &initial_prob_mat,
                     uvec &excluded_vox,
                     int poly_degree,
                     double poly_a,
                     double poly_b,
                     uvec &radius_candidats,
                     umat &partition,
                     bool oc_theta = true,
                     bool oc_alpha = true,
                     bool oc_delta = true,
                     bool oc_z = true,
                     bool oc_sigma = true,
                     bool oc_a = true,
                     double initial_sigma_sq = 1,
                     double threshold = 0.5,
                     double prior_max = 1,
                     double prior_min = 0,
                     int mcmc_sample = 500,
                     int burnin = 5000,
                     int thinning = 10,
                     int radius_update_thinning = 50,
                     int max_iter = 1000,
                     int verbose = 5000,
                     int save_profile = 1)
{
  std::cout << "Starting Clock... \n";
  wall_clock timer;
  timer.tic();
  BE_MSE model;
  
  std::cout << "Loading Data... \n";
  model.load_data(y, X, 
                  poly_degree,
                  poly_a,
                  poly_b,
                  partition,
                  radius_candidats,
                  excluded_vox,
                  prior_max = prior_max,
                  prior_min = prior_min);
  
  std::cout << "Setting Gibbs Control... \n";
  
  model.set_gibbs_control(mcmc_sample,
                          burnin,
                          thinning,
                          verbose,
                          save_profile,
                          radius_update_thinning);
  
  std::cout << "Initializing Parameters... \n";
  model.set_paras_initial_values(initial_sigma_sq,
                                 initial_a,
                                 initial_alpha,
                                 initial_delta,
                                 initial_z,
                                 initial_theta,
                                 initial_prob_mat);
  
  std::cout << "Running Gibbs Sampler... \n";
  model.run_gibbs(oc_sigma,
                  oc_a,
                  oc_alpha,
                  oc_delta,
                  oc_theta,
                  oc_z,
                  poly_degree,
                  poly_a,
                  poly_b);
  double elapsed = timer.toc();
  
  List output;
  output = List::create(Named("post_mean") = model.get_gibbs_post_mean(threshold),
                        Named("mcmc") = model.get_gibbs_sample(),
                        Named("trace") = model.get_gibbs_trace(),
                        Named("mcmc_control") = model.get_gibbs_control(),
                        Named("elapsed") = elapsed);
  return output;
}
