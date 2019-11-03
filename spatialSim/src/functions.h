// rmultinom
array<double> rmultinom(double n, matrix<double> p) {
  int nk = p.size();
  array<double> samples(nk);
  double pused = 0;
  // sample from conditional binomals to obtain multinomial sample
  for (int k=0; k<nk-1; k++) {
    samples(k) = rbinom(n, p(k)/(1 - pused));
    n -= samples(k);
    pused += p(k);
  }
  samples(nk-1) = n;
  return samples;
}

// dlnorm
double dlnorm(double x, double mean, double sd, bool loged = true) {
  double logres = dnorm(log(x), mean, sd, true) - log(x);
  if (loged) {
    return logres;
  } else {
    return exp(logres);
  }
}

// define lognormal_uniform CDF
double plnu(double x, double mu, double sigma, double omega) {
  double pr;
  if((x - omega) <= 0) {
    pr = (x*pnorm(log(x), mu, sigma) -
      exp(mu + 0.5*pow(sigma,2)) *
      pnorm((log(x) - mu - pow(sigma,2))/sigma, double(0), double(1)))/omega;
  }
  else{
    pr = (x*pnorm(log(x), mu, sigma) -
      exp(mu + 0.5*pow(sigma,2)) *
      pnorm((log(x) - mu - pow(sigma,2))/sigma, double(0), double(1)) -

      (x-omega)*pnorm(log(x-omega), mu, sigma) +
      exp(mu + 0.5*pow(sigma,2)) *
      pnorm((log(x-omega) - mu - pow(sigma,2))/sigma, double(0), double(1)))/omega;
  }
  return pr;
}
  
// define lognormal_uniform CDF between a and b
double plnu_ab(double mu, double sigma, double a, double b, double omega, 
             double (*plnu)(double, double, double, double)) {
  double pr;
  if (a <= 0) {
    pr = (*plnu)(b, mu, sigma, omega);
  }
  else{
    pr = (*plnu)(b, mu, sigma, omega) - (*plnu)(a, mu, sigma, omega);
  }
  return pr;
}

// logistic curve
double logistic(double l, double l50, double l95) {
  return 1/(1+exp(-log(19)*((l-l50)/(l95-l50))));
}

// decreasing logistic growth curve
double dec_logistic(double l1, double delta_max, double l50, double l95) {
  return delta_max/(1+exp(log(19)*((l1-l50)/l95)));
}

double vonBert(double l, double linf, double beta) {
  return (1 - exp(-beta)) * (linf - l);
}

// function for detecting NAs
bool isNA(double x) {
  return R_IsNA(asDouble(x));
}

// function to insure variable remains positive and differentiable
double posfun(double x, double eps, double &pen) {
  pen += CppAD::CondExpLt(x, eps, double(0.01) * pow(x-eps,2), double(0));
  return CppAD::CondExpGe(x, eps, x, eps/(double(2)-x/eps));
}

// inverse of scaled and translated logit function
double invlogitST(double x, double a, double b) {
  return a + (b - a)*invlogit(x);
}

array<double> invlogit_array(array<double> x){
  array<double> a = x;
  a = exp(x) / (double(1) + exp(x));
  return a;
}

// backtransform random field
array<double> invRF(array<double> rf, double mu, double tau, std::string trsfrm) {
  array<double> ans = rf;
  if (trsfrm == "log") {
    ans = exp(rf/tau + mu);
  }
  else if (trsfrm == "logit") {
    ans = exp(rf/tau + mu) / (double(1.0) + exp(rf/tau + mu));
  } else {
    error("invalid RF transformation");
  }
  return ans;
}  

double invRF(double rf, double mu, double tau, std::string trsfrm) {
  if (trsfrm == "log") {
    return exp(rf/tau + mu);
  }
  else if (trsfrm == "logit") {
    return exp(rf/tau + mu) / (1 + exp(rf/tau + mu));
  } else {
    error("invalid RF transformation");
  }
}  

// calculate growth transition matrix
void calc_G(matrix<double>& G, vector<double> lmid, double omega,
                    double delta_max, double l50_tag, 
                    double l95_tag, double sigma_tag,
                    double (*g)(double, double, double, double),
                    double (*plnu)(double, double, double, double),
                    int nl, int np) {
  vector<double> delta(nl);
  // calculate G (Gy = annual, Gp = inter-annual)
  for (int r=0; r<nl; r++) {
    // expected annual growth increment given initial size r
    delta(r) = g(lmid(r), delta_max, l50_tag, l95_tag);
    // prob growing to each size class c
    for (int c=0; c<nl; c++) {
      if (r <= c) {
        if((c-r)*omega <= 0) {
          G(r,c) = plnu(omega+(c-r)*omega,
                        log(delta(r)/np) - pow(sigma_tag,2)/2,
                        sigma_tag,
                        omega);
        } else{
          G(r,c) = plnu(omega+(c-r)*omega,
                        log(delta(r)/np) - pow(sigma_tag,2)/2,
                        sigma_tag,
                        omega) - 
                   plnu((c-r)*omega,
                        log(delta(r)/np) - pow(sigma_tag,2)/2,
                        sigma_tag,
                        omega);
        }
      }
      if (r > c) G(r,c) = double(0);
    }
    G.row(r) = G.row(r)/G.row(r).sum();
  }
}



// recruitment calculations
void calc_recruit(matrix<double>& R_cst_l, 
                  double R0,
                  double SSB0_c,
                  matrix<double> SSB_ct_s,
                  vector<int> s_cs_sstar,                  
                  double h, 
                  double E_cst, 
                  matrix<double> psi_l, 
                  double psi_p,
                  double psi_d) {
  int nsstar = s_cs_sstar.size();
  double SSBavail = 0;
  for (int sstar = 0; sstar < nsstar; sstar++) {
    SSBavail += SSB_ct_s(s_cs_sstar(sstar));
  }
  double pSSB = SSBavail/nsstar / SSB0_c;
  R_cst_l = R0 * pSSB / (1 - ((5 * h - 1)/(4 * h))*(1 - pSSB)) * E_cst * psi_p*psi_d*psi_l;
}


          
// biomass calculations
void calc_biomass(double& B_cst, double& SSB_cst, 
                  double& B_ct, double& SSB_ct, 
                  matrix<double> N_l, matrix<double> selectivity_l, 
                  matrix<double> pmat_l, matrix<double> weight_1l) {
  B_cst = (weight_1l * (N_l.cwiseProduct(selectivity_l)))(0,0);
  SSB_cst = (weight_1l * (N_l.cwiseProduct(pmat_l)))(0,0);
  B_ct += B_cst;
  SSB_ct += SSB_cst;
}


// exploitation and survival rates
void calc_survival(double& SFpos_cst, 
                   double& ERate_cst, 
                   matrix<double>& SF_cst_l, 
                   double x, 
                   double eps, 
                   double &pen,
                   matrix<double> selectivity, 
                   matrix<double> v_l,
                   double (*posfun)(double, double, double&)) {
  SFpos_cst = (*posfun)(x, eps, pen);
  ERate_cst = 1 - SFpos_cst;
  SF_cst_l = v_l - (selectivity*ERate_cst);
}

// convert class AD --> int (has the same effect as floor)
matrix<int> asInt(matrix<double> x) {
  int nc = x.cols();
  int nr = x.rows();
  matrix<int> x_int(nr,nc);
  for (int r=0; r<nr; r++) {
    for (int c=0; c<nc; c++) {
      x_int(r,c) = CppAD::Integer(x(r,c));
    } 
  }
  return x_int;
}


//// template<class double>
//double std_dev() = ((s1.rowwise() - s1.colwise().mean()).square().colwise().sum()/(M-1)).sqrt();



// calculate haversine great circle distance between lat long points

double distance(matrix<double> lonlat1, matrix<double> lonlat2) {
  double lon1 = lonlat1(0,0) * M_PI / 180.0;
  double lat1 = lonlat1(0,1) * M_PI / 180.0;
  double lon2 = lonlat2(0,0) * M_PI / 180.0;
  double lat2 = lonlat2(0,1) * M_PI / 180.0;
  double d_lat = abs(lat1 - lat2);
  double d_lon = abs(lon1 - lon2);
  double a = pow(sin(d_lat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(d_lon / 2), 2);
  double d_sigma = 2 * asin(sqrt(a));
  return d_sigma * 6371000;
}



// A function to generate a random permutation of v 
void randomize(vector<int>& v) { 
    int n = v.size(); 
    // srand( time(NULL) ); // randomize seed
    for (int i = n-1; i > 0; i--) { 
        // Pick a random index from 0 to i 
        int j = rand() % (i+1); 
        // Swap v(i) with the element at random index j
        int temp = v(i); 
        v(i) = v(j); 
        v(j) = temp; 
    }   
}
