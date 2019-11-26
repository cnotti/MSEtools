#define TMB_LIB_INIT R_init_OM

#include <TMB.hpp>
#include "functions.h"

#define DATA_DMATRIX(name)                                               \
matrix<double> name(asMatrix<double>(                                    \
getListElement(TMB_OBJECTIVE_PTR -> data, #name, &Rf_isMatrix)));

#define DATA_DVECTOR(name) vector<double> name(asVector<double>(         \
getListElement(TMB_OBJECTIVE_PTR -> data, #name, &Rf_isNumeric)));

#define DATA_DSCALAR(name)                                               \
double name(asVector<double>(getListElement(TMB_OBJECTIVE_PTR -> data,   \
#name,&isNumericScalar))[0]);

#define DATA_DARRAY(name)                                                \
tmbutils::array<double> name(tmbutils::asArray<double>(                  \
getListElement(TMB_OBJECTIVE_PTR -> data, #name, &Rf_isArray)));

#define DATA_SPARSE_DMATRIX(name)                                        \
Eigen::SparseMatrix<double> name(tmbutils::asSparseMatrix<double>(       \
getListElement(TMB_OBJECTIVE_PTR -> data,                                \
#name, &isValidSparseMatrix)));

// load libs
using namespace Eigen;
  
// define method to extract matrix rows/cols using indexing vectors
template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functor {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const ColIndexType &m_colIndices;
public:
  typedef Matrix<typename ArgType::Scalar,
                 RowIndexType::SizeAtCompileTime,
                 ColIndexType::SizeAtCompileTime,
                 ArgType::Flags&RowMajorBit?RowMajor:ColMajor,
                 RowIndexType::MaxSizeAtCompileTime,
                 ColIndexType::MaxSizeAtCompileTime> MatrixType;
  indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
  {}
  const typename ArgType::Scalar& operator() (Index row, Index col) const {
    return m_arg(m_rowIndices[row], m_colIndices[col]);
  }
};
template <class ArgType, class RowIndexType, class ColIndexType>
CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
mat_indexing(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
  typedef indexing_functor<ArgType,RowIndexType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
}


// data-struct: list of integer vectors 
template<class Type>
struct IVECTORlist_t : vector< vector<int> > {
  IVECTORlist_t(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for (int i=0; i<LENGTH(x); i++) {
      SEXP sm = VECTOR_ELT(x, i);
      if (!Rf_isNumeric(sm)) Rf_error("Not a vector");
      (*this)(i) = asVector<int>(sm);
    }
  }
};


// data-struct: list of integer matrices 
template<class Type>
struct IMATRIXlist_t : vector< matrix<int> > {
  IMATRIXlist_t(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      if(!Rf_isMatrix(sm)) Rf_error("Not a matrix");
      (*this)(i) = asMatrix<int>(sm);
    }
  }
};


// data-struct: list of integer matrices 
template<class Type>
struct DMATRIXlist_t : vector< matrix<double> > {
  DMATRIXlist_t(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      if(!Rf_isMatrix(sm)) Rf_error("Not a matrix");
      (*this)(i) = asMatrix<double>(sm);
    }
  }
};


// cumulative summation function
vector<double> cumsum(vector<double> x) {
  int n = x.size();
 	vector<double> res(n);
 	std::partial_sum(x.data(), x.data() + n, res.data());
 	return res;
}


// median function. Note, approx. for even numbers.
double median(vector<double> v) {
  // calculate the median from all the points
  size_t n = v.size() / 2;
  std::nth_element(v.data(), v.data() + n, v.data() + v.size());
  return v[n];
}


//  Resevoir sample algorithm (Efraimidis and Spirakis 2006)
vector<int> sample(vector<int> x, int size, vector<double> prob) {
  if (x.size() != prob.size()) error("incongruent vector sizes");
  
  // define comparison struct for use within partial_sort
  struct Comp{
    Comp(const vector<double>& v ) : _v(v) {}
    bool operator ()(int a, int b) { return _v[a] < _v[b]; }
    const vector<double>& _v;
  };
  
  // define vars
  int n = x.size();
  vector<double> one_n(n); one_n.setOnes();
  vector<double> rnd = rexp(one_n) / prob;
  vector<int> vx(n);
  // set indices and sample
  for(int i=0; i<n; ++i) vx[i] = i;
  std::partial_sort(vx.data(), vx.data() + size, vx.data() + vx.size(), Comp(rnd));
  return vx.segment(0, size).unaryExpr(x);
}


// remove elements in x that match elments in y
vector<int> drop_y_from_x(vector<int> x, vector<int> y) {
  int nx = x.size();
  int ny = y.size();
  int nNew = nx - ny;
  vector<int> res(nNew);
  const int negInf = -std::numeric_limits<int>::max();
  vector<int> indexes(ny);
  
  for (int i=0; i<ny; i++) {
    auto it = std::find(x.data(), x.data() + nx, y(i));
    indexes(i) = std::distance(x.data(), it);
    x(indexes(i)) = negInf;
  }
  /*std::partial_sort(x.data(), 
                    x.data() + ny, 
                    x.data() + nx);*/
  std::sort(x.data(), x.data() + nx);
  return x.tail(nNew);
}


// function equiv. to the R function which(x == y)
int which_x_equal_y(vector<int> x, int y) {
  auto it = std::find(x.data(), x.data() + x.size(), y);
  return std::distance(x.data(), it);
}


// operating model
template<class Type>
Type objective_function<Type>::operator() () {

  PARAMETER(dummyParameter);   // dummy par required to compile
  
  // ================== ================== ==================
  // Data declaration
  // ================== ================== ==================  
  
  // indexing
  DATA_INTEGER(tstart);             // first year of population dynamics (always 0 in first run)
  DATA_INTEGER(tstop);              // final year of population dynamics (max is nt)
  DATA_INTEGER(nb);                 // number of time-steps in burn in period
  DATA_INTEGER(ny);                 // number of years (should equal nt in estimating model)
  DATA_INTEGER(np);                 // number of inter-annual periods/time-steps
  DATA_INTEGER(nt);                 // total number of time-steps (ny * np)
  DATA_INTEGER(nc);                 // number of species
  DATA_INTEGER(ns);                 // number of vertices in SPDE mesh
  DATA_INTEGER(nx);                 // number of survey tows
  DATA_INTEGER(nl);                 // number of size intervals
  DATA_INTEGER(nf);                 // number of fishable area-swept nodes
  DATA_IVECTOR(p_t);               // period at each time point
  DATA_IMATRIX(f_xq);        	      // association of each survey tow with a given fishing location f
  DATA_IVECTOR(t_catch);            // periods with available(1)/unavailable(0) catch_t data
  DATA_INTEGER(nfishp);             // number of periods in fishing season
  DATA_STRUCT(f_c_fstar, IVECTORlist_t);  // species specific harvest locations f 
  DATA_STRUCT(s_cs_sstar, IVECTORlist_t);  // nodes within spawning range to calculate SSB_s 
  DATA_STRUCT(f_ct_fsurv, IVECTORlist_t); // species specific survey locations f at relevant time points 
  DATA_IMATRIX(ct_ct);
  
  // growth
  DATA_STRUCT(G_c_ll, DMATRIXlist_t);
  
  // natural mortality
  DATA_DARRAY(M_csbl);
  DATA_DARRAY(M_cstl);
  
  // fishing mortality
  DATA_DMATRIX(selectivityF_cl);
  DATA_DVECTOR(limit_c);
  DATA_DSCALAR(areaOmegastar);       // area of spatial domain Omegastar
  DATA_DSCALAR(areadredge);          // area-swept by dredge 
  DATA_SPARSE_DMATRIX(A_fs);         // projection matrix s -> f 
  DATA_DARRAY(ncatch_cstl);
  DATA_DSCALAR(ptarget);             // site sampling probability given B_i > x 
  DATA_DSCALAR(F_intensity);         // proportion of an area that will be fished at time t
  DATA_IVECTOR(F_settings);          // set whether muF is defined by mean or median 
        
  // recruitment
  DATA_DVECTOR(R0_c);
  DATA_DVECTOR(h_c);
  DATA_DMATRIX(pmat_cl);
  DATA_DVECTOR(psi_p);
  DATA_DMATRIX(psi_l);
  DATA_DMATRIX(psi_cs);              // prob of recruitment given environmental conditions at node s
  DATA_DARRAY(E_csb);                // burn-in recruitment innovations
  DATA_DARRAY(E_cst);                // RF on recruitment
  
  // population size structure
  //DATA_DVECTOR(omega_c);
  DATA_DMATRIX(lmid_cl);
  DATA_DMATRIX(weight_cl);
  
  // exit/re-entry R-harvest variables
  DATA_DARRAY(N_csl);                // used for re starting population dynamics at t=tstart
  DATA_DARRAY(SSB0_c);               // used for re starting population dynamics at t=tstart
  DATA_DARRAY(SSB_sc);               // used for re starting population dynamics at t=tstart
  DATA_IVECTOR(ffirst_c);            // grid cell where fishing will start at time t+1
  
  // ================== ================== ==================
  // debugging checks
  // ================== ================== ==================
  
  // debug floating point errors
  //feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW);
  
  //if (nt != ny*np) error("total time periods nt not equal to ny times np");
    
  
  // ================== ================== ==================
  // Model set-up
  // ================== ================== ==================
  
  // additional indicies
  vector<int> nfstar_c(nc);
  for (int c=0; c<nc; c++) {
    nfstar_c(c) = f_c_fstar(c).size();
  }
  
  // variables to initialize and create containers
  vector<int> d_cst(3); d_cst << nc,ns,nt;    // dimension species space time
  vector<int> d_csb(3); d_csb << nc,ns,nb;    // dimension species space time_bip
  vector<int> d_cs(2); d_cs << nc,ns;         // dimension species space
  vector<int> d_ct(2); d_ct << nc,nt;         // dimension species space
  vector< matrix<double> > V_cst(nc*ns*nt);     // vector container to be rearranged into array
  vector< matrix<double> > V_csb(nc*ns*nb);     // vector container to be rearranged into array
  vector< matrix<double> > V_cs(nc*ns);         // vector container to be rearranged into array
  vector< matrix<double> > V_ct(nc*nt);         // vector container to be rearranged into array
  matrix<double> m_l(nl,1); m_l.setOnes();      // matrix to initialize nlx1 structures
  matrix<double> m_1l(1,nl); m_1l.setOnes();    // matrix to initialize 1xnl structures
  vector<int> v_f(nf); v_f.setOnes();         // vector to initialize nf structures
  matrix<double> m_s(ns,1); m_s.setOnes();      // matrix to initialize nsx1 structures
  matrix<double> m_sl(ns,nl); m_sl.setZero();   // matrix to initialize nsxnl structures
  
  //------------------------//
  // initialization section //
  //------------------------//
  
  // numbers at-length
  array< matrix<double> > N_csb_l(V_csb, d_csb);
  array< matrix<double> > N_cst_l(V_cst, d_cst);
  array< matrix<double> > N_ct_sl(V_ct, d_ct);
  array<double> N_cstl(nc,ns,nt,nl);
  
  // recruitment
  matrix<double> R_l(nl,1);
  array< matrix<double> > R_cst_l(V_cst, d_cst);
  array<double> R_ct(nc,nt); R_ct;
  vector< matrix<double> > pmat_c_l(nc);              // proportion mature
  vector<double> maxpsis_c(nc);
  
  // biomass
  array<double> B_ct(nc,nt); B_ct.setZero();
  array<double> B_cst(nc,ns,nt); B_cst.setZero();
  matrix< matrix<double> > B_ct_s(nc,nt);
  matrix< matrix<double> > SSB_cb_s(nc,nb);
  matrix< matrix<double> > SSB_ct_s(nc,nt);
  array<double> SSB_ct(nc,nt); SSB_ct.setZero();      // ssb summed over space
  vector< matrix<double> > weight_c_l(nc);
  vector< matrix<double> > weight_c_1l(nc);
  
  // survival/mortality  
  matrix<double> S_l(nl,1);
  array< matrix<double> > SF_cst_l(V_cst, d_cst);     // survival from fishing
  vector< matrix<double> > selectivityF_c_l(nc);      // fishing selectivity
  vector< matrix<double> > selectivityF_c_1l(nc);     // fishing selectivity
  
  // area of each node
  double areas = areaOmegastar / ns;
  
  // variables used in harvest
  matrix<double> Ad_fs = matrix<double>(A_fs);        // dense proj matrix
  matrix<double> proj_1s(1,ns);
  vector<int> cols_s(ns); std::iota(cols_s.data(), cols_s.data() + ns, 0);
  vector<int> cols_l(nl); std::iota(cols_l.data(), cols_l.data() + nl, 0);
  matrix<double> ncatch_il;
  vector<double> catch_i, prob_f;
  vector<double> B_f(nf);
  matrix<double> N_fl(nf,nl);
  double muBf;
  int niNew, nsamps, nsites, nsites1, nsites2;
  int ni = 0;
  vector<double> catch_n(1); catch_n(0) = 0;
  vector<int> f_n(1);
  vector<int> f_sample, c_i, f_i, t_i;
  matrix<double> scale_cl(nc,nl);
  vector<double> limitp_c = limit_c/nfishp;
  
  vector<int> ftarg;
  vector<int> fsurv;
  int ntarg;
  int nsurv;
  int fone;
  
  
  //--------------------------------------------------//
  // calculate variables to be used in state dynamics //
  //--------------------------------------------------//
  
  // calculate maturity-at-size probabilities and dredge selectivity
  for (int c=0; c<nc; c++) {
    weight_c_l(c) = weight_cl.row(c).transpose();
    weight_c_1l(c) = weight_cl.row(c);
    pmat_c_l(c) = pmat_cl.row(c).transpose();
    selectivityF_c_l(c) = selectivityF_cl.row(c).transpose();
    selectivityF_c_1l(c) = selectivityF_cl.row(c);
    scale_cl.row(c) = selectivityF_c_1l(c) / areas * areadredge;
    maxpsis_c(c) = psi_cs.row(c).maxCoeff();
  }
  
  
  //---------------------------//
  // population initialization //
  //---------------------------//
  
  int p = 0;     // p = annual period
  int cs = 0;
  array<int> n0_c(nc); n0_c.setZero();
  array<double> bswitch_c(nc); bswitch_c.setZero();   // switch used to calculate SSB0_c 
  
  if (tstart == 0) {
    for (int b=0; b<nb; b++) {
      cs = 0;
      for (int c=0; c<nc; c++) {     
        SSB_cb_s(c,b) = m_s;
        for (int s=0; s<ns; s++) {
          // initialize population N_csb_l = 1
          if (b == 0) {
            N_csb_l(c,s,b) = m_l;
          }
          // during first half of initialization calculate SSB0_c
          else if (b < 2*nb/3) {
            // here, recruitment is based off SR = 1 
            R_l = R0_c(c) * areas * E_csb(c,s,b-1) * psi_p(p) * psi_cs(c,s) * psi_l;
            // natural mortality
            for (int l=0; l<nl; l++) {
              S_l(l) = exp(-M_csbl(c,s,b-1,l));
            }
            // state equation
            N_csb_l(c,s,b) = G_c_ll(c) * ((N_csb_l(c,s,b-1) + R_l).cwiseProduct(S_l));
            // SSB estimates
            SSB_cb_s(c,b)(s) = (weight_c_1l(c) * (N_csb_l(c,s,b).cwiseProduct(pmat_c_l(c))))(0,0);
            if (b > nb/3 & psi_cs(c,s) == maxpsis_c(c)) {
              SSB0_c(c) += SSB_cb_s(c,b)(s); // sum over all locations with optimal recruitment conditions
              n0_c(c) += 1;
            }
          }
          // during second half of initialization calculate N_{c,s,t-1,l}
          else {
            // turn off switch indicating first pass of 'else' clause and calc SSB0_c 
            if (bswitch_c(c) == 0) {
              bswitch_c(c) = 1;
              SSB0_c(c) /= n0_c(c);
            }
            calc_recruit(R_l, R0_c(c) * areas, SSB0_c(c), SSB_cb_s(c,b-1), s_cs_sstar(cs),
                         h_c(c), E_csb(c,s,b-1), psi_l, psi_p(p), psi_cs(c,s));      
            // natural mortality
            for (int l=0; l<nl; l++) {
              S_l(l) = exp(-M_csbl(c,s,b-1,l));
            }
            // state equation
            N_csb_l(c,s,b) = G_c_ll(c) * ((N_csb_l(c,s,b-1) + R_l).cwiseProduct(S_l));
            // SSB estimate
            SSB_cb_s(c,b)(s) = (weight_c_1l(c) * (N_csb_l(c,s,b).cwiseProduct(pmat_c_l(c))))(0,0);
          }
          cs += 1;
        }
      }
      p += 1;
      if (p == np) {
        p = 0;
      }
    }
  } else {
    for (int c=0; c<nc; c++) {
      for (int s=0; s<ns; s++) {
        // initialize object
        N_cst_l(c,s,tstart-1) = m_l;
        for (int l=0; l<nl; l++) {
          N_cst_l(c,s,tstart-1)(l) = N_csl(c,s,l);
        }
      }
    }
  }
  
  int debug = 0;
  int DEBUGfone, DEBUGffirst;
  vector<int> DEBUGftarg;

  //----------------------//
  // model state dynamics //
  //----------------------//

  // state dynamics using innovations parameterization
  for (int t=tstart; t<=tstop; t++) {
    cs = 0;
    for (int c=0; c<nc; c++) {
      N_ct_sl(c,t) = m_sl;
      B_ct_s(c,t) = m_s;
      for (int s=0; s<ns; s++) {
        if (t == 0) {
          // numbers at beginning of time t subject to fishing, growth, and natural mortality
          for (int l=0; l<nl; l++) {
            S_l(l) = exp(-M_csbl(c,s,nb-1,l));
          }
          calc_recruit(R_l, R0_c(c) * areas, SSB0_c(c), SSB_cb_s(c,nb-1), s_cs_sstar(cs),
                       h_c(c), E_csb(c,s,nb-1), psi_l, psi_p(p_t(t)), psi_cs(c,s));          
          N_cst_l(c,s,t) = G_c_ll(c) * ( (N_csb_l(c,s,nb-1) + R_l).cwiseProduct(S_l) );
        }
        if (t >= 1) {
          SF_cst_l(c,s,t-1) = m_l;
          SSB_ct_s(c,t-1) = m_s;
          // survival rates
          for (int l=0; l<nl; l++) {
            if (ncatch_cstl(c,s,t-1,l) <= N_cst_l(c,s,t-1)(l)) {
              SF_cst_l(c,s,t-1)(l) = 1 - ncatch_cstl(c,s,t-1,l) / (N_cst_l(c,s,t-1)(l) + 1e-6);
            } else {
              SF_cst_l(c,s,t-1)(l) = 1e-12;
            }
            
            S_l(l) = exp(-M_cstl(c,s,t-1,l));
          }
          // remove catch from t-1
          N_cst_l(c,s,t-1) = N_cst_l(c,s,t-1).cwiseProduct(SF_cst_l(c,s,t-1));
          // biomass calculations
          calc_biomass(B_cst(c,s,t-1), SSB_ct_s(c,t-1)(s), B_ct(c,t-1), SSB_ct(c,t-1), 
                       N_cst_l(c,s,t-1), selectivityF_c_l(c), pmat_c_l(c), weight_c_1l(c));
          // numbers at beginning of time t subject to fishing, growth, and natural mortality
          calc_recruit(R_cst_l(c,s,t-1), R0_c(c) * areas, SSB0_c(c), SSB_ct_s(c,t-1), s_cs_sstar(cs),
                       h_c(c), E_cst(c,s,t-1), psi_l, psi_p(p_t(t)), psi_cs(c,s));
          N_cst_l(c,s,t) = G_c_ll(c) * (N_cst_l(c,s,t-1) + R_cst_l(c,s,t-1)).cwiseProduct(S_l);
        }
        B_ct_s(c,t)(s) = (weight_c_1l(c) * N_cst_l(c,s,t).cwiseProduct(selectivityF_c_l(c)))(0,0) / areas * areadredge;
        // reshape and scale N for projection
        N_ct_sl(c,t).row(s) = (N_cst_l(c,s,t).cwiseProduct(selectivityF_c_l(c))).transpose() / areas * areadredge;
        // iterate cs
        cs += 1;
      }

      
      
      // <><><><><><><><> begin harvest <><><><><><><><>
      if (t_catch(t) == 0) {

        if (f_ct_fsurv(ct_ct(c,t))(0) >= 0) {
          // remove survey sites from commercial domain
          nsurv = f_ct_fsurv(ct_ct(c,t)).size();
          fsurv.resize(nsurv);
          fsurv = f_ct_fsurv(ct_ct(c,t));
          ntarg = nfstar_c(c) - nsurv;
          ftarg.resize(ntarg);
          ftarg = drop_y_from_x(f_c_fstar(c), fsurv);
        } else {
          // no surveys
          nsurv = 0;
          ntarg = nfstar_c(c);
          ftarg.resize(ntarg);
          ftarg = f_c_fstar(c);
        }
        
        // set fone as closest target site to ffirst_c which(targ ~= ffirst_c) 
        (ftarg - ffirst_c(c)).abs().minCoeff(&fone);
        
        if(debug < 1) {
          DEBUGfone = fone;
          DEBUGftarg = ftarg;
          DEBUGffirst = ffirst_c(c);
          debug = 2;
        }
        
        // abundance at each grid cell
        N_fl = A_fs * N_ct_sl(c,t);
        B_f = A_fs * B_ct_s(c,t);

        // calculate site selection probs
        if (F_settings(0) == 0) {
          muBf = B_f.mean();
        } else {
          muBf = median(B_f);
        }
        
        nsamps = std::min(ceil(limitp_c(c) / muBf), ceil(ntarg * 0.9));
        nsites = std::min(ceil(nsamps / F_intensity), ceil(ntarg));
        // if F_intensity > 1 adjust nsamps
        if (nsamps > nsites) {
          nsamps = nsites;
        }
        
        if (nsites < ntarg - fone) {
          // the case where fishing can continue without issue
          f_sample.resize(nsites);
          f_sample = ftarg.segment(fone, nsites);
          fone += nsites;
        } else {
          // the case where fishing must restart part-way through time-step at location 1
          nsites1 = ntarg - fone;
          nsites2 = nsites - nsites1;
          f_sample.resize(nsites);
          f_sample.tail(nsites1) = ftarg.segment(fone, nsites1);
          f_sample.head(nsites2) = ftarg.segment(0, nsites2);
          fone = nsites2;
        }
        
        prob_f.resize(nsites);
        for (int f=0; f < nsites; f++) {
          if (B_f(f_sample(f)) > muBf) {
            prob_f(f) = ptarget;
          } else {
            prob_f(f) = 1 - ptarget;
          }
        }
        
        while (catch_n.sum() < limitp_c(c) & nsamps <= nsites) {
          // resize containers
          f_n.resize(nsamps);
          catch_n.resize(nsamps);
          // sample catch
          f_n = sample(f_sample, nsamps, prob_f);
          catch_n = f_n.unaryExpr(B_f);
          nsamps *= 1.2; // + nsites - nsamps
        }

        // select nsize ref closest to catch limit
        (cumsum(catch_n) - limitp_c(c)).abs().minCoeff(&nsamps);
        ++nsamps;
        
        // save vector of catch and book-keeping vars
        niNew = ni + nsamps + nsurv;
        // resize containers
        ncatch_il.conservativeResize(niNew,nl);
        catch_i.conservativeResize(niNew);
        f_i.conservativeResize(niNew);
        c_i.conservativeResize(niNew);
        t_i.conservativeResize(niNew);
        // fill re-sized areas with sampled data
        if (nsurv > 0) {
          catch_i.segment(ni, nsurv) = fsurv.unaryExpr(B_f);
          ncatch_il.block(ni,0,nsurv,nl) = mat_indexing(N_fl, fsurv, cols_l);
          f_i.segment(ni, nsurv) = fsurv;
        }
        ncatch_il.bottomRows(nsamps) = mat_indexing(N_fl, f_n.head(nsamps), cols_l);
        catch_i.tail(nsamps) = catch_n.head(nsamps);
        f_i.tail(nsamps) = f_n.head(nsamps);
        std::fill(c_i.data() + ni, c_i.data() + niNew, c);
        std::fill(t_i.data() + ni, t_i.data() + niNew, t);

        // map catch back to triangulation nodes
        proj_1s = mat_indexing(Ad_fs, f_i.tail(nsurv + nsamps), cols_s).colwise().sum();
        for (int s=0; s<ns; ++s) {
          for (int l=0; l<nl; ++l) {
            ncatch_cstl(c,s,t,l) = N_ct_sl(c,t)(s,l) * proj_1s(s);
          }
        }

        // update variables for harvest in next time period
        ffirst_c(c) = ftarg(fone);
        ni = niNew;
        catch_n.setZero();
      }
      // <><><><><><><><> end harvest <><><><><><><><>
    }
  }
  
  // Reported variables
  for (int c=0; c<nc; c++) { 
    for (int s=0; s<ns; s++) {
      for (int l=0; l<nl; l++) {
        N_csl(c,s,l) = N_cst_l(c,s,tstop)(l);
      }
    }
  }
  
  REPORT(DEBUGfone);
  REPORT(DEBUGffirst);
  REPORT(DEBUGftarg);

  // harvest data
  REPORT(SSB0_c);
  REPORT(N_csl);
  REPORT(B_cst);
  REPORT(ncatch_cstl);
  REPORT(c_i);
  REPORT(f_i);
  REPORT(t_i);
  REPORT(catch_i);

  return 0;
}
