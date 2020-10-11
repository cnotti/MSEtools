#define TMB_LIB_INIT R_init_SIM
#include <TMB.hpp>
#include "functions.h"

// model
template<class Type>
Type objective_function<Type>::operator() () {

  // load libs
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;

  PARAMETER(dummyParameter);   // dummy par required to compile
  
  // ================== ================== ==================
  // Data declaration
  // ================== ================== ==================  
  
  // options
  DATA_STRING(g_I);                 // log | logit
  DATA_STRING(g_H);                 // log | logit
  DATA_STRING(g_E);                 // log | logit
 
  // indexing
  DATA_INTEGER(nb);                 // number of time-steps in burn in period
  DATA_INTEGER(nt);                 // total number of time-steps (ny * np)
  DATA_INTEGER(nc);                 // number of species
  DATA_INTEGER(ns);                 // number of vertices in SPDE mesh
  DATA_INTEGER(nl);                 // number of size intervals
  DATA_INTEGER(np);                 // number of periods
  DATA_IVECTOR(splus_s);            // identfier for splus nodes within depth domain

  // natural mortality
  DATA_VECTOR(alpha_c); 
  DATA_SCALAR(mugI);                // mean of gI
  DATA_SCALAR(mugH);                // mean of gH
  DATA_SCALAR(sigmagI);             // sd of gI
  DATA_SCALAR(taugH);               // scales marginal SD of gH (g[Eta])
  DATA_SCALAR(kappaH);              // range of spatial variation (Eta)
  DATA_SCALAR(phiIt);               // temporal corelation (pr Eta)
  DATA_SCALAR(phiHt);               // temporal corelation (Eta)
  DATA_SCALAR(phiHl);               // corelation between sizes (Eta)
  
  // recruitment
  DATA_SCALAR(mugE);                // mean of gE
  DATA_SCALAR(taugE);               // scales marginal SD of gE (g[Epsilon])
  DATA_SCALAR(kappaE);              // range of spatial variation (Epsilon)
  DATA_SCALAR(phiEt);               // temporal corelation (Epsilon)
  DATA_MATRIX(Rho_cc);              // species recruitment correlation matrix
  
  // SPDE objects
  DATA_STRUCT(spde, spde_aniso_t);  // spde objects
  DATA_VECTOR(lnaniso);             // anisotrpy
  

  // ================== ================== ==================
  // Simulation procedure for fields
  // ================== ================== ==================
  
  // anisotropy
  matrix<Type> aniso(2,2);
  aniso(0,0) = exp(lnaniso(0));
  aniso(1,0) = lnaniso(1);
  aniso(0,1) = lnaniso(1);
  aniso(1,1) = (1+lnaniso(1)*lnaniso(1))/exp(lnaniso(0));

  // SPDE spatial precision matrices
  SparseMatrix<Type> QH = Q_spde(spde, kappaH, aniso);
  SparseMatrix<Type> QE = Q_spde(spde, kappaE, aniso);
  
  
  // number of spatial nodes + additional nodes on boundary
  int nsplus = QH.rows();
  
  // simulate random effects as innovations and complete program execution
  SIMULATE{
    
    MVNORM_t<Type> MVN_c(Rho_cc);
    
    array<Type> gE_csn(nc,nsplus,nb+nt);
    array<Type> E_csb(nc,ns,nb);
    array<Type> E_cst(nc,ns,nt);
    Type Emean = 0;
    
    array<Type> gI_n(nb+nt);
    array<Type> I_b(nb);
    array<Type> I_t(nt);
    Type Imean = 0;
    
    array<Type> gH_snl(nsplus,nb+nt,nl);
    array<Type> H_sbl(ns,nb,nl);
    array<Type> H_stl(ns,nt,nl);
    Type Hmean = 0;
		
    array<Type> M_csbl(nc,ns,nb,nl);
    array<Type> M_cstl(nc,ns,nt,nl);
		Type Mmean = 0;
    
    SEPARABLE(AR1(phiEt),
              SEPARABLE(GMRF(QE), MVN_c)).simulate(gE_csn);
    
    AR1(phiIt).simulate(gI_n);
    
    SEPARABLE(AR1(phiHl), 
              SEPARABLE(AR1(phiHt), GMRF(QH))).simulate(gH_snl);

    for (int c=0; c<nc; c++) {
      for (int s=0; s<ns; s++) {  
        for (int t=0; t<(nb+nt); t++) {
          if (t < nb) {
            if (c == 0) {
              if (s == 0) {
                I_b(t) = invlogit(gI_n(t) * sigmagI + mugI);
                Imean += I_b(t);
              }
              for (int l=0; l<nl; l++) {
                H_sbl(s,t,l) = exp(gH_snl(splus_s(s),t,l) / taugH + mugH);
                Hmean += H_sbl(s,t,l);
								Mmean += H_sbl(s,t,l) * I_b(t);
              }
            }
            E_csb(c,s,t) = invlogit(gE_csn(c,splus_s(s),t) / taugE + mugE);
            Emean += E_csb(c,s,t);
          } else {
            if (c == 0) {
              if (s == 0) {
                I_t(t-nb) = invlogit(gI_n(t) * sigmagI + mugI);
                Imean += I_t(t-nb);
              }
              for (int l=0; l<nl; l++) {
                H_stl(s,t-nb,l) = exp(gH_snl(splus_s(s),t,l) / taugH + mugH);
                Hmean += H_stl(s,t-nb,l);
								Mmean += H_stl(s,t-nb,l) * I_t(t-nb);
              }
            }
            E_cst(c,s,t-nb) = invlogit(gE_csn(c,splus_s(s),t) / taugE + mugE);
            Emean += E_cst(c,s,t-nb);
          }
					
        }
      }
    }
    
    Emean /= (nc*ns*(nb+nt));
    Imean /= (nb+nt);
    Hmean /= (ns*(nb+nt)*nl);
		Mmean /= (ns*(nb+nt)*nl);
    
    for (int c=0; c<nc; c++) {
      for (int s=0; s<ns; s++) {  
        for (int l=0; l<nl; l++) {
          for (int b=0; b<nb; b++) {        
            M_csbl(c,s,b,l) = alpha_c(c)/np * I_b(b)*H_sbl(s,b,l)/Mmean; //(Imean*Hmean);
          }
          for (int t=0; t<nt; t++) {
            M_cstl(c,s,t,l) = alpha_c(c)/np * I_t(t)*H_stl(s,t,l)/Mmean; //(Imean*Hmean);
          }
        }
      }
    }
    
    REPORT(gE_csn);
    REPORT(E_csb);
    REPORT(E_cst);
    REPORT(Emean);
    
    REPORT(I_b);
    REPORT(I_t);
    REPORT(Imean);
    
    REPORT(gH_snl);
    REPORT(H_sbl);
    REPORT(H_stl);
    REPORT(Hmean);
    
    REPORT(M_csbl);
    REPORT(M_cstl);
    
  }
  // finish program execution
  return 0;
}
