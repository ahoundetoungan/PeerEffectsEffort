// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;
using namespace std;


// Function for foptim
//[[Rcpp::export]]
List fsigma2eps(const double& tau,
                const double& rho,
                const arma::vec& resids,
                List& W,
                List& WW,
                List& I,
                const arma::vec& nst,
                const int& sumn,
                const int& S){
  int n1(0), n2(-1);
  double sig2eps(0);
  List lOm(S);
  for(int s(0); s < S; ++ s){
    n1            = n2 + 1;
    n2            = n1 + nst(s) - 1;
    arma::vec es  = resids.subvec(n1, n2);
    arma::mat Ws  = W[s];
    arma::mat WWs = WW[s];
    arma::mat Is  = I[s];
    arma::mat Om  = Is + tau*tau*WWs + rho*tau*(Ws + arma::trans(Ws));
    lOm[s]        = Om;
    sig2eps      += arma::cdot(es, arma::solve(Om, es));
  }
  
  return List::create(Named("se2") = sig2eps/sumn, Named("Om") = lOm);
}

//[[Rcpp::export]]
double floglike(const arma::vec& x,
                const arma::vec& resids,
                List& W,
                List& WW,
                List& I,
                const arma::vec& nst,
                const int& sumn,
                const int& S){
  double tau      = exp(x(0));
  double rho      = (exp(x(1)) - 1)/(exp(x(1)) + 1);
  //Rprintf("tau: %f\n", tau);
  //Rprintf("rho: %f\n", rho);
  
  double sig2eps, llh(0), val, sign;
  List lOm;
  
  {List tmp1      = fsigma2eps(tau, rho, resids, W, WW, I, nst, sumn, S);
    double tmp2   = tmp1[0];
    List tmp3     = tmp1[1];
    sig2eps       = tmp2;
    lOm           = tmp3;}
  
  
  //Rprintf("sigma^2_epsilon: %f\n", sig2eps);
  for(int s(0); s < S; ++ s){
    arma::mat Om  = lOm[s];
    log_det(val, sign, Om); 
    llh          += (-0.5*nst(s)*(log(sig2eps) + 1) - 0.5*(val + log(sign)));
  }
  //Rprintf("log(likelihood): %f\n", llh);
  //Rprintf("*****************\n");
  return -llh;
}

//[[Rcpp::export]]
List fdataFs(List& J){
  int S           = J.length();
  arma::vec eigval;
  arma::uvec ind;
  arma::mat eigvec;
  List F(S);
  for(int s(0); s < S; ++ s){
    //Rprintf("Find Eigenvalues and Eigenvectors: %i/%i\n", s + 1, S);
    arma::mat Js = J[s];
    arma::eig_sym(eigval, eigvec, Js, "std");
    ind          = arma::find(eigval > 0.999);
    eigvec       = eigvec.cols(ind);
    F[s]         = eigvec.t();
  }
  return F;
}


//[[Rcpp::export]]
List ftoolsml(const arma::vec& resids, 
              List& F,
              List& network,
              const double& lambda,
              const int& S){
  arma::vec es, n(S);
  arma::uvec ind;
  arma::mat Is, Ws, WWs, As;
  int n1(0), n2(-1), ns, nss;
  List I(S), W(S), WW(S), lresid(S);
  for(int s(0); s < S; ++ s){
    arma::mat Gs = network[s];
    ns           = Gs.n_rows;
    As           = arma::eye(ns, ns) - lambda*Gs;
    arma::mat Fs = F[s];
    nss          = Fs.n_rows;
    Is           = arma::eye(nss, nss);
    Ws           = Fs*As;
    WWs          = Ws*Ws.t();
    Ws           = Ws*Fs.t();
    n1           = n2 + 1;
    n2           = n1 + Fs.n_rows - 1;
    //cout<<n2<<endl;
    es           = resids.subvec(n1, n2);
    n(s)         = nss;
    I[s]         = Is;
    W[s]         = Ws;
    WW[s]        = WWs;
    lresid[s]    = es;
  }
  return List::create(Named("I") = I, Named("W") = W, Named("WW") = WW, Named("lresid") = lresid, Named("n") = n);
}



//////////////////////////////////////// Tobit model
// log of determinant of SPD matrices
double logdetSPD(const Eigen::MatrixXd& St) {
  Eigen::LLT<Eigen::MatrixXd> lltS(St);
  if (lltS.info() != Eigen::Success) {
    throw std::runtime_error("LLT failed: matrix may be semi-definite or ill-conditioned.");
  }
  const Eigen::MatrixXd& L = lltS.matrixL();
  return 2.0 * L.diagonal().array().log().sum();
}

// log of determinant of SPD matrices from their L form
double logdetLform(const Eigen::MatrixXd& L) {
  return 2.0 * L.diagonal().array().log().sum();
}

// Simulate from normal
Eigen::VectorXd fnormal(const int& k,
                        std::mt19937& rng) {
  Eigen::VectorXd theta(k);
  std::normal_distribution<> dist(0.0, 1.0);
  for(int i = 0; i < k; ++i) {
    theta(i) = dist(rng);
  }
  return theta;
}

// Simulate from inverse W
Eigen::MatrixXd rinvwishart(const int& k,
                            const int& df, 
                            const Eigen::MatrixXd& S, 
                            std::mt19937& rng) {
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(k, k);
  
  // Fill A using Bartlett decomposition
  std::normal_distribution<> normal(0.0, 1.0);
  
  for(int i = 0; i < k; ++i){
    std::chi_squared_distribution<> chi = std::chi_squared_distribution<>(df - i);
    A(i,i) = std::sqrt(chi(rng));
    for(int j=0; j<i; ++j){
      A(i,j) = normal(rng);
    }
  }
  
  Eigen::MatrixXd L(Eigen::LLT<Eigen::MatrixXd>(S).matrixL());
  return L * (A * A.transpose()).llt().solve(L.transpose());
}


/* Simulation from truncated normal [low, +inf]
 * following Robert, C. P. (1995). Simulation of truncated normal variables.
 Statistics and computing, 5(2), 121-125.*/
double tnorm(const double& low, 
             std::mt19937& rng) {
  double z(0);
  bool repeat(true);
  
  if(low <= 0){
    std::normal_distribution<> Norm(0.0, 1.0);
    while (repeat) {
      z = Norm(rng);
      repeat = (z <= low);
    }
  } else {
    double alpha(0.5 * (low + sqrt(pow(low, 2.0) + 4.0))), e(0), rho(0), u(0);
    std::uniform_real_distribution<> Unif(0.0, 1.0);
    std::exponential_distribution<> Exp(1.0);
    while (repeat) {
      e = Exp(rng);
      z = low + e / alpha ;
      
      rho = exp(-pow(alpha - z, 2.0) / 2) ;
      u   = Unif(rng);
      if (u <= rho) {
        repeat = false;
      }
    }
  }
  return z;
}

// update theta = [beta, gamma]
void updtheta(Eigen::VectorXd& parm,
              Eigen::VectorXd& Vtheta,
              const Eigen::VectorXd& Az,
              const std::vector<Eigen::LLT<Eigen::MatrixXd>>& lltS,
              const Eigen::MatrixXd& V,
              const int Ktheta,
              const int& M,
              const Eigen::ArrayXi& nvec,
              const Eigen::ArrayXi& igroup,
              const Eigen::VectorXd& m0theta,
              const Eigen::MatrixXd& iS0theta,
              std::mt19937& rng,
              const int& nthreads) {
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif
  std::vector<Eigen::MatrixXd> lViSV(nthreads);
  std::vector<Eigen::VectorXd> lViSzt(nthreads);
#pragma omp parallel 
{
  int tid     = omp_get_thread_num();
  lViSV[tid]  = Eigen::MatrixXd(Ktheta, Ktheta);
  lViSzt[tid] = Eigen::VectorXd(Ktheta);
  lViSV[tid].setZero();
  lViSzt[tid].setZero();
#pragma omp for schedule(static)
  for (int m = 0; m < M; ++m) {
    int nm(nvec(m));
    Eigen::MatrixXd Vm(V.block(igroup(m), 0, nm, Ktheta));
    lViSV[tid]  += (Vm.transpose() * lltS[m].solve(Vm));
    lViSzt[tid] += (Vm.transpose() * lltS[m].solve(Az.segment(igroup(m), nm)));
  }
}

Eigen::MatrixXd ViSV(lViSV[0]);
Eigen::VectorXd ViSzt(lViSzt[0]);
for (int k = 1; k < nthreads; ++k) {
  ViSV  += lViSV[k];
  ViSzt += lViSzt[k];
}

// posterior mean and variance
Eigen::MatrixXd iShat(ViSV + iS0theta);
Eigen::LLT<Eigen::MatrixXd> lltiS(iShat);
if (lltiS.info() != Eigen::Success) {
  throw std::runtime_error("LLT failed: matrix may be semi-definite or ill-conditioned.");
}
Eigen::VectorXd mhat(lltiS.solve(ViSzt + iS0theta * m0theta));

// simulation
parm.segment(1, Ktheta) = lltiS.matrixU().solve(fnormal(Ktheta, rng)) + mhat;
Vtheta = V * parm.segment(1, Ktheta);
}

// update seta, sepsilon, rho
void updSigma(Eigen::VectorXd& parm,
              std::vector<Eigen::MatrixXd>& S,
              std::vector<Eigen::LLT<Eigen::MatrixXd>>& lltS,
              const Eigen::VectorXd& Az,
              const std::vector<Eigen::MatrixXd>& A,
              const std::vector<Eigen::MatrixXd>& AAT,
              const Eigen::VectorXd& Vtheta,
              const int Ktheta,
              const int& M,
              const Eigen::ArrayXi& nvec,
              const Eigen::ArrayXi& igroup,
              const double& df0,
              const Eigen::MatrixXd scale0, 
              std::mt19937& rng,
              const int& nthreads) {
  // Extract seeds for each thread
  std::vector<std::uint32_t> seeds(nthreads);
  std::uniform_int_distribution<std::uint32_t> Unif;
  for (int t = 0; t < nthreads; ++t) {
    seeds[t] = Unif(rng);
  }
  
  // residual
  Eigen::VectorXd u(Az - Vtheta);
  
  // (eta, epsilon)
  Eigen::MatrixXd ee(nvec.sum(), 2);
  
  // variance components
  double seta(parm(Ktheta + 1)), sepsilon(parm(Ktheta + 2)), rho(parm(Ktheta + 3));
  
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif
#pragma omp parallel 
{
  int tid     = omp_get_thread_num();
  std::mt19937 rngt(seeds[tid]);
#pragma omp for schedule(static)
  for (int m = 0; m < M; ++m) {
    int nm(nvec(m));
    // residual
    Eigen::VectorXd um(u.segment(igroup(m), nm));
    
    // S11
    Eigen::MatrixXd S11(Eigen::MatrixXd::Identity(nm, nm));
    S11.diagonal() *= (seta * seta);
    
    // S21
    Eigen::MatrixXd S21(Eigen::MatrixXd::Identity(nm, nm));
    S21.diagonal() *= rho * seta * sepsilon;
    S21 +=  seta * seta * A[m]; 
    
    // posterior mean and variance of eta
    Eigen::VectorXd mud(S21.transpose() * lltS[m].solve(um));
    Eigen::MatrixXd Sd(S11 - S21.transpose() *  lltS[m].solve(S21));
    
    // simulation of eta
    Eigen::LLT<Eigen::MatrixXd> lltSd1(Sd);
    if (lltSd1.info() != Eigen::Success) {
      throw std::runtime_error("LLT failed: matrix may be semi-definite or ill-conditioned.");
    }
    ee.block(igroup(m), 0, nm, 1) = lltSd1.matrixL() * fnormal(nm, rngt) + mud;
    
    // computation of epsilon
    ee.block(igroup(m), 1, nm, 1) = um - A[m] * ee.block(igroup(m), 0, nm, 1);
  }
}

// covarianve matrix
Eigen::MatrixXd Sigma = rinvwishart(2, df0 + nvec.sum(), scale0 + ee.transpose() * ee, rng);

// update seta, sepsilon, and rho
seta     = std::sqrt(Sigma(0, 0));
sepsilon = std::sqrt(Sigma(1, 1));
rho      = Sigma(0, 1) / (seta * sepsilon);
parm(Ktheta + 1) = seta;
parm(Ktheta + 2) = sepsilon;
parm(Ktheta + 3) = rho;

// update S
#pragma omp parallel for schedule(static)
for (int m = 0; m < M; ++m) {
  S[m] = seta * seta * AAT[m] + rho * seta * sepsilon * (A[m] + A[m].transpose());
  S[m].diagonal().array() += sepsilon * sepsilon; 
  lltS[m] = S[m].llt();
  if (lltS[m].info() != Eigen::Success) {
    throw std::runtime_error("LLT failed: S[m] is not positive definite.");
  }
}
}


// update zeta, where lambda = (exp(zeta) - 1)/(exp(zeta) + 1)
void updzeta(Eigen::VectorXd& parm,
             std::vector<Eigen::MatrixXd>& S,
             std::vector<Eigen::LLT<Eigen::MatrixXd>>& lltS,
             Eigen::VectorXd& Az,
             std::vector<Eigen::MatrixXd>& A,
             std::vector<Eigen::MatrixXd>& AAT,
             std::vector<Eigen::MatrixXd>& LAAT,
             unsigned int& zetaaccept,
             const Eigen::VectorXd& z,
             const std::vector<Eigen::MatrixXd>& G,
             const Eigen::VectorXd& Vtheta,
             const int Ktheta,
             const int& M,
             const Eigen::ArrayXi& nvec,
             const Eigen::ArrayXi& igroup,
             const double& jumpzeta,
             const double& m0zeta,
             const double& iS0zeta, 
             std::mt19937& rng,
             const int& nthreads) {
  double zeta(std::log(parm(0) + 1) - log(1 - parm(0)));
  std::normal_distribution<> Norm(zeta, jumpzeta);
  std::uniform_real_distribution<> Unif(0, 1);
  double zetast(Norm(rng));
  double lambdast((exp(zetast) - 1.0) / (exp(zetast) + 1.0));
  double seta(parm(Ktheta + 1)), sepsilon(parm(Ktheta + 2)), rho(parm(Ktheta + 3));
  Eigen::ArrayXd rate(M);
  
  // Potential objects to update
  std::vector<Eigen::MatrixXd> St(M);
  std::vector<Eigen::LLT<Eigen::MatrixXd>> lltSt(M);
  Eigen::VectorXd Azt(nvec.sum());
  std::vector<Eigen::MatrixXd> At(M);
  std::vector<Eigen::MatrixXd> AATt(M);
  std::vector<Eigen::MatrixXd> LAATt(M);
  
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for schedule(static)
  for (int m = 0; m < M; ++m) {
    int nm(nvec(m));
    At[m] = -lambdast * G[m];
    At[m].diagonal().array() += 1;
    Azt.segment(igroup(m), nm) = At[m] * z.segment(igroup(m), nm);
    Eigen::VectorXd res(Az.segment(igroup(m), nm) - Vtheta.segment(igroup(m), nm));
    Eigen::VectorXd rest(Azt.segment(igroup(m), nm) - Vtheta.segment(igroup(m), nm));
    AATt[m] = At[m] * At[m].transpose();
    St[m] = seta * seta * AATt[m] + rho * seta * sepsilon * (At[m] + At[m].transpose());
    St[m].diagonal().array() += sepsilon * sepsilon; 
    lltSt[m] = St[m].llt();
    if (lltSt[m].info() != Eigen::Success) {
      throw std::runtime_error("LLT failed: S[m] is not positive definite.");
    }
    Eigen::LLT<Eigen::MatrixXd> lltAATt(AATt[m]);
    if (lltAATt.info() != Eigen::Success) {
      throw std::runtime_error("LLT failed: S[m] is not positive definite.");
    }
    LAATt[m] = lltAATt.matrixL();
    rate(m) = 0.5 * (logdetLform(lltS[m].matrixL()) - logdetLform(lltSt[m].matrixL())) + 
      0.5 * (logdetLform(LAATt[m]) - logdetLform(LAAT[m])) + 
      0.5 * (res.dot(lltS[m].solve(res)) - rest.dot(lltSt[m].solve(rest))) + 
      0.5 * (std::pow(zeta - m0zeta, 2) / iS0zeta - std::pow(zetast - m0zeta, 2) / iS0zeta);
  }
  
  if(std::log(Unif(rng)) < rate.sum()){
    parm(0) = lambdast;
    zetaaccept +=1;     //Increase acceptance number to 1 
    S = St;
    lltS = lltSt;
    Az = Azt;
    A = At;
    AAT = AATt;
    LAAT = LAATt;
  }
}


// update z: data augmentation
void updz(Eigen::VectorXd& Az,
          Eigen::VectorXd& z,
          const Eigen::VectorXd& parm,
          const Eigen::VectorXd& y,
          const double miny,
          const double maxy,
          const double& lambda,
          const std::vector<Eigen::MatrixXd>& S,
          const std::vector<Eigen::LLT<Eigen::MatrixXd>>& lltS,
          const std::vector<Eigen::MatrixXd>& A,
          const Eigen::VectorXd& Vtheta,
          const int& Ktheta,
          const int& M,
          const Eigen::ArrayXi& nvec,
          const Eigen::ArrayXi& igroup,
          std::mt19937& rng,
          const int& nthreads) {
  // Extract seeds for each thread
  std::vector<std::uint32_t> seeds(nthreads);
  std::uniform_int_distribution<std::uint32_t> Unif;
  for (int t = 0; t < nthreads; ++t) {
    seeds[t] = Unif(rng);
  }
  
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif
#pragma omp parallel 
{
  int tid     = omp_get_thread_num();
  std::mt19937 rngt(seeds[tid]);
#pragma omp for schedule(static)
  for (int m = 0; m < M; ++m) {
    int nm(nvec(m));
    // mean of z
    Eigen::PartialPivLU<Eigen::MatrixXd> luA(A[m]);
    Eigen::VectorXd mz(luA.solve(Vtheta.segment(igroup(m), nvec(m))));
    // Inverse of variance of z
    Eigen::MatrixXd iSz(A[m].transpose() * lltS[m].solve(A[m]));
    // simulate z
    for (int i = 0; i < nm; ++i) {
      // conditional standard dev
      double szi(std::sqrt(1/iSz(i, i)));
      // conditional mean
      double mzi(z(igroup(m) + i) - iSz.row(i).transpose().dot(z.segment(igroup(m), nm) - mz) / iSz(i, i));
      if (y(igroup(m) + i) <= miny) {
        z(igroup(m) + i) = -tnorm((mzi - miny) / szi, rngt) * szi + mzi;
      }
      else if (y(igroup(m) + i) >= maxy) {
        z(igroup(m) + i) = tnorm((maxy - mzi) / szi, rngt) * szi + mzi;
      }
    }
    Az.segment(igroup(m), nm) = (A[m] * z.segment(igroup(m), nm));
  }
}
}

//[[Rcpp::export]]
Rcpp::List fTobit(const Eigen::VectorXd& y,
                  const Eigen::MatrixXd& V,
                  const std::vector<Eigen::MatrixXd>& G,
                  const int sim = 1000,
                  const int& nthreads = 1,
                  const Rcpp::Nullable<double>& lby = R_NilValue,
                  const Rcpp::Nullable<double>& uby = R_NilValue,
                  const double& target  = 0.4,
                  const double& jumpmin = 1e-3,
                  const double& jumpmax = 1,
                  const Rcpp::Nullable<Eigen::VectorXd>& parm0 = R_NilValue,
                  const Rcpp::Nullable<Eigen::VectorXd>& z0 = R_NilValue,
                  const Rcpp::Nullable<unsigned long long> seed = R_NilValue) {
  // parameter;
  int Ktheta(V.cols());
  Eigen::VectorXd parm(Ktheta + 4);
  if (parm0.isNotNull()) {
    parm = Rcpp::as<Eigen::VectorXd>(parm0);
    if (parm.size() != (Ktheta + 4)) {
      Rcpp::stop("The initial value `parm0` does not suit the model: The size is not compatible.");
    }
    if (std::abs(parm(0)) >= 1) {
      Rcpp::stop("The initial value `parm0` does not suit the model: peer effects are large.");
    }
    if ((parm(Ktheta + 1) <= 0) || (parm(Ktheta + 2) <= 0)) {
      Rcpp::stop("The initial value `parm0` does not suit the model: the standard errors are negative or zero.");
    }
    if (std::abs(parm(Ktheta + 3)) >= 1) {
      Rcpp::stop("The initial value `parm0` does not suit the model: The correlation is larger than one.");
    }
  } else {
    parm << 0.1, Eigen::VectorXd::Zero(Ktheta), 1, 1, 0;
  }
  
  // prior
  double m0zeta(0.5), iS0zeta(1);
  Eigen::VectorXd m0theta(Eigen::VectorXd::Zero(Ktheta));
  Eigen::MatrixXd iS0theta(Eigen::MatrixXd::Identity(Ktheta, Ktheta) / 100);
  double df0(4); //p + 2
  Eigen::MatrixXd scale0 = Eigen::MatrixXd::Identity(2, 2);
  
  // min and max y
  double miny(y.minCoeff());
  double maxy(y.maxCoeff());
  if (lby.isNotNull()) {
    miny = Rcpp::as<double>(lby);
  }
  if (uby.isNotNull()) {
    maxy = Rcpp::as<double>(uby);
  }
  if (miny >= maxy) {
    Rcpp::stop("The lower bound of y is greater than or equal to the upper bound of y.");
  }
  
  // seed
  unsigned long long base_seed;
  if (seed.isNotNull()) {
    base_seed = Rcpp::as<unsigned long long>(seed);
  } else {
    std::random_device rd;
    base_seed = rd();
  }
  std::mt19937 rng(base_seed);
  
  // group
  int M(G.size());
  Eigen::ArrayXi nvec(M), igroup(M + 1);
  igroup(0) = 0;
  for (int m = 0; m < M; ++m) {
    nvec(m)       = G[m].rows();
    igroup(m + 1) = igroup(m) + nvec(m);
  }
  int N(nvec.sum());
  
  // zaccept and jumping scale
  unsigned int zetaaccept(0);
  double jumpzeta(0.1);
  
  // Objects to update
  Eigen::VectorXd z(y);
  if (z0.isNotNull()) {
    z = Rcpp::as<Eigen::VectorXd>(z0);
  }
  z = ((y.array() > miny) && (y.array() < maxy)).select(y, z);
  if (z.size() != N) {
    Rcpp::stop("`z0` is not an N-vector.");
  }
  Eigen::VectorXd Vtheta(V * parm.segment(1, Ktheta));
  Eigen::VectorXd Az(N);
  std::vector<Eigen::MatrixXd> A(M), S(M), AAT(M), LAAT(M);
  std::vector<Eigen::LLT<Eigen::MatrixXd>> lltS(M);
  double lambda(parm(0)), seta(parm(Ktheta + 1)), sepsilon(parm(Ktheta + 2)), rho(parm(Ktheta + 3));
  
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for schedule(static)
  for (int m = 0; m < M; ++m) {
    int nm(nvec(m));
    A[m]    = Eigen::MatrixXd::Identity(nm, nm) - lambda * G[m];
    Az.segment(igroup(m), nm) = A[m] * z.segment(igroup(m), nm);
    AAT[m]  = A[m] * A[m].transpose();
    S[m]    = seta * seta * AAT[m] + rho * seta * sepsilon * (A[m] + A[m].transpose());
    S[m].diagonal().array() += sepsilon * sepsilon; 
    lltS[m] = S[m].llt();
    Eigen::LLT<Eigen::MatrixXd> lltAAT(AAT[m]);
    if (lltAAT.info() != Eigen::Success) {
      throw std::runtime_error("LLT failed: S[m] is not positive definite.");
    }
    LAAT[m] = lltAAT.matrixL();
  }
  
  // storing
  Eigen::MatrixXd Posterior(Ktheta + 4, sim);
  Eigen::MatrixXd zsave(N, std::min(1000, sim));
  
  // acceptance rate 
  double arate;
  
  // MCMC
  Rcpp::Rcout << "Iteration 0/" << sim << std::endl;
  Rcpp::Rcout << "    [lambda, seta, sepsilon, rho] = ["
              << parm(0) << ", " << parm(Ktheta + 1) << ", "
              << parm(Ktheta + 2) << ", " << parm(Ktheta + 3)
              << "]" << std::endl;
  Rcpp::Rcout << std::endl;
  for (int t = 0; t < sim; ++t) {
    // update theta = [beta, gamma]
    updtheta(parm, Vtheta, Az, lltS, V, Ktheta, M, nvec, igroup, m0theta,
             iS0theta, rng, nthreads);
    
    // update lambda = (exp(zeta) - 1)/(exp(zeta) + 1)
    updzeta(parm, S, lltS, Az, A, AAT, LAAT, zetaaccept, z, G, Vtheta, Ktheta,
            M, nvec, igroup, jumpzeta, m0zeta, iS0zeta, rng, nthreads);
    
    // update jumping scale
    arate     = (zetaaccept + 0.0)/(t + 1.0);
    jumpzeta += (arate - target) / pow(t + 1.0, 0.6);
    jumpzeta  = std::clamp(jumpzeta, jumpmin, jumpmax);
    
    // update seta, sepsilon, rho
    updSigma(parm, S, lltS, Az, A, AAT, Vtheta, Ktheta, M, nvec, igroup, df0,
             scale0, rng, nthreads);
    
    // update z: data augmentation
    updz(Az, z, parm, y, miny, maxy, lambda, S, lltS, A, Vtheta, Ktheta, M, nvec, igroup, rng,
         nthreads);
    
    Posterior.col(t)    = parm;
    zsave.col(t % 1000) = z;
    // print progression
    Rcpp::Rcout << "Iteration " << t + 1<< "/" << sim << std::endl;
    Rcpp::Rcout << "    [lambda, seta, sepsilon, rho] = ["
                << parm(0) << ", " << parm(Ktheta + 1) << ", "
                << parm(Ktheta + 2) << ", " << parm(Ktheta + 3)
                << "]" << std::endl;
    Rcpp::Rcout << "    Acceptance rate = "<< arate << std::endl;
    Rcpp::Rcout << std::endl;
  }
  return Rcpp::List::create(Rcpp::_["parms"] = Posterior.transpose(), 
                            Rcpp::_["z"]     = zsave);
}