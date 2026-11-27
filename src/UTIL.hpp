
#ifndef UTIL_HPP
#define UTIL_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <sys/time.h>
#include <fstream>

const static std::unordered_map<std::string,int> string_to_case{
   {"best_guess",1},
   {"mean",2},
   {"minor",3}
};

double getWeights(std::string t_kernel, 
                  double t_freq, 
                  arma::vec t_wBeta);

void imputeGeno(arma::vec& GVec, 
                double freq, 
                std::vector<uint32_t> posMissingGeno);

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat);


bool imputeGenoAndFlip(arma::vec& t_GVec,
                       double & t_altFreq,
		       double & t_altCount,
                       std::vector<uint32_t> &  t_indexForMissing,
                       std::string t_impute_method,
                       double t_dosage_zerod_cutoff,
                       double t_dosage_zerod_MAC_cutoff,
                       double & t_MAC,
		       std::vector<uint> & t_indexZero,
                       std::vector<uint> & t_indexNonZero);

arma::vec getTime();

void printTime(arma::vec t1, arma::vec t2, std::string message);

double getinvStd(double t_freq);

// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
// void set_seed(unsigned int seed) {
//   Rcpp::Environment base_env("package:base");
//   Rcpp::Function set_seed_r = base_env["set.seed"];
//   set_seed_r(seed);  
// };

// 
// arma::vec nb(int n){
//   return(Rcpp::rbinom(n,1,0.5));
// }

arma::vec nb(unsigned int n);

double sum_arma1(arma::vec& X);

double add_logp(double p1, double p2);

double calculatePearsonCorrelation(arma::vec& x, arma::vec& y);

// Flushes, closes and clears the streams; returns the first of t_names whose write, flush or close failed, else "".
std::string close_outfiles(std::vector<std::ofstream> & t_streams, std::vector<std::string> & t_names);

// Removes a temp file. A failure is reported on stderr, not raised; a file already gone is no failure.
void remove_temp_file(const std::string & t_file);

// Removes the recorded temp files on scope exit, so also when the call throws.
struct TempFileGuard {
  std::vector<std::string> files;
  ~TempFileGuard(){ for(const std::string & file : files) remove_temp_file(file); }
};

// Saves a matrix to a temp file recorded in t_guard; a failed save raises an R error naming the file.
template <typename T>
void save_temp_matrix(const T & t_mat, const std::string & t_file, TempFileGuard & t_guard){
  t_guard.files.push_back(t_file);
  if(!t_mat.save(t_file)) Rcpp::stop("Could not write the temp file " + t_file);
}
#endif

