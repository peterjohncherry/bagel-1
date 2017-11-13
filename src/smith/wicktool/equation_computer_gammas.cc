#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/equation_computer.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets the gammas in tensor format. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<Tensor_<double>>>>
Equation_Computer::Equation_Computer::get_gammas(int II, int JJ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::Equation_Computer::get_gammas"  << endl;
  shared_ptr<Tensor_<double>> New_gamma = make_shared<Tensor_<double>>(); 
  shared_ptr<vector<shared_ptr<Tensor_<double>>>> gamma_vec = make_shared<vector<shared_ptr<Tensor_<double>>>>(0);
    
  dvec_sigma_map = make_shared<std::map< std::string, std::shared_ptr<Dvec>>>();

  cout << "A " << endl;
  shared_ptr<const Determinants> det_IBra_orig = cc_->data(II)->det();
  shared_ptr<const Determinants> det_JKet_orig = cc_->data(JJ)->det(); 

  cout << "B " << endl;

  shared_ptr<Determinants> det_IBra = make_shared<Determinants>(det_IBra_orig->norb(), det_IBra_orig->nelea(), det_IBra_orig->neleb(), false, /*mute=*/true);
  shared_ptr<Determinants> det_JKet = make_shared<Determinants>(det_JKet_orig->norb(), det_JKet_orig->nelea(), det_JKet_orig->neleb(), false, /*mute=*/true);

  cout << "C " << endl;

  compute_gamma12( cc_->data(II), cc_->data(JJ), det_IBra, det_JKet, to_string(II), to_string(JJ));

  return gamma_vec;                              
  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Equation_Computer::Equation_Computer::compute_gamma12( shared_ptr<const Civec> IBra, shared_ptr<const Civec> JKet,
                                                       shared_ptr<Determinants> det_IBra, shared_ptr<Determinants> det_JKet,
                                                       string II_name, string JJ_name ) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "compute_gamma12_from_civec" << endl;

  shared_ptr<Dvec> sigma_ij_IBra = make_shared<Dvec>(det_IBra, det_IBra->norb()*det_IBra->norb());

  cout << "D1" <<endl;
  sigma_2a1(IBra, sigma_ij_IBra, det_IBra);
  cout << "D1" <<endl;
  sigma_2a2(IBra, sigma_ij_IBra, det_IBra);
  cout << "D3" <<endl;

  cout << " sigma_II_name = " << II_name << endl;  
  cout << " sigma_JJ_name = " << JJ_name << endl;  

  dvec_sigma_map->emplace( II_name, sigma_ij_IBra ) ;
  cout << "put into map" << endl;
 
  shared_ptr<Dvec> sigma_kl_JKet;
  if (IBra != JKet) {

  cout << "D4" <<endl;
    sigma_kl_JKet = make_shared<Dvec>(det_JKet, norb_*norb_);
  cout << "D5" <<endl;

    sigma_2a1(JKet, sigma_kl_JKet, det_JKet);
  cout << "D6" <<endl;
    sigma_2a2(JKet, sigma_kl_JKet, det_JKet);
  cout << "D7" <<endl;

    dvec_sigma_map->emplace( JJ_name, sigma_kl_JKet ) ;
  cout << "D8" <<endl;

  } else {
    sigma_ij_IBra = sigma_kl_JKet;
  }

  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Taken directly from src/ci/fci/knowles_compute.cc         
//////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::sigma_2a1(shared_ptr<const Civec> cvec, shared_ptr<Dvec> sigma, shared_ptr<Determinants> dets  ) const {
//////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "sigma_2a1" << endl;
  const int lb = sigma->lenb();
  const int ij = sigma->ij();
  const double* const source_base = cvec->data();

  for (int ip = 0; ip != ij; ++ip) {
    double* const target_base = sigma->data(ip)->data();

    //DetMap(size_t t, int si, size_t s, int o) : target(t), sign(si), source(s), ij(o) {}
    for (auto& iter : dets->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      double* const target_array = target_base + iter.source*lb;
      blas::ax_plus_y_n(sign, source_base + iter.target*lb, lb, target_array);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Taken directly from src/ci/fci/knowles_compute.cc         
// I'm sure there was a version which used transposition of the civector; this looks slow.
///////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::sigma_2a2(shared_ptr<const Civec> cvec, shared_ptr<Dvec> sigma, shared_ptr<Determinants> dets) const {
///////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "sigma_2a2" << endl;
  const int la = sigma->lena();
  const int ij = sigma->ij();

  for (int i = 0; i < la; ++i) {
    const double* const source_array0 = cvec->element_ptr(0, i);

    for (int ip = 0; ip != ij; ++ip) {
      double* const target_array0 = sigma->data(ip)->element_ptr(0, i);

      for (auto& iter : dets->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}

 
#endif
