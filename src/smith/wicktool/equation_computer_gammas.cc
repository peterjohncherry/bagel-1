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
Equation_Computer::Equation_Computer::get_gamma(string gamma_name){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::Equation_Computer::get_gamma"  << endl;
  shared_ptr<Tensor_<double>> New_gamma = make_shared<Tensor_<double>>(); 
  shared_ptr<vector<shared_ptr<Tensor_<double>>>> gamma_vec = make_shared<vector<shared_ptr<Tensor_<double>>>>(0);
 
  cout << "gamma_name = " <<  gamma_name << endl;

  shared_ptr<GammaInfo> gamma_info = GammaMap->at(gamma_name);

  cout << "got " << gamma_name << " info "<< endl;

  //note this has reverse iterators!
  if (gamma_info->id_ranges->size() > 2 ) { 
    
    get_gamma(gamma_info->predecessor_gamma_name);
 
  } else {

    string IBra_name = gamma_info->Bra_info->name();
    string JKet_name = gamma_info->Ket_info->name();

    get_wfn_data( gamma_info->Bra_info );
    get_wfn_data( gamma_info->Ket_info );
   
    compute_sigma2( IBra_name, JKet_name );
    get_gamma2_from_sigma2_and_civec( IBra_name, JKet_name );
 }

  cout << "D" << endl;

  return gamma_vec;                              
  
}
   
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::get_wfn_data( shared_ptr<CIVecInfo<double>>  cvec_info ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string cvec_name = cvec_info->name();
  auto cvec_det_loc = det_old_map->find(cvec_name); 
  if( cvec_det_loc == det_old_map->end()){
    int II = cvec_info->state_num();
    shared_ptr<const Determinants> det_cvec_orig = cc_->data(II)->det();
    shared_ptr<Determinants> det_cvec = make_shared<Determinants>(det_cvec_orig->norb(), det_cvec_orig->nelea(), det_cvec_orig->neleb(), false, /*mute=*/true);
    det_old_map->emplace( cvec_name, det_cvec);
    cvec_old_map->emplace( cvec_name, make_shared<Civec>(*(cc_->data(II))) );
  }  

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::get_gamma2_from_sigma2_and_civec( string IBra_name, string sigma2_Ket_name){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "get_gamma2_from_sigma2_and_civec" << endl; 

  shared_ptr<Dvec>         sigma2_Ket = dvec_sigma_map->at(sigma2_Ket_name); 
  shared_ptr<Civec>        IBra       = cvec_old_map->at(IBra_name);
  shared_ptr<Determinants> dets       = det_old_map->at(sigma2_Ket_name); 
  
  int norb = dets->norb();

  unique_ptr<double[]> gamma2_data(new double[ norb*norb]);

  for ( int ii = 0 ; ii != norb ; ii++) 
    for ( int jj = 0 ; jj != norb ; jj++) 
      *(gamma2_data.get() + (ii*norb+jj)) =  ddot_( sigma2_Ket->det()->size(),  sigma2_Ket->data(ii*norb+jj)->data(), 1, IBra->data(), 1); 

  cout << "gamma2_data " << endl;
  for ( int ii = 0 ; ii != norb ; ii++) {
    
    for ( int jj = 0 ; jj != norb ; jj++) { 
      cout << *(gamma2_data.get()+ii*norb+jj) << " " ;cout.flush();
    }
    cout <<endl;
  }

  return; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::compute_sigma2( string IBra_name, string JKet_name )  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "compute_sigma2" << endl;

  if (IBra_name == JKet_name) { 

    shared_ptr<Determinants> JKet_det = det_old_map->at( JKet_name ); 
    shared_ptr<Civec>        JKet = cvec_old_map->at( JKet_name );
    
    shared_ptr<Dvec> sigma2_JKet = make_shared<Dvec>( JKet_det, JKet_det->norb()*JKet_det->norb() );
    
    sigma_2a1( JKet, sigma2_JKet, JKet_det );
    sigma_2a2( JKet, sigma2_JKet, JKet_det );
    
    dvec_sigma_map->emplace( JKet_name, sigma2_JKet );
    dvec_sigma_map->emplace( IBra_name, sigma2_JKet );

  } else {

    cout << "spin transitions sigmas not implemented yet " << endl;  

  }
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Taken directly from src/ci/fci/knowles_compute.cc         
//////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::sigma_2a1(shared_ptr<const Civec> cvec, shared_ptr<Dvec> sigma, shared_ptr<Determinants> dets  ) const {
//////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "sigma_2a1" << endl;
  const int lb = dets->lenb();
  const int ij = dets->norb()*dets->norb();
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
  const int la = dets->lena();
  const int ij = dets->norb()*dets->norb();

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

///////////////////////////////////////////////////////////////////////////////////////////////
//   shared_ptr<Dvec> sigma_kl_JKet;
//  if (IBra != JKet) {
//
//    sigma_kl_JKet = make_shared<Dvec>(det_JKet, norb_*norb_);
//
//    sigma_2a1(JKet, sigma_kl_JKet, det_JKet);
//    sigma_2a2(JKet, sigma_kl_JKet, det_JKet);
//    dvec_sigma_map->emplace( JJ_name, sigma_kl_JKet ) ;
//
//  } else {
//   sigma_ij_IBra = sigma_kl_JKet;
//  }


#endif
