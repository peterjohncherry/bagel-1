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
  string IBra_name = gamma_info->Bra_info->name();

  cout << "got " << gamma_name << " info "<< endl;

  //note this has reverse iterators!
  if (gamma_info->order > 2 ) { 

    get_gamma(gamma_info->predecessor_gamma_name);
    compute_sigmaN( gamma_info->predecessor_gamma_name, gamma_info->name );

  } else {

    string JKet_name = gamma_info->Ket_info->name();
    get_wfn_data( gamma_info->Bra_info );
    get_wfn_data( gamma_info->Ket_info );

    compute_sigma2_test( IBra_name, JKet_name );
    compute_sigma2( IBra_name, JKet_name, gamma_info->name );
    get_gamma2_from_sigma2_and_civec( IBra_name, gamma_name );

 }


  cout << "D" << endl;

  return gamma_vec;                              
  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  Equation_Computer::Equation_Computer::compute_sigma2_test( string Bra_name, string Ket_name ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "Equation_Computer::compute_sigma2_test" <<endl;

  shared_ptr<Determinants> Ket_det = det_old_map->at( Ket_name ); 
  shared_ptr<Civec> Ket_A = make_shared<Civec>(* (cvec_old_map->at( Ket_name )));
  shared_ptr<Civec> Ket_B = make_shared<Civec>(* (cvec_old_map->at( Ket_name )));
  int norb = Ket_det->norb();
  
  shared_ptr<Dvec> sigma2_A = make_shared<Dvec>( Ket_det, (norb*norb) );
  shared_ptr<Dvec> sigma2_B = make_shared<Dvec>( Ket_det, (norb*norb) );


  double* Ket_B_ptr  = Ket_B->data();
  double* sigma2_B_ptr  = sigma2_B->data(0)->data();

  sigma_2a1(Ket_A, sigma2_A, Ket_det );
  sigma_2a2(Ket_A, sigma2_A, Ket_det );

  sigma_2a1( Ket_B_ptr, sigma2_B_ptr, Ket_det );
  sigma_2a2( Ket_B_ptr, sigma2_B_ptr, Ket_det );

  for (int ii = 0 ; ii != norb*norb*Ket_det->size() ; ii++ ) 
    if ( sigma2_A->data(0)->data(ii) - sigma2_B->data(0)->data(ii) != 0.0) 
      cout << "2a2_diff[" << ii <<  "] = " << sigma2_A->data(0)->data(ii) - sigma2_B->data(0)->data(ii) << endl;

  cout << "gamma2_test orig routines" << endl;
  for ( int ii = 0 ; ii != norb; ii++ ){ 
    for ( int jj = 0 ; jj != norb; jj++ ){
       cout << ddot_( Ket_det->size(), sigma2_A->data(norb*ii + jj)->data(), 1, Ket_A->data(), 1) << " "; cout.flush(); 
     } 
    cout <<endl;
  }

  cout << "gamma2_test ptr_routines" << endl;
  for ( int ii = 0 ; ii != norb; ii++ ){ 
    for ( int jj = 0 ; jj != norb; jj++ ){
       cout << ddot_( Ket_det->size(), sigma2_B->data(norb*ii + jj)->data(), 1, Ket_B->data(), 1) << " "; cout.flush(); 
     } 
    cout <<endl;
  }
  cout <<endl; cout << " ham " <<endl; 
  return;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::compute_sigmaN( string predecessor_gamma_name, string sigmaN_name)  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::compute_sigmaN" << endl;

  shared_ptr<GammaInfo>  predecessor_info = GammaMap->at(predecessor_gamma_name);
  shared_ptr<GammaInfo>  sigmaN_info = GammaMap->at(sigmaN_name);

  string  Bra_name = sigmaN_info->Bra_info->name();
  string  Ket_name = sigmaN_info->Ket_info->name();
  

  if (Bra_name == Ket_name) { 
    cout << "(Bra_name == Ket_name)  =  " <<  Bra_name << " == " << Ket_name << endl;
    shared_ptr<Civec>  IBra = cvec_old_map->at( Bra_name );          cout << "get dets" << endl;

    int sorder = sigmaN_info->order; 
    shared_ptr<Determinants> Ket_det = det_old_map->at( Ket_name ); 
    shared_ptr<Dvec>      pred_sigma = dvec_sigma_map->at( predecessor_gamma_name );

    int orb_dim = pow(Ket_det->norb(), sorder-2);
    int orb2    = Ket_det->norb()*Ket_det->norb();

    shared_ptr<Dvec> sigmaN = make_shared<Dvec>( Ket_det, orb_dim*orb2  );
    cout << "orb_dim = " << orb_dim << endl;
    for ( int  ii = 0; ii != orb_dim; ii++) {
      sigma_2a1( pred_sigma->data(ii)->data(), sigmaN->data(ii*orb2)->data(), Ket_det );
      sigma_2a2( pred_sigma->data(ii)->data(), sigmaN->data(ii*orb2)->data(), Ket_det );
    }

    dvec_sigma_map->emplace( sigmaN_name, sigmaN );
    unique_ptr<double[]> gammaN( new double[orb_dim*orb2]);
    std::fill_n(gammaN.get(), orb_dim*orb2, 0.0);

    for ( int  ii = 0; ii != orb_dim*orb2; ii++) 
      gammaN[ii] = ddot_( Ket_det->size(), sigmaN->data(ii)->data(), 1, IBra->data(), 1); 
     
    int norb = Ket_det->norb();
    cout << "--------------- gamma " << sorder << " ------------------- (from recursive)" << endl;
    for ( int ii = 0; ii != orb_dim*orb2 ; ii++ ) {
      if ( !( ii % norb )  ) {
        cout << endl;
        if (! (ii % orb2) ) {
          cout <<  " [ " ; cout.flush();
          int ii_tmp = ii;
          for (int jj = sorder-1 ; jj != 1 ; jj-- ) {   
            int id_pos = ii_tmp / ((int)pow( norb , jj ));
            cout << id_pos  << " " ;  
            ii_tmp -= id_pos*(int)pow(norb, jj );
          }
          cout << "] " << endl; 
        }
      }
      cout << gammaN[ii] << " " ; cout.flush();
    }

    if (sorder == 4 ) { 
     int n2 = norb*norb;
     int n3 = n2*norb;

     unique_ptr<double[]> gamma2( new double[n2]);
     std::fill_n(gamma2.get(), n2, 0.0);
     for ( int  ii = 0; ii != n2; ii++) 
       gamma2[ii] = ddot_( Ket_det->size(), pred_sigma->data(ii)->data(), 1, IBra->data(), 1); 
 
     cout << "------------ tmp gamma2 -------------" << endl;
     for ( int ii = 0 ; ii != norb; ii++ ) {
       for ( int jj = 0 ; jj != norb; jj++ ) {
          cout << gamma2[ii*norb+jj] << " " ; cout.flush();
       } 
       cout <<endl;
     }
     cout <<endl;


     for ( int ii = 0 ; ii != norb; ii++ ) {
       for ( int jj = 0 ; jj != norb; jj++ ) {
         for ( int kk = 0 ; kk != norb; kk++ ) {
           if ( jj == kk ) {
             for ( int ll = 0 ; ll != norb; ll++ ) 
               gammaN[ ii*n3+ jj*n2 + kk*norb + ll ] -= gamma2[ii*norb+ll] ;
           }
         }
       }
     }

    cout << "--------------- -rdm4 " << sorder << " ------------------- (from recursive)" << endl;
    for ( int ii = 0; ii != orb_dim*orb2 ; ii++ ) {
      if ( !( ii % norb )  ) {
        cout << endl;
        if (! (ii % orb2) ) {
          cout << endl <<  " [ " ; cout.flush();
          int ii_tmp = ii;
          for (int jj = sorder-1 ; jj != 1 ; jj-- ) {   
            int id_pos = ii_tmp / ((int)pow( norb , jj ));
            cout << id_pos  << " " ;  
            ii_tmp -= id_pos*(int)pow(norb, jj );
          }
          cout << "] " << endl; 
        }
      }
      cout << gammaN[ii] << " " ; cout.flush();
    }

   }


  } else {
  
    cout << "spin transitions sigmas not implemented yet " << endl;  
  
  }
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::sigma_2a1(double* cvec_ptr, double* sigma_ptr, shared_ptr<Determinants> dets  )  {
//////////////////////////////////////////////////////////////////////////////////////////////                               
  //cout << "sigma_2a1" << endl;                                                                                          
  const int lb = dets->lenb();                                                                                            
  const int ij = dets->norb()*dets->norb();                                                                               
                                                                                                                          
  for (int ip = 0; ip != ij; ++ip) {                                                                                      
    double* target_base = sigma_ptr+dets->size()*ip;                                                                              
    for (auto& iter : dets->phia(ip)) {                                                                                   
      const double sign = static_cast<double>(iter.sign);                                                                 
      double* const target_array = target_base + iter.source*lb;                                                          
      blas::ax_plus_y_n(sign, cvec_ptr + iter.target*lb, lb, target_array);                                               
    }                                                                                                                     
  }                                                                                                                       
}                                                                                                                         
                                                                                                                          
///////////////////////////////////////////////////////////////////////////////////////////////                             
void Equation_Computer::Equation_Computer::sigma_2a2( double* cvec_ptr, double* sigma_ptr, shared_ptr<Determinants> dets) { 
///////////////////////////////////////////////////////////////////////////////////////////////                             
//  cout << "sigma_2a2" << endl;
  const int la = dets->lena();
  const int lb = dets->lenb();
  const int ij = dets->norb()*dets->norb();

  for (int i = 0; i < la; ++i) {
    double* source_array0 = cvec_ptr+i*lb;

    for (int ip = 0; ip != ij; ++ip) {
      double* target_array0 = sigma_ptr + (ip * dets->size()) + i*lb;

      for (auto& iter : dets->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::get_gamma2_from_sigma2_and_civec( string IBra_name, string sigma2_Ket_name){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "get_gamma2_from_sigma2_and_civec" << endl; 

  cout << "sigma2_Ket_name = " << sigma2_Ket_name << endl;
  shared_ptr<Dvec>         sigma2_Ket = dvec_sigma_map->at( sigma2_Ket_name );  cout << "IBra_name = " << IBra_name << endl;
  shared_ptr<Civec>        IBra       = cvec_old_map->at( IBra_name );          cout << "get dets" << endl;
  shared_ptr<Determinants> dets       = det_old_map->at( IBra_name ); 
 
  cout << "ham " << endl; 
  int norb = dets->norb();

  unique_ptr<double[]> gamma2_data(new double[norb*norb]);

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
void Equation_Computer::Equation_Computer::compute_sigma2( string IBra_name, string JKet_name, string sigma2_JKet_name)  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "compute_sigma2" << endl;

  if (IBra_name == JKet_name) { 

    shared_ptr<Determinants> JKet_det = det_old_map->at( JKet_name ); 
    shared_ptr<Civec>        JKet = cvec_old_map->at( JKet_name );
    
    shared_ptr<Dvec> sigma2_JKet = make_shared<Dvec>( JKet_det, JKet_det->norb()*JKet_det->norb() );
    
    sigma_2a1( JKet, sigma2_JKet, JKet_det );
    sigma_2a2( JKet, sigma2_JKet, JKet_det );
    
    dvec_sigma_map->emplace( sigma2_JKet_name, sigma2_JKet );

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
