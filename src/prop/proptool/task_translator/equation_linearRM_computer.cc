
#include <src/prop/proptool/task_translator/equation_computer.h>
#include <src/prop/proptool/task_translator/equation_linearRM_computer.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_LinearRM<DataType>::solve_equation(){  
//////////////////////////////////////////////////////////////////////////////////////////////////

  // \sum{NN} < L | X_{MM, NN} (H - e_{QQ} -e_{s}) \X_{MM,NN} | N >
  for (int QQ = 0; QQ != ref_space_dim; ++QQ)  
    for ( int MM = 0; MM != ref_space_dim; ++MM )  
      shared_ptr<MultiTensor_<double>> residues = this->get_tensop_vector( "residual", make_pair("T1", MM ), make_pair("e_shift", QQ) ) ;

  // solve linear equation for t amplitudes
  bool converged = true;
  for (int LL = 0; LL != ref_space_dim; LL++) { 
    bool conv = false;
    shared_ptr<MultiTensor_<DataType>> residues = this->get_tensop_vector( "residual", make_pair("bra", LL) ) ;
    shared_ptr<MultiTensor_<DataType>> pt_amps  = this->get_tensop_vector( "pt_amps", make_pair("bra", LL) );

    shared_ptr<Tensor_<DataType>> pt_term = this->get_tensop( "pt_term" , make_pair("state", LL ) );
    if ( pt_term->rms() < 1.0e-15) {
      if (LL+1 != ref_space_dim) cout << endl;
      continue;
    } else {
//    Both of these have a sum over NN 
//    update_amplitude(pt_amps, pt_term);
    }
    
//    auto solver = make_shared<LinearRM<MultiTensor>>( max_davidson_subspace_size_, pt_term[LL]);
    int max_iter = 5;
    for (int iter = 0; iter != max_iter; ++iter) {
      residues->zero();
    
      //TODO this step is present in the original formulation but I do not quite see why it must be normalized yet...
      DataType norm = pt_amps->norm();
      pt_amps->scale(1.0/norm);
    
      // Calculate residuals
      for (int MM = 0; MM != ref_space_dim; ++MM)  // MM bra vector
        for (int NN = 0; NN != ref_space_dim; ++NN)  // NN ket vector
          this->evaluate_expression( "residual", make_pair("T1",MM), make_pair("T2", NN), make_pair("bra", LL)  );
        
      // these are Bagel routines, nothe that the compute residual here is _not_ the residual form the ci equations
      // residues = solver->compute_residual(pt_amps, residues);
      // pt_amps = solver->civec();
    
      // energy is now the Hylleraas energy
      //energy_[LL] = detail::real(dot_product_transpose(proj_op[LL], pt_amps));
      //energy_[LL] += detail::real(dot_product_transpose(residues[LL], pt_amps));
    
      // compute rms for state i
      // error = residues->norm() / pow(residues->size(), 0.25); // Why is this power to a 1/4?
      // conv = error < conv_threshold;
    
      // compute delta t2 and update amplitude
      if (!conv) {
        pt_amps->zero();
        //update_amplitude(pt_amps, residues);
      }
      if (conv) break;
    }
    if (LL+1 != ref_space_dim) cout << endl;
    converged &= conv;
  }

  //This is the last term in the effective Hamiltonian    
  for (int NN = 0; NN != ref_space_dim; ++NN) {
    DataType norm = this->get_scalar_result( "norm", make_pair("state",NN ));
    norm = (DataType)(0.0);
    for (int LL = 0; LL != ref_space_dim; ++LL)  // bra
      for (int MM = 0; MM != ref_space_dim; ++MM)  // ket
        this->evaluate_expression( "norm", make_pair("bra", LL), make_pair("ket", MM ) );
   
    DataType pt_energy = this->get_scalar_result( "pt_energy", make_pair("state", NN) );
    //pt2energy = energy_[NN]+(xms_ham)(NN,NN) - e_shift_*norm;
  }

//  shared_ptr<Matrix> h_eff = make_shared<Matrix>( ref_space_dim, ref_space_dim );
 
  // Loop through projected ham terms 
  for (int QQ = 0; QQ != ref_space_dim; ++QQ) {
    for (int MM = 0; MM != ref_space_dim; ++MM){ 
      this->evaluate_expression( "norm", make_pair("ket", QQ), make_pair("bra", MM ) );
      this->evaluate_expression( "proj_ham", make_pair("ket", QQ), make_pair("bra", MM ) );
    }
    for (int LL = 0; LL != ref_space_dim; ++LL) {
      if (QQ == LL) {
        (*h_eff)(QQ, QQ);// = pt2energy_[QQ];
      } else {
        //H^{ref}_{LQ} + 1/2 \sum_{M} [ < M  | T^{\dagger}_{LM}  H  |  Q > + < L |H T_{QM} | M > ]
        (*h_eff)( LL , QQ );// = dot_product_transpose(proj_ham, pt_amps) + (xms_ham)( LL, QQ );
      }
    }
  }
}
///////////////////////////////////////////////////////////////
template class Equation_Computer_LinearRM<double>;
///////////////////////////////////////////////////////////////
