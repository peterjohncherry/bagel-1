//
// Newint - Parallel electron correlation program.
// Filename: geometry.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <src/scf/geometry.h>
#include <src/scf/atommap.h>
#include <src/stackmem.h>

using namespace std;

typedef std::shared_ptr<Shell> RefShell;
typedef std::shared_ptr<Atom> RefAtom;

extern StackMem* stack;

Geometry::Geometry(const std::shared_ptr<const InputData> inpt)
  : spherical_(true), input_(""), lmax_(0), level_(0) {

  multimap<string, string> geominfo = inpt->get_input("molecule");

  schwarz_thresh_ = read_input<double>(geominfo, "schwarz_thresh", 1.0e-12); 
  const double thresh_overlap = read_input<double>(geominfo, "thresh_overlap", 1.0e-8); 

  // cartesian or not.
  const bool cart = read_input<bool>(geominfo, "cartesian", false); 
  if (cart) {
    cout << "  Cartesian basis functions are used" << endl;
    spherical_ = false;
  }

  // basis
  basisfile_ = read_input<string>(geominfo, "basis", "");
  if (basisfile_ == "") throw runtime_error("There is no basis specification");

  auxfile_ = read_input<string>(geominfo, "df_basis", "");

  // symmetry
  symmetry_ = read_input<string>(geominfo, "symmetry", "c1");

  const bool angstrom = read_input<bool>(geominfo, "angstrom", false);

  AtomMap tmp;
  map<string, int> amap = tmp.atommap;

  // read geometry
  nbasis_ = 0;
  naux_ = 0;
  nele_ = 0;
  nfrc_ = 0;

  pair<multimap<string,string>::const_iterator, multimap<string,string>::const_iterator> bound = geominfo.equal_range("atom");
  for (auto iter = bound.first; iter != bound.second; ++iter) { 
    boost::smatch what;
    const boost::regex atom_reg("\\(\\s*([A-Za-z]+),\\s*([0-9\\.-]+),\\s*([0-9\\.-]+),\\s*([0-9\\.-]+)\\s*\\)");
    auto start = iter->second.begin();
    auto end = iter->second.end();
    if (regex_search(start, end, what, atom_reg)) {
      const string aname(what[1].first, what[1].second);
      const string x_str(what[2].first, what[2].second);
      const string y_str(what[3].first, what[3].second);
      const string z_str(what[4].first, what[4].second);
      vector<double> positions;
      const double prefac = angstrom ? 1.0/0.529177249 : 1.0 ;
      positions.push_back(boost::lexical_cast<double>(x_str)*prefac);
      positions.push_back(boost::lexical_cast<double>(y_str)*prefac);
      positions.push_back(boost::lexical_cast<double>(z_str)*prefac);

      {
        RefAtom catom(new Atom(spherical_, aname, positions, basisfile_));

        map<string, int>::const_iterator aiter = amap.find(aname);
        assert(aiter != amap.end());
        nele_ += aiter->second;

        const vector<RefShell> tmp = catom->shells(); 
        int cc = 0;
        vector<int> coffsets;
        for (vector<RefShell>::const_iterator iter = tmp.begin(); iter != tmp.end(); ++iter) {
          coffsets.push_back(nbasis_ + cc);
          const int ang = (*iter)->angular_number();
          const int angsize = spherical_ ? (2 * ang + 1) : (ang + 1) * (ang + 2) / 2;
          cc += angsize * (*iter)->num_contracted();
        }

        lmax_ = max(lmax_, catom->lmax());
        nbasis_ += catom->nbasis();
        offsets_.push_back(coffsets);
        atoms_.push_back(catom); 
      }
      if (!auxfile_.empty()){
        RefAtom catom(new Atom(spherical_, aname, positions, auxfile_));

        map<string, int>::const_iterator aiter = amap.find(aname);
        assert(aiter != amap.end());

        const vector<RefShell> tmp = catom->shells(); 
        int cc = 0;
        vector<int> coffsets;
        for (vector<RefShell>::const_iterator iter = tmp.begin(); iter != tmp.end(); ++iter) {
          coffsets.push_back(naux_ + cc);
          const int ang = (*iter)->angular_number();
          const int angsize = spherical_ ? (2 * ang + 1) : (ang + 1) * (ang + 2) / 2;
          cc += angsize * (*iter)->num_contracted();
        }

        aux_lmax_ = max(lmax_, catom->lmax());
        naux_ += catom->nbasis();
        aux_offsets_.push_back(coffsets);
        aux_atoms_.push_back(catom); 
      }
    } else {
      throw runtime_error("One of the atom lines is corrupt");
    } 
  }

  print_atoms();
  nuclear_repulsion_ = compute_nuclear_repulsion();

  // symmetry set-up
  plist_ = std::shared_ptr<Petite>(new Petite(atoms_, symmetry_));
  nirrep_ = plist_->nirrep();
  // Misc
  aux_merged_ = false;

  cout << endl;
  cout << "  Number of basis functions: " << setw(8) << nbasis() << endl;
  cout << "  Number of electrons      : " << setw(8) << nele() << endl << endl;
  if (!auxfile_.empty()) cout << "  Number of auxiliary basis functions: " << setw(8) << naux() << endl << endl;

  if (!auxfile_.empty()) {
    cout << "  Since a DF basis is specified, we compute 2- and 3-index integrals:" << endl;
    cout << "    o Being stored without compression. Storage requirement is "
         << setprecision(3) << static_cast<size_t>(naux_)*nbasis()*nbasis()*8.e-9 << " GB" << endl;
    const int t = ::clock();
    df_ = form_fit<ERIFit>(thresh_overlap, true); // true means we construct J^-1/2

    cout << "        elapsed time:  " << setw(10) << setprecision(2) << (::clock() - t)/static_cast<double>(CLOCKS_PER_SEC) << " sec." << endl << endl; 
  }
}

static vector<double> vec3(const double i, const double j, const double k) {
  vector<double> out(3); out[0] = i; out[1] = j; out[2] = k; return out; 
};

Geometry::Geometry(const Geometry& o, const std::shared_ptr<GradFile> displ, const std::shared_ptr<const InputData> inpt)
  : spherical_(o.spherical_), input_(o.input_), aux_merged_(o.aux_merged_), nbasis_(o.nbasis_), nele_(o.nele_), nfrc_(o.nfrc_), naux_(o.naux_),
    lmax_(o.lmax_), aux_lmax_(o.aux_lmax_), offsets_(o.offsets_), aux_offsets_(o.aux_offsets_), level_(o.level_), basisfile_(o.basisfile_),
    auxfile_(o.auxfile_), symmetry_(o.symmetry_), schwarz_thresh_(o.schwarz_thresh_), gamma_(o.gamma_) { 

  // first construct atoms using displacements
  double* disp = displ->data();
  for (auto i = o.atoms_.begin(), j = o.aux_atoms_.begin(); i != o.atoms_.end(); ++i, ++j, disp += 3) {
    vector<double> cdispl = vec3(*disp, *(disp+1), *(disp+2));
    atoms_.push_back(shared_ptr<Atom>(new Atom(**i, cdispl)));
    aux_atoms_.push_back(shared_ptr<Atom>(new Atom(**j, cdispl)));
  }

  plist_ = std::shared_ptr<Petite>(new Petite(atoms_, symmetry_));
  nirrep_ = plist_->nirrep();
  if (nirrep_ != o.nirrep_) throw logic_error("opt with symmetry is not tested well yet"); 
  nuclear_repulsion_ = compute_nuclear_repulsion();

  if (!auxfile_.empty()) {
    const double thresh_overlap = read_input<double>(inpt->get_input("molecule"), "thresh_overlap", 1.0e-8); 
    cout << "  Since a DF basis is specified, we compute 2- and 3-index integrals:" << endl;
    cout << "    o Being stored without compression. Storage requirement is "
         << setprecision(3) << static_cast<size_t>(naux_)*nbasis()*nbasis()*8.e-9 << " GB" << endl;
    const int t = ::clock();
    df_ = form_fit<ERIFit>(thresh_overlap, true); // true means we construct J^-1/2

    cout << "        elapsed time:  " << setw(10) << setprecision(2) << (::clock() - t)/static_cast<double>(CLOCKS_PER_SEC) << " sec." << endl << endl; 
  }
}


Geometry::Geometry(const string s, const int levl)
  : spherical_(true), input_(s), lmax_(0), level_(levl) {

  // open input file
  ifstream ifs;
  ifs.open(input_.c_str());
  assert(ifs.is_open());

  boost::regex frozen_str("frozen"); 
  bool frozen = false;
  while(!ifs.eof()) {
    string sline;
    getline(ifs, sline);
    if(sline.empty()) continue;
    string::const_iterator start = sline.begin();
    string::const_iterator end = sline.end();
    boost::smatch what;
    if (regex_search(start, end, what, frozen_str)) {
      frozen = true; 
      break;
    }
  }
  ifs.clear(); 
  ifs.seekg(0);

  boost::regex gamma_str("gamma");
  boost::regex gamma_num("[0-9\\.eE\\+-]+");
  double gamma = 1.5;
  while(!ifs.eof()) {
    string sline;
    getline(ifs, sline);
    if(sline.empty()) continue;
    string::const_iterator start = sline.begin();
    string::const_iterator end = sline.end();
    boost::smatch what;
    if (regex_search(start, end, what, gamma_str)) {
      start = what[0].second;
      if (regex_search(start, end, what, gamma_num)) {
        const string gamma_str(what[0].first, what[0].second);
        gamma = boost::lexical_cast<double>(gamma_str);
        break;
      }
    }
  }
  gamma_ = gamma;
  ifs.clear();
  ifs.seekg(0); 

  // read basis file
  // this will be read from the input
  string basis_name("Basis");
  for (int i = 0; i != level_; ++i) basis_name = "G" + basis_name;
  basis_name = "\\b" + basis_name;
  const boost::regex basis_reg(basis_name);

  while(!ifs.eof()) {
    string sline;
    getline(ifs, sline);
    if(sline.empty()) continue;
    string::const_iterator start = sline.begin();
    string::const_iterator end = sline.end();
    boost::smatch what;
    if (regex_search(start, end, what, basis_reg)) {
      start = what[0].second;
      boost::regex car_reg("cartesian");
      if (regex_search(start, end, what, car_reg)) spherical_ = false; 
      if (!spherical_) cout << "  Cartesian basis functions are used" << endl;
      getline(ifs, sline);
      if (sline.empty()) continue;
      start = sline.begin();
      end = sline.end();
      const boost::regex reg("(\\S+)");
      const bool found = boost::regex_search(start, end, what, reg);
      const string tmpstr(what[1].first, what[1].second);
      basisfile_ = tmpstr; 

      start = what[0].second;
      if(regex_search(start, end, what, reg)) {
        const string auxstr(what[1].first, what[1].second);
        auxfile_ = auxstr;
      }
      
      break;
    }
  }

  assert(!basisfile_.empty());
  ifs.clear();
  ifs.seekg(0);

  AtomMap tmp;
  map<string, int> amap = tmp.atommap;

  // read geometry
  nbasis_ = 0;
  naux_ = 0;
  nele_ = 0;
  nfrc_ = 0;

  const boost::regex mole_reg("Molecule");
  while(!ifs.eof()){
    string sline;
    getline(ifs, sline); 
    if (sline.empty()) continue; 
    string::const_iterator start = sline.begin();
    string::const_iterator end = sline.end();
    boost::smatch what;
    if (regex_search(start, end, what, mole_reg)) {
      start = what[0].second;
      boost::regex sym_reg("[cdCD][1-2]?[vdshVDSHiI]?");
      if (regex_search(start, end, what, sym_reg)) {
        string symtmp(what[0].first, what[0].second);
        transform(symtmp.begin(), symtmp.end(), symtmp.begin(),(int (*)(int))std::tolower);
        symmetry_ = symtmp; 
      } else {
        string symtmp("c1");
        symmetry_ = symtmp;
      }
      const boost::regex atom_reg("Atom\\s*\\(\\s*([A-Za-z]+),\\s*([0-9\\.-]+),\\s*([0-9\\.-]+),\\s*([0-9\\.-]+)\\s*\\)");
      while (1) {
        string atomline; 
        getline(ifs, atomline); 
        start = atomline.begin();
        end = atomline.end();
        if (regex_search(start, end, what, atom_reg)) {
          const string aname(what[1].first, what[1].second);
          const string x_str(what[2].first, what[2].second);
          const string y_str(what[3].first, what[3].second);
          const string z_str(what[4].first, what[4].second);
          vector<double> positions;
          positions.push_back(boost::lexical_cast<double>(x_str));
          positions.push_back(boost::lexical_cast<double>(y_str));
          positions.push_back(boost::lexical_cast<double>(z_str));

          {
            RefAtom catom(new Atom(spherical_, aname, positions, basisfile_));

            map<string, int>::const_iterator aiter = amap.find(aname);
            assert(aiter != amap.end());
            nele_ += aiter->second;
            if (aiter->second > 2 && frozen) nfrc_ += 2;

            const vector<RefShell> tmp = catom->shells(); 
            int cc = 0;
            vector<int> coffsets;
            for (vector<RefShell>::const_iterator iter = tmp.begin(); iter != tmp.end(); ++iter) {
              coffsets.push_back(nbasis_ + cc);
              const int ang = (*iter)->angular_number();
              const int angsize = spherical_ ? (2 * ang + 1) : (ang + 1) * (ang + 2) / 2;
              cc += angsize * (*iter)->num_contracted();
            }

            lmax_ = max(lmax_, catom->lmax());
            nbasis_ += catom->nbasis();
            offsets_.push_back(coffsets);
            atoms_.push_back(catom); 
          }
          if (!auxfile_.empty()){
            RefAtom catom(new Atom(spherical_, aname, positions, auxfile_));

            map<string, int>::const_iterator aiter = amap.find(aname);
            assert(aiter != amap.end());

            const vector<RefShell> tmp = catom->shells(); 
            int cc = 0;
            vector<int> coffsets;
            for (vector<RefShell>::const_iterator iter = tmp.begin(); iter != tmp.end(); ++iter) {
              coffsets.push_back(naux_ + cc);
              const int ang = (*iter)->angular_number();
              const int angsize = spherical_ ? (2 * ang + 1) : (ang + 1) * (ang + 2) / 2;
              cc += angsize * (*iter)->num_contracted();
            }

            aux_lmax_ = max(lmax_, catom->lmax());
            naux_ += catom->nbasis();
            aux_offsets_.push_back(coffsets);
            aux_atoms_.push_back(catom); 
          }
        } else {
          break; 
        }
      }
      break;
    } 
  }
  ifs.close();

  if (level_ == 0) print_atoms();
  nuclear_repulsion_ = compute_nuclear_repulsion();

  // symmetry set-up
  std::shared_ptr<Petite> tmpp(new Petite(atoms_, symmetry_));
  plist_ = tmpp;
  nirrep_ = plist_->nirrep();
  // Misc
  aux_merged_ = false;

  cout << endl;
  cout << "  Number of basis functions: " << setw(8) << nbasis() << endl;
  cout << "  Number of electrons      : " << setw(8) << nele() << endl << endl;
  if (!auxfile_.empty()) cout << "  Number of auxiliary basis functions: " << setw(8) << naux() << endl << endl;
  cout << endl;
}


Geometry::~Geometry() {
}


const double Geometry::compute_nuclear_repulsion() {
  double out = 0.0;
  for (vector<RefAtom>::const_iterator iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    const vector<double> tmp = (*iter)->position();
    const double c = static_cast<double>((*iter)->atom_number());
    for (vector<RefAtom>::const_iterator titer = iter + 1; titer != atoms_.end(); ++titer) {
      const RefAtom target = *titer;   
      const double x = target->position(0) - tmp[0];
      const double y = target->position(1) - tmp[1];
      const double z = target->position(2) - tmp[2];
      const double dist = ::sqrt(x * x + y * y + z * z);
      const double charge = c * static_cast<double>(target->atom_number()); 
      out += charge / dist;
    }
  }
  return out;
}


void Geometry::print_atoms() const {
  cout << "  === Geometry ===" << endl << endl;
  cout << "  Symmetry: " << symmetry() << endl;
  cout << endl;

  for (vector<RefAtom>::const_iterator iter = atoms_.begin(); iter != atoms_.end(); ++iter)
    (*iter)->print();
  cout << endl;

}


int Geometry::num_count_ncore() {
  int out = 0;
  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    if ((*iter)->atom_number() >= 2 && (*iter)->atom_number() <= 10) out += 2; 
    if ((*iter)->atom_number() > 10) throw logic_error("needs to modify Geometry::count_num_ncore for atoms beyond Ne"); // TODO
  }
  nfrc_ = out;
  return out;
}

int Geometry::num_count_full_valence_nocc() {
  int out = 0;
  for (auto iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
    if ((*iter)->atom_number() < 2) out += 1;
    if ((*iter)->atom_number() >= 2 && (*iter)->atom_number() <= 10) out += 5; 
    if ((*iter)->atom_number() > 10) throw logic_error("needs to modify Geometry::num_count_full_valence_nocc for atoms beyond Ne"); // TODO
  }
  return out;
};


const vector<double> Geometry::compute_grad_vnuc() const {
  // the derivative of Vnuc
  vector<double> grad(natom()*3);
  fill(grad.begin(), grad.end(), 0.0);
  auto giter = grad.begin();
  for (auto aiter = atoms_.begin(); aiter != atoms_.end(); ++aiter, giter+=3) {
    const double ax = (*aiter)->position(0);
    const double ay = (*aiter)->position(1);
    const double az = (*aiter)->position(2);
    const double ac = (*aiter)->atom_number();
    for (auto biter = atoms_.begin(); biter != atoms_.end(); ++biter) {
      if (aiter == biter) continue;
      const double bx = (*biter)->position(0);
      const double by = (*biter)->position(1);
      const double bz = (*biter)->position(2);
      const double c = (*biter)->atom_number() * ac;
      const double dist = sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz));
      *(giter+0) += c*(bx-ax)/(dist*dist*dist);
      *(giter+1) += c*(by-ay)/(dist*dist*dist);
      *(giter+2) += c*(bz-az)/(dist*dist*dist);
    }
  }
  return grad;
}


void Geometry::merge_obs_aux() {
  aux_merged_ = true;
  atoms_.insert(atoms_.end(), aux_atoms_.begin(), aux_atoms_.end());
  for (std::vector<std::vector<int> >::iterator iter = aux_offsets_.begin(); iter != aux_offsets_.end(); ++iter) {
    for (std::vector<int>::iterator citer = iter->begin(); citer != iter->end(); ++citer) {
      *citer += nbasis_;
    }
  }
  offsets_.insert(offsets_.end(), aux_offsets_.begin(), aux_offsets_.end());
  nbasis_ += naux_;
}


shared_ptr<GradFile> Geometry::gradfile() const {
  vector<double> tmp;
  for (auto i = atoms_.begin(); i != atoms_.end(); ++i) {
    tmp.push_back((*i)->position(0));
    tmp.push_back((*i)->position(1));
    tmp.push_back((*i)->position(2));
  }
  return shared_ptr<GradFile>(new GradFile(tmp));
}
