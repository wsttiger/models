/*
 * heisen.h
 *
 *  Created on: Aug 22, 2011
 *      Author: wsttiger
 */

#ifndef HEISEN_H_
#define HEISEN_H_

#include <iostream>
#include <fstream>
#include <bitset>
#include <map>
#include <cassert>

#include "matrix.h"

using std::cout;
using std::endl;
using std::map;
using std::pair;
using std::string;
using std::bitset;
using std::vector;
using std::complex;
using std::ifstream;

#define NSIZE 4 

typedef bitset<NSIZE> ketT;

//*****************************************************************************
void print_state(const ketT& ket) {
  printf("|%s>\n",ket.to_string().c_str());
}
//*****************************************************************************

//*****************************************************************************
int cntbits(int s1, int s2, ketT state) {
  if (s1 == s2) return 0;

  int sum = 0;
  if (s1 > s2) {
    for (int i = s2;i < s1; i++) {
      if (state[i]) sum++;
    }
  }
  else if (s1 < s2) {
    for (int i = s1;i < s2; i++) {
      if (state[i]) sum++;
    }
  }
  return sum;
}
//*****************************************************************************

//*****************************************************************************
class States {
public:
  States () {}

  States(int nup, int nsites)
    : nsites(nsites), nup(nup) {

    printf("states constructor ... begin\n");
    int tnstates = 1 << nsites;

    // create states for up and down spins
    for (int i = 0; i < tnstates; i++) {
      ketT tstate(i);
      if (nup < 0) {
        states.push_back(tstate);
      }
      else if (cntbits(0,nsites, tstate) == nup) {
        states.push_back(tstate);
      }
    }

    for (unsigned int i = 0; i < states.size(); i++) {
      umap.insert(pair<unsigned long,unsigned long>(states[i].to_ulong(),i));
    }

    nst= states.size();
    printf("states constructor ... end\n");
  }

  int get_nst() const {return nst;} 

  ketT get_state(int indx) const {
    return states[indx];
  }

  unsigned long get_state_index(const ketT s) const {
    return umap.find(s.to_ulong())->second;
  }

  int get_nup() {return nup;}

  int get_ndown() {return nsites-nup;}

  void print() {
    printf("nsites:  %d     nup:  %d\n", nsites, nup);
    for (unsigned int ist = 0; ist < states.size(); ist++) {
      printf("%d   |%s>    %ld\n",ist,
          states[ist].to_string().c_str(),states[ist].to_ulong());
    }
  }

private:
  vector< ketT > states;
  vector< ketT > dstates;
  map<unsigned long,unsigned long> umap;

  int nsites;
  int nup;
  int nst;

};
//*****************************************************************************

//*****************************************************************************
struct Bond {
  int i, j;
  double J;
  Bond() {}
  Bond(int i, int j, double J) 
   : i(i), j(j), J(J) {}
};
typedef vector<Bond> latticeT;
//*****************************************************************************

//*****************************************************************************
class Sector {
public:
  Sector() {}

  Sector(latticeT lattice, int sites, int nup) 
   : lattice(lattice) {
    states = States(nup, sites);
    int nst = states.get_nst();
  }

  virtual vector<double> apply_works_no_debug(const vector<double>& svec) {
    vector<double> rvec(svec.size());
    for (auto ist = 0; ist < svec.size(); ist++) {
      ketT st = states.get_state(ist);
      double val = svec[ist];
      if (std::abs(val) > 1e-15) {
        // loop over lattice sites
        for (auto is = 0; is < lattice.size(); is++) {
          Bond bond = lattice[is];
          int i = bond.i; int j = bond.j; double J = bond.J;
          // diagonal
          if (st[i] == st[j])
            rvec[ist] += -0.25*J*val;
          // off-diagonal
          else {
            rvec[ist] += 0.25*J*val;
            ketT stij = ketT(st);
            stij[i] = !st[i]; stij[j] = !st[j];
            unsigned long lstij = states.get_state_index(stij);
            rvec[lstij] -= 0.5*J*val;
          }
        } 
      }
    }
    return rvec;
  }

  virtual vector<double> apply_doesnt_work(const vector<double>& svec) {
    vector<double> rvec(svec.size());
    for (auto ist = 0; ist < svec.size(); ist++) {
      ketT st = states.get_state(ist);
      double val = svec[ist];
      if (std::abs(val) > 1e-15) {
        // loop over lattice sites
        for (auto is = 0; is < lattice.size(); is++) {
          Bond bond = lattice[is];
          int i = bond.i; int j = bond.j; double J = bond.J;
          // off-diagonal
          ketT spmop(0);
          spmop[i] = true; spmop[j] = true;
          int matchbits = cntbits(0,NSIZE,spmop & st);
          if (matchbits == 1) {
            rvec[ist] += 0.25*J*val;
            ketT stij = ketT(st);
            stij ^= spmop;
            unsigned long lstij = states.get_state_index(stij);
            rvec[lstij] -= 0.5*J*val;
          }
          else if (matchbits == 2) {
            rvec[ist] += -0.25*J*val;
          }
        } 
      }
    }
    return rvec;
  }

  virtual vector<double> apply_works(const vector<double>& svec) {
    vector<double> rvec(svec.size());
    for (auto ist = 0; ist < svec.size(); ist++) {
      ketT st = states.get_state(ist);
      double val = svec[ist];
      printf("apply(): ist = %d\n",ist);
      if (std::abs(val) > 1e-15) {
        // loop over lattice sites
        for (auto is = 0; is < lattice.size(); is++) {
          Bond bond = lattice[is];
          int i = bond.i; int j = bond.j; double J = bond.J;
          // off-diagonal
          printf("Bond (i,j): %d     %d     J = %4.2f\n",i,j,J);
          if (st[i] != st[j]) {
            ketT stij = ketT(st);
            stij[i] = !st[i]; stij[j] = !st[j];
            unsigned long lstij = states.get_state_index(stij);
            std::cout << lstij << "|" << st << ">    ----->       | " << stij << ">" << std::endl;
            rvec[lstij] -= 0.5*J*val;
            printf("rvec[%ld]:  %5.3f\n", lstij, rvec[lstij]);
          }
          // diagonal
          printf("J:  %10.5f     val:  %10.5f\n", J, val);
          std::cout << "(before)   |" << st << ">     --> " << rvec[ist] << std::endl;
          rvec[ist] += (st[i] == st[j]) ? -0.25*J*val : 0.25*J*val;
          printf("Diagonal: \n");
          std::cout << "(after)    |" << st << ">     --> " << rvec[ist] << std::endl;
          printf("\n\n");
        } 
      }
    }
    return rvec;
  }

  vector<double> make_matrix()
  {
    int nst= states.get_nst();
    vector<double> rmat(nst*nst,0.0);

    for (int i = 0; i < nst; i++) {
      vector<double> v1(nst,0.0);
      v1[i] = 1.0;
      for (int j = 0; j < nst; j++) {
        vector<double> v2(nst,0.0);
        v2[j] = 1.0;

        vector<double> rvec = apply_works_no_debug(v1);
        double t1 = dotvec(rvec,v2);
        rmat[i*nst+j] = t1;
      }
    }

    // make sure that the matrix is symmetric
    double tol = 1e-10;
    for (int i = 0; i < nst; i++) {
      for (int j = 0; j < i; j++) {
        if (fabs(rmat[i*nst+j]-rmat[j*nst+i]) > tol)
        {
          printf("ERROR: (make_matrix()) elements [%3d,%3d] != [%3d,%3d]\n", i,j,j,i);
          exit(EXIT_FAILURE);
        }
      }
    }

    return rmat;
  }

  int get_nst() const {return states.get_nst();}
private:
  States states;
  latticeT lattice;
};
//*****************************************************************************

//*****************************************************************************
class HeisenCalculation {
public:
  
//  HeisenCalculation(const string& sfile) {
//    ifstream astream(sfile.c_str());
//    int nup;
//    double J_tmp;
//    int nx;
//    int ny;
//    int latticetype;
//    std::string lattice_file;
//
//    while(!astream.eof()) {
//      string s1;
//      astream >> s1;
//      s1.erase(remove_if(s1.begin(), s1.end(), isspace), s1.end());
//      if (s1.compare("") == 0) {}
//      else if (s1.compare("J") == 0) {
//        astream >> J_tmp;
//      }
//      else if (s1.compare("nup") == 0) {
//        astream >> nup;
//      }
//      else if (s1.compare("latticetype") == 0) {
//        astream >> latticetype;
//      }
//      else if (s1.compare("nx") == 0) {
//        astream >> nx;
//      }
//      else if (s1.compare("ny") == 0) {
//        astream >> ny;
//      }
//      else {
//        printf("ERROR reading input file: %s\n\n", s1.c_str());
//        exit(EXIT_FAILURE);
//      }
//    }
//
//    if (latticetype == 0) {
//      if (ny == 1)
//      {
//        //create_1d_lattice();
//      }
//      else {
//        //create_2d_square_lattice();
//      }
//    }
//    else {
//      printf("lattice file not implemented yet!\n\n");
//      assert(false);
//    }
//
//    states = States(nup,nx*ny);
//
//    astream.close();
//  }

  HeisenCalculation(int sites) {
    create_ring(sites);
  }

  void create_ring(int sites) {
//    double J = 1.0;
//    for (int i = 0; i < sites-1; i++) {
//      lattice.push_back(Bond(i,i+1,J));
//    }
//    lattice.push_back(Bond(sites-1,0,J));
//    sectors.push_back(Sector(lattice,sites,-1));

    double J = 1.0;
    for (int i = 0; i < sites-1; i++) {
      lattice.push_back(Bond(i,i+1,J));
    }
    lattice.push_back(Bond(sites-1,0,J));
    for (int i = 0; i < sites; i++) {
      sectors.push_back(Sector(lattice,sites,i));
    }
  }

  vector<double> eigenvalues() {
    vector<double> e;
    for (unsigned int is = 0; is < sectors.size(); is++) {
      Sector s = sectors[is];
      int nst = s.get_nst();
      vector<double> mat = s.make_matrix(); 
      vector<double> es(nst,0.0);
      vector<double> ev(nst*nst,0.0);
      diag_matrix(mat,nst,es,ev);
      std::copy(es.begin(),es.end(),std::back_inserter(e));
    }
    std::sort(e.begin(),e.end(),[](const double& a, const double& b) {return a < b;});
    return e;
  }

private:
  latticeT lattice;
  std::vector<Sector> sectors;
};
//*****************************************************************************

#endif
