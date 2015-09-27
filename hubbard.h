/*
 * hubbard.h
 *
 *  Created on: Aug 22, 2011
 *      Author: wsttiger
 */

#ifndef HUBBARD_H_
#define HUBBARD_H_

#include <iostream>
#include <fstream>
#include <bitset>
#include <omp.h>

#include "matrix.h"
#include "lanczos.h"

using std::cout;
using std::endl;
using std::map;
using std::pair;
using std::string;
using std::bitset;
using std::vector;
using std::complex;
using std::ifstream;

#define NSIZE 26

//*****************************************************************************
int cntbits(int s1, int s2, bitset<NSIZE> state)
{
  if (s1 == s2) return 0;

  int sum = 0;
  if (s1 > s2)
  {
    for (int i = s2;i < s1; i++)
    {
      if (state[i]) sum++;
    }
  }
  else if (s1 < s2)
  {
    for (int i = s1;i < s2; i++)
    {
      if (state[i]) sum++;
    }
  }
  return sum;
}
//*****************************************************************************

//*****************************************************************************
double total_U_for_state(const bitset<NSIZE>& state1, const bitset<NSIZE>& state2,
    const std::vector<double> U)
{
  bitset<NSIZE> state = state1 & state2;
  double totalU = 0.0;
  for (unsigned int i = 0; i < U.size(); i++)
  {
    if (state[i]) totalU += U[i];
  }
  return totalU;
}
//*****************************************************************************

//*****************************************************************************
class States
{
public:
  States(int nup, int ndown, int nsites)
  {
    this->nup = nup;
    this->ndown = ndown;

    int tnstates = 1 << nsites;
    int usum = 0;
    int dsum = 0;
    nup2 = nchoosek(nsites,nup);
    ndown2 = nchoosek(nsites,ndown);

    // create states for up and down spins
    for (int i = 0; i < tnstates; i++)
    {
      bitset<NSIZE> tstate(i);
      if (cntbits(0,nsites, tstate) == nup)
      {
        ustates.push_back(tstate);
        usum++;
      }
      if (cntbits(0,nsites,tstate) == ndown)
      {
        dstates.push_back(tstate);
        dsum++;
      }
      if (usum > nup2)
      {
        cout << "[[Warning]] create_states: usum > nup2" << endl;
        cout << "            usum: " << usum << "  nup2: " << nup2 << endl;
        exit(EXIT_FAILURE);
      }
      if (dsum > ndown2)
      {
        cout << "[[Warning]] create_states: dsum > ndown2" << endl;
        cout << "            dsum: " << usum << "  ndown2: " << nup2 << endl;
        exit(EXIT_FAILURE);
      }
    }

    // create map where the first item is the state of a single spin state
    // i.e. 00000000000000000000000011101000 and give the associated index
    for (unsigned int i = 0; i < ustates.size(); i++)
    {
      umap.insert(pair<unsigned long,unsigned long>(ustates[i].to_ulong(),i));
    }
    for (unsigned int i = 0; i < dstates.size(); i++)
    {
      dmap.insert(pair<unsigned long,unsigned long>(dstates[i].to_ulong(),i));
    }

    nstates = ustates.size()*dstates.size();
  }

  int get_nstates() {return nstates;}

  bitset<NSIZE> get_ustate(int indx)
  {
    return ustates[indx];
  }

  bitset<NSIZE> get_dstate(int indx)
  {
    return dstates[indx];
  }

  unsigned long get_ustate_index(const bitset<NSIZE> us)
  {
    return umap.find(us.to_ulong())->second;
  }

  unsigned long get_dstate_index(const bitset<NSIZE> ds)
  {
    return dmap.find(ds.to_ulong())->second;
  }

  int get_nup() {return nup;}

  int get_ndown() {return ndown;}

  int get_nup2() {return nup2;}

  int get_ndown2() {return ndown2;}

  void print_states2()
  {
    printf("\nStates: %d up states     %d down states\n",nup,ndown);
    for (unsigned int ist = 0; ist < ustates.size()*dstates.size(); ist++)
    {
      int ust = ist / ndown2;
      int dst = ist % ndown2;
      printf("%d    %d    %d    [%s , %s]\n",ist,ust,dst,
          ustates[ust].to_string().c_str(),dstates[dst].to_string().c_str());
    }
  }


private:
  vector< bitset<NSIZE> > ustates;
  vector< bitset<NSIZE> > dstates;
  map<unsigned long,unsigned long> umap;
  map<unsigned long,unsigned long> dmap;

  int nup;
  int ndown;
  int nup2;
  int ndown2;
  int nstates;

};
//*****************************************************************************

////*****************************************************************************
//class IModel
//{
//public:
//  virtual vector<double> apply(const vector<double>& v) = 0;
//  virtual vector<double> apply(const vector<double>& v, States* s) = 0;
//  virtual int get_nstates() const = 0;
//};
////*****************************************************************************

//*****************************************************************************
//class HubbardCalculation : public IModel
class HubbardCalculation 
{
private:
  int nx;
  int ny;
//  double U;
//  double t;
  States* states;
  double tol;
  bool bdebug;

  // Ground state vector
  vector<double> gsv;

  struct neighbor_2d
  {
    int me;
    int n1;

    neighbor_2d(int me_, int n1_)
    {
      me = me_;
      n1 = n1_;
    }
  };

  vector<neighbor_2d> lattice;
  vector<double> t;
  vector<double> U;

public:

  HubbardCalculation(const string& sfile)
  {
    ifstream astream(sfile.c_str());
    int nup;
    int ndown;
    double t_tmp;
    double U_tmp;
    double nx_tmp;
    double ny_tmp;
    int latticetype;
    std::string lattice_file;

    while(!astream.eof())
    {
      string s1;
      astream >> s1;
      s1.erase(remove_if(s1.begin(), s1.end(), isspace), s1.end());
      if (s1.compare("") == 0) {}
      else if (s1.compare("t") == 0)
      {
        astream >> t_tmp;
      }
      else if (s1.compare("U") == 0)
      {
        astream >> U_tmp;
      }
      else if (s1.compare("nup") == 0)
      {
        astream >> nup;
      }
      else if (s1.compare("ndown") == 0)
      {
        astream >> ndown;
      }
      else if (s1.compare("latticetype") == 0)
      {
        astream >> latticetype;
      }
      else if (s1.compare("nx") == 0)
      {
        astream >> nx_tmp;
      }
      else if (s1.compare("ny") == 0)
      {
        astream >> ny_tmp;
      }
//      else if (s1.compare("") == 0)
//      {
//        astream >> ny_tmp;
//      }
      else
      {
        printf("ERROR reading input file: %s\n\n", s1.c_str());
        exit(EXIT_FAILURE);
      }
    }

    if (latticetype == 0)
    {
      nx = nx_tmp;
      ny = ny_tmp;

      if (ny == 1)
      {
        create_1d_lattice();
      }
      else
      {
        create_2d_square_lattice();
      }

      for (unsigned int il = 0; il < lattice.size(); il++)
      {
        t.push_back(t_tmp);
      }
      for (unsigned int isite = 0; isite < nx*ny; isite++)
      {
        U.push_back(U_tmp);
      }
    }
    else
    {
      printf("lattice file not implemented yet!\n\n");
      assert(false);
    }

    tol = 1e-8;
    bdebug = false;

    states = new States(nup,ndown,nx*ny);
    //states->print_states2();

    astream.close();
  }

  HubbardCalculation(int nx, int ny, int nup, int ndown, double U_, double t_ = 1.0)
   : nx(nx), ny(ny)
  {
    if (ny == 1)
    {
      create_1d_lattice();
    }
    else
    {
      create_2d_square_lattice();
    }
    tol = 1e-8;
    bdebug = false;

    states = new States(nup,ndown,nx*ny);

  }

  ~HubbardCalculation()
  {
    delete states;
  }

  virtual States* get_states()
  {
    return states;
  }

  virtual int get_nstates() const
  {
    return states->get_nstates();
  }

//  void print_states()
//  {
//    printf("\nStates: %d up states     %d down states\n",nup,ndown);
//    int i = 0;
//    for (unsigned int ust = 0; ust < ustates.size(); ust++)
//    {
//      for (unsigned int dst = 0; dst < dstates.size(); dst++,i++)
//      {
//        printf("%5d [%s , %s]\n",i,ustates[ust].to_string().c_str(),dstates[dst].to_string().c_str());
//      }
//    }
//  }

//  void print_states2()
//  {
//    printf("\nStates: %d up states     %d down states\n",nup,ndown);
//    for (unsigned int ist = 0; ist < ustates.size()*dstates.size(); ist++)
//    {
//      int ust = ist / ndown2;
//      int dst = ist % ndown2;
//      printf("%d    %d    [%s , %s]\n",ust,dst,
//          ustates[ust].to_string().c_str(),dstates[dst].to_string().c_str());
//    }
//  }

  void print_state(const vector<double>& vec)
  {
    printf("\nState:\n");
    for (unsigned int ist = 0; ist < states->get_nstates(); ist++)
    {
      unsigned int ndown2 = states->get_ndown2();
      unsigned int ust = ist / ndown2;
      unsigned int dst = ist % ndown2;
      printf("%10.5f  [%s , %s]\n",vec[ist],
          states->get_ustate(ust).to_string().c_str(),
          states->get_dstate(ust).to_string().c_str());
    }
  }

  void print_lattice()
  {
    printf("\nLattice: \n");
    for (unsigned int i = 0; i < lattice.size(); i++)
    {
      int me = lattice[i].me;
      int n1 = lattice[i].n1;
      printf("%10d%10d\n", me, n1);
    }
  }

  void print_matrix_element(int il1, int il2)
  {
    int nstates = states->get_nstates();
    vector<double> v1(nstates,0.0);
    v1[il1] = 1.0;
    vector<double> v2(nstates,0.0);
    v2[il2] = 1.0;
    vector<double> rvec = apply(v2);
    double me = dotvec(rvec,v1);
    printf("Matrix element for [%5d,%5d] = %7.3f\n",il1,il2,me);
  }

  vector<double> make_matrix()
  {
    return make_matrix(this->states);
  }

  vector<double> make_matrix(States* s)
  {
    int nstates = s->get_nstates();
    vector<double> rmat(nstates*nstates,0.0);

    for (int i = 0; i < nstates; i++)
    {
      vector<double> v1(nstates,0.0);
      v1[i] = 1.0;
      for (int j = 0; j < nstates; j++)
      {
        vector<double> v2(nstates,0.0);
        v2[j] = 1.0;

        // DEBUG
        //printf("V1:  ");
        //for (int iv = 0; iv < v1.size(); iv++) printf("%5.3f  ",v1[iv]);
        //printf("\n\n");
        //printf("V2:  ");
        //for (int iv = 0; iv < v2.size(); iv++) printf("%5.3f  ",v2[iv]);
        //printf("\n\n");
 
        vector<double> rvec = apply(v1);
        double t1 = dotvec(rvec,v2);
        rmat[i*nstates+j] = t1;
      }
    }

    // make sure that the matrix is symmetric
    for (int i = 0; i < nstates; i++)
    {
      for (int j = 0; j < i; j++)
      {
        printf("%15.8f%15.8f\n",rmat[i*nstates+j],rmat[j*nstates+i]);
        if (fabs(rmat[i*nstates+j]-rmat[j*nstates+i]) > tol)
        {
          printf("ERROR: (make_matrix()) elements [%3d,%3d] != [%3d,%3d]\n", i,j,j,i);
          exit(EXIT_FAILURE);
        }
      }
    }

    return rmat;
  }

  void print_apply_U(bitset<NSIZE> us, bitset<NSIZE> ds, double uval)
  {
    printf("%s    %s     U = %10.3f\n", us.to_string().c_str(),
        ds.to_string().c_str(),uval);
  }

  void print_apply_thop(bitset<NSIZE> bs1, bitset<NSIZE> bs2, int me, int n1, double tval)
  {
    printf("%s ----> %s     [%4d%4d]  t = %10.3f\n", bs1.to_string().c_str(),
        bs2.to_string().c_str(),me,n1,tval);
  }

  void create_2d_square_lattice()
  {
    // loop over lattice sites (delegate to lattice object?)
    for (int si = 0; si < nx; si++)
    {
      for (int sj = 0; sj < ny; sj++)
      {
        // get index
        int me = si*ny+sj;

        // the (0,-1) case
        int si2 = si;
        int sj2 = (sj==0) ? (ny-1) : sj-1;
        int n1 = si2*ny+sj2;
        lattice.push_back(neighbor_2d(me,n1));

        // the (0,+1) case
        si2 = si;
        sj2 = (sj==(ny-1)) ? 0 : sj+1;
        int n2 = si2*ny+sj2;
        lattice.push_back(neighbor_2d(me,n2));

        // the (-1,0) case
        si2 = (si==0) ? (nx-1) : si-1;
        sj2 = sj;
        int n3 = si2*ny+sj2;
        lattice.push_back(neighbor_2d(me,n3));

        // the (+1,0) case
        si2 = (si==(nx-1)) ? 0 : si+1;
        sj2 = sj;
        int n4 = si2*ny+sj2;
        lattice.push_back(neighbor_2d(me,n4));
      }
    }
  }

  void create_1d_lattice()
  {
    // loop over lattice sites (delegate to lattice object?)
    for (int si = 0; si < nx; si++)
    {
      // the -1 case
      int si2 = (si==0) ? (nx-1) : si-1;
      lattice.push_back(neighbor_2d(si,si2));

      // the +1 case
      si2 = (si==(nx-1)) ? 0 : si+1;
      lattice.push_back(neighbor_2d(si,si2));
    }
  }

  virtual vector<double> apply_w_openmp(const vector<double>& svec1)
  {
    vector<double> svec2(svec1.size(),0.0);

    // loop of coeffs of state
    int nsites = nx*ny;
    int ndown2 = states->get_ndown2();
    #pragma omp parallel
    {
      vector<double> _svec2(svec1.size(),0.0);
      #pragma omp for
      for (unsigned ist = 0; ist < svec1.size(); ist++)
      {
        int ust = ist / ndown2;
        int dst = ist % ndown2;
        bitset<NSIZE> us = states->get_ustate(ust);
        bitset<NSIZE> ds = states->get_dstate(dst);
        string ustr = us.to_string();
        string dstr = ds.to_string();

        double val1 = svec1[ist];
        // double occupancies
        if (fabs(val1) > tol)
        {
          double totalU = total_U_for_state(us, ds, U);
          _svec2[ust*ndown2+dst] += totalU*val1;
        }

        for (unsigned int li = 0; li < lattice.size(); li++)
        {
          neighbor_2d nghbr = lattice[li];
          int me = nghbr.me;
          int n1 = nghbr.n1;

          // up spin
          // destination state (do hopping)
          if ((fabs(val1) > tol) && !us[me] && us[n1])
          {
            bitset<NSIZE> us2 = thop(me,n1,us);
            string us2str = us2.to_string();
            // get the index of the destination state
            unsigned long du1 = states->get_ustate_index(us2);
            // sign
            double ut1 = -t[li]*computephase(me,n1,us);
            // write destination
            _svec2[du1*ndown2+dst] += ut1*val1;
          }

          // down spin
          // destination state (do hopping)
          if ((fabs(val1) > tol) && !ds[me] && ds[n1])
          {
            bitset<NSIZE> ds2 = thop(me,n1,ds);
            string ds2str = ds2.to_string();
            // get the index of the destination state
            unsigned long dd1 = states->get_dstate_index(ds2);
            // sign
            double dt1 = -t[li]*computephase(me,n1,ds);
            // write destination
            _svec2[ust*ndown2+dd1] += dt1*val1;
          }
        }
      }
      #pragma omp critical (update_sum)
      {
        for (unsigned int i = 0; i < svec2.size(); i++) svec2[i] += _svec2[i];
      }
    }
    return svec2;
  }

  virtual vector<double> apply(const vector<double>& svec1)
  {
    return apply(svec1, states);
  }

  virtual vector<double> apply(const vector<double>& svec1, States* s)
  {
    // output vector
    vector<double> svec2(svec1.size(),0.0);

    // loop of coeffs of state
    int ndown2 = s->get_ndown2();
    for (unsigned int ist = 0; ist < svec1.size(); ist++)
    {
      // get up and down spins portions
      int ust = ist / ndown2;
      int dst = ist % ndown2;
      bitset<NSIZE> us = s->get_ustate(ust);
      bitset<NSIZE> ds = s->get_dstate(dst);
      string ustr = us.to_string();
      string dstr = ds.to_string();

      double val1 = svec1[ist];
      // double occupancies
      if (fabs(val1) > tol)
      {
        double totalU = total_U_for_state(us, ds, U);
        svec2[ust*ndown2+dst] += totalU*val1;
      }

      for (unsigned int li = 0; li < lattice.size(); li++)
      {
        neighbor_2d nghbr = lattice[li];
        int me = nghbr.me;
        int n1 = nghbr.n1;

        // up spin
        // destination state (do hopping)
        if ((fabs(val1) > tol) && !us[me] && us[n1])
        {
          bitset<NSIZE> us2 = thop(me,n1,us);
          string us2str = us2.to_string();

          // get the index of the destination state
          unsigned long du1 = s->get_ustate_index(us2);
          // sign
          double ut1 = -t[li]*computephase(me,n1,us);
          // write destination
          svec2[du1*ndown2+dst] += ut1*val1;
          
          // DEBUG
          //int ist2 = du1*ndown2+dst;
          //if ((ist == 0 && ist2 == 2) || (ist == 2) && (ist2 == 0))
          //{
          //  printf("[[[ ist = %d      ist2 = %d ]]]\n",ist,ist2);
          //  printf("*****************************************************\n");
          //  cout << "ustr: " << ustr << endl;
          //  cout << "dstr: " << dstr << endl;
          //  cout << endl;

          //  printf("\nme = %d     n1 = %d\n",me,n1);
          //  cout << "ustr2: " << us2str << endl;
          //  printf("*****************************************************\n\n");
          //}

        }

        // down spin
        // destination state (do hopping)
        if ((fabs(val1) > tol) && !ds[me] && ds[n1])
        {
          bitset<NSIZE> ds2 = thop(me,n1,ds);
          string ds2str = ds2.to_string();

          // get the index of the destination state
          unsigned long dd1 = s->get_dstate_index(ds2);
          // sign
          double dt1 = -t[li]*computephase(me,n1,ds);
          // write destination
          svec2[ust*ndown2+dd1] += dt1*val1;

          // DEBUG
          //int ist2 = dd1*ndown2+dst;
          //if ((ist == 0 && ist2 == 2) || (ist == 2) && (ist2 == 0))
          //{
          //  printf("[[[ ist = %d      ist2 = %d ]]]\n",ist,ist2);
          //  printf("*****************************************************\n");
          //  cout << "ustr: " << ustr << endl;
          //  cout << "dstr: " << dstr << endl;
          //  cout << endl;

          //  printf("\nme = %d     n1 = %d\n",me,n1);
          //  cout << "dstr2: " << ds2str << endl;
          //  printf("*****************************************************\n\n");
          //}

        }
      }
    }
    return svec2;
  }

  void compute_1p_greens_function_matrix(double omega0,
                          double omega1, int nomega, double eta)
  {
    // ground state
    int lsize = 100;
    Lanczos<HubbardCalculation> l(this,lsize);
    l.run();
    vector<double> gs = l.lowstate();

    // create enlarged Hilbert space (up electrons)
    int nsites = nx*ny;
    States s(states->get_nup()+1,states->get_ndown(), nsites);

    // energy points
    vector< complex<double> > gf1p(nomega,0.0);
    double domega = (omega1-omega0)/nomega;

    // create Hamiltonian matrix (with enlarged Hilbert space)
    vector<double> H = make_matrix(&s);

    // We're going to reexpress ground state in a larger Hilbert space
    // where we've created a particle at site isite.
    for (int isite = 0; isite < 1; isite++)
    {
      vector<double> gs2(s.get_nstates(),0.0);
      int ndwn = states->get_ndown2();
      for (unsigned int ist = 0; ist < gs.size(); ist++)
      {
        // if gs[ist] is not big enough to bother
        if (fabs(gs[ist]) > tol)
        {
          int ust = ist / ndwn;
          int dst = ist % ndwn;
          // need to make copy so we use the copy constructor
          bitset<NSIZE> us(states->get_ustate(ust));
          bitset<NSIZE> ds(states->get_dstate(dst));
          // get string for debugging purposes
          string ustr = us.to_string();
          string dstr = ds.to_string();

          // Now, we'ed like to create a particle at isite
          // Is there already an electron at isite with up spin?
          if (!us[isite])
          {
            // create "up" particle at isite, and translate this
            // "new" state to gs2
            us[isite] = true;
            unsigned long ust2 = s.get_ustate_index(us);
            unsigned int ist2 = ust2*ndwn+dst;
            gs2[ist2] = gs[ist];
          }
        }
      }

      // loop over energy grid
      for (int iw = 0; iw < nomega; iw++)
      {
        // create matrix
        int nstates = s.get_nstates();
        vector< complex<double> > zmat(nstates*nstates,0.0);
        double omega = omega0 + iw * domega;
        for (int iz = 0; iz < nstates; iz++)
        {
          zmat[iz*nstates+iz] = complex<double>(omega-H[iz*nstates+iz],eta);
          for (int jz = 0; jz < iz; jz++)
          {
            zmat[iz*nstates+jz] = complex<double>(-H[iz*nstates+jz],0.0);
          }
        }

        // invert matrix
        vector< complex<double> > izmat = invert_matrix(zmat,nstates);
        // get first element
        complex<double> t1 = izmat[0];
        double t2 = 1/(double)nsites;
        gf1p[iw] += t1 * t2;
      }
    }

    for (int iw = 0; iw < nomega; iw++)
    {
      double omega = omega0 + iw * domega;
      char tmpstr[256];
      sprintf(tmpstr, "%15.8f   %15.8f   %15.8f",omega,real(gf1p[iw]), imag(gf1p[iw]));
      cout << tmpstr << std::endl;
    }

  }

  void compute_1p_greens_function(double omega0,
                          double omega1, int nomega, double eta)
  {
    // ground state
    int lsize = 100;
    Lanczos<HubbardCalculation> l(this,lsize);
    l.run();
    vector<double> gs = l.lowstate();

    // create enlarged Hilbert space (up electrons)
    int nsites = nx*ny;
    States s(states->get_nup()+1,states->get_ndown(), nsites);

    // energy points
    vector< complex<double> > gf1p(nomega,0.0);
    double domega = (omega1-omega0)/nomega;

    // We're going to reexpress ground state in a larger Hilbert space
    // where we've created a particle at site isite.
    for (int isite = 0; isite < 1; isite++)
    {
      vector<double> gs2(s.get_nstates(),0.0);
      int ndwn = states->get_ndown2();
      for (unsigned int ist = 0; ist < gs.size(); ist++)
      {
        // if gs[ist] is not big enough to bother
        if (fabs(gs[ist]) > tol)
        {
          int ust = ist / ndwn;
          int dst = ist % ndwn;
          // need to make copy so we use the copy constructor
          bitset<NSIZE> us(states->get_ustate(ust));
          bitset<NSIZE> ds(states->get_dstate(dst));
          // get string for debugging purposes
          string ustr = us.to_string();
          string dstr = ds.to_string();

          // Now, we'ed like to create a particle at isite
          // Is there already an electron at isite with up spin?
          if (!us[isite])
          {
            // create "up" particle at isite, and translate this
            // "new" state to gs2
            us[isite] = true;
            unsigned long ust2 = s.get_ustate_index(us);
            unsigned int ist2 = ust2*ndwn+dst;
            gs2[ist2] = gs[ist];
          }
        }
      }

      // swap enlarged Hilbert space and smaller
      States* ststmp = states;
      states = &s;

      // we want to construct an approximation for H (i.e. the diagonal
      // and off diagonal coefficients) in the Lanczos basis, but we
      // want to apply to our N+1 ground state or other
      // words c_i^{+} | GS >
      Lanczos<HubbardCalculation> l2(this,lsize,gs2);
      l2.run();

      // loop over energy grid
      for (int iw = 0; iw < nomega; iw++)
      {
        // create matrix
        vector< complex<double> > zmat(lsize*lsize,0.0);
        double omega = omega0 + iw * domega;
        for (int iz = 0; iz < lsize; iz++)
        {
          zmat[iz*lsize+iz] = complex<double>(omega-l._a[iz],eta);
          if (iz < lsize-1)
          {
            zmat[(iz+1)*lsize+iz] = complex<double>(-l._b[iz],0.0);
            zmat[iz*lsize+iz+1] = complex<double>(-l._b[iz],0.0);
          }
        }
        // invert matrix
        vector< complex<double> > izmat = invert_matrix(zmat,lsize);
        // get first element
        complex<double> t1 = izmat[0];
        double t2 = 1/(double)nsites;
        gf1p[iw] += t1 * t2;
      }

      // return states
      states = ststmp;
    }

    std::ofstream fstr("gf1p.txt");
    for (int iw = 0; iw < nomega; iw++)
    {
      double omega = omega0 + iw * domega;
      char tmpstr[256];
      sprintf(tmpstr, "%15.8f   %15.8f   %15.8f",omega,real(gf1p[iw]), imag(gf1p[iw]));
      cout << tmpstr << std::endl;
    }
    fstr.close();

  }

  bitset<NSIZE> thop(const int& s1, const int& s2, const bitset<NSIZE>& state)
  {
    // tij --> i == s1 and j == s2
    if (!state[s1] && state[s2])
    {
      bitset<NSIZE> state2(state);
      state2.flip(s1);
      state2.flip(s2);
      return state2;
    }
  }

  double computephase(int s1, int s2, bitset<NSIZE> state)
  {
    if (s1 == s2) return 0.0;

    int sum = 0;
    if (s1 > s2)
    {
      for (int i = s2+1;i < s1; i++)
      {
        if (state[i]) sum++;
      }
    }
    else if (s1 < s2)
    {
      for (int i = s1+1;i < s2; i++)
      {
        if (state[i]) sum++;
      }
    }
    return ((sum % 2) == 0) ? 1.0 : -1.0;

  }

//  vector<double> test_apply(const vector<double>& v)
//  {
//    int n = v.size();
//    vector<double> A(n*n,0.0);
//    vector<double> v2(n,0.0);
////    for (int j = 0; j < n; j++)
////    {
////      for (int i = 0; i < j+1; i++)
////      {
////        A[i*n+j] = i*n+j;
////        A[j*n+i] = i*n+j;
////      }
////    }
//    A = make_matrix();
//    for (int i = 0; i < n; i++)
//    {
//      for (int j = 0; j < n; j++)
//      {
//        v2[i] += A[i*n+j] * v[j];
//      }
//    }
//
//    return v2;
//  }

  void lanczos(int n = 100)
  {
    // create random initial vector
    int nstates = states->get_nstates();
    vector<double> v = random_vector(nstates);

//    // save initial vector
//    for (int i = 0; i < nstates; i++)
//    {
//      vinit[i] = v[i];
//    }

    vector<double> vold(nstates,0.0);
    vector<double> a(n,0.0);
    vector<double> b(n-1,0.0);
    vector<double> v2(nstates,0.0);
    for (int i = 0; i < n; i++)
    {
      if (i > 0)
      {
        v2 = wst_daxpy(1.0,apply(v),-b[i-1],vold);
      }
      else
      {
        v2 = apply(v);
      }
      a[i] = dotvec(v,v2);
      if (i < (n-1))
      {
        wst_daxpy_inplace(1.0,v2,-a[i],v);
        b[i] = norm2(v2);
        for (int iv = 0; iv < nstates; iv++)
        {
          vold[iv] = v[iv];
          v[iv] = v2[iv]/b[i];
        }
      }
    }

    vector<double> mat(n*n,0.0);
    for (int i = 0; i < n-1; i++)
    {
      mat[i*n+i] = a[i];
      mat[i*n+i+1] = b[i];
      mat[(i+1)*n+i] = b[i];
    }
    mat[n*n-1] = a[n-1];

    vector<double> e(n,0.0);
    vector<double> ev(n*n, 0.0);
    diag_matrix(mat,n,e,ev);
    printf("Lanczos: lowest eigenvalue is %15.8f\n\n", e[0]);
  }
};
//*****************************************************************************



#endif /* HUBBARD_H_ */
