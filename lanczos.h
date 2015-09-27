/*
 * lanczos.h
 *
 *  Created on: Aug 22, 2011
 *      Author: wsttiger
 */

#ifndef LANCZOS_H_
#define LANCZOS_H_

#include <iostream>
#include <fstream>
#include <bitset>
#include <omp.h>
#include <ctime>

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

static double ttt, sss;
void START_TIMER() {
//    ttt=wall_time(); sss=cpu_time();
}

void END_TIMER(const char* msg) {
//    ttt=wall_time()-ttt; sss=cpu_time()-sss; printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}

//*****************************************************************************
struct MatrixLanczos
{
  vector<double> _H;

  unsigned int _nstates; // is the bigger dimension (full size)
  unsigned int _nsize; // is the size of the approximate matrix in the Lanczos basis

  vector<double> _vinit;
  vector<double> _a;
  vector<double> _b;
  vector<double> _e;
  vector<double> _ev;

  MatrixLanczos(const vector<double>& H, const int& nsize)
  {
    _H = vector<double>(H);
    _nstates = int(floor(sqrt(H.size())));
    printf("MatrixLanczos constructor: H.size() = %5d    _nstates = %5d\n\n", int(H.size()), _nstates);
    _nsize = nsize;
    create_init_random_vector();
  }

  MatrixLanczos(const vector<double>& H, const int& nsize, const vector<double>& vinit)
  {
    _H = vector<double>(H);
    _nstates = int(floor(sqrt(H.size())));
    printf("MatrixLanczos constructor: H.size() = %5d    _nstates = %5d\n\n", int(H.size()), _nstates);
    _nsize = nsize;
    assert(vinit.size() == _nstates);
    _vinit = vector<double>(vinit);
  }

  void create_init_random_vector()
  {
    _vinit = random_vector(_nstates);
  }

  // run() computes the _a and _b parameters (matrix diagonal
  // and off-diagonal) in the Lanczos basis from an initial vector
  // Initial vector can be generated or provided
  void run()
  {
    vector<double> v(_nstates,0.0);
    for (unsigned int i = 0; i < _nstates; i++)
    {
      v[i] = _vinit[i];
    }
    vector<double> vold(_nstates,0.0);
    _a = vector<double>(_nsize,0.0);
    _b = vector<double>(_nsize-1,0.0);
    vector<double> v2(_nstates,0.0);
    for (unsigned int i = 0; i < _nsize; i++)
    {
      if (i > 0)
      {
        v2 = wst_daxpy(1.0,matrix_vector_mult(_H,v),-_b[i-1],vold);
      }
      else
      {
        v2 = matrix_vector_mult(_H,v);
      }
      _a[i] = dotvec(v,v2);
      if (i < (_nsize-1))
      {
        wst_daxpy_inplace(1.0,v2,-_a[i],v);
        _b[i] = norm2(v2);
        for (unsigned int iv = 0; iv < _nstates; iv++)
        {
          vold[iv] = v[iv];
          v[iv] = v2[iv]/_b[i];
        }
      }
    }

    vector<double> mat(_nsize*_nsize,0.0);
    for (unsigned int i = 0; i < _nsize-1; i++)
    {
      mat[i*_nsize+i] = _a[i];
      mat[i*_nsize+i+1] = _b[i];
      mat[(i+1)*_nsize+i] = _b[i];
    }
    mat[_nsize*_nsize-1] = _a[_nsize-1];

    _e = vector<double>(_nsize,0.0);
    _ev = vector<double>(_nsize*_nsize, 0.0);
    diag_matrix(mat,_nsize,_e,_ev);
    printf("Lanczos: lowest eigenvalue is %15.8f\n\n", _e[0]);
  }

};
//*****************************************************************************

//*****************************************************************************
template <typename T>
struct Lanczos
{
  T* _model;
  unsigned int _nstates; // is the bigger dimension (full size)
  unsigned int _nsize; // is the size of the approximate matrix in the Lanczos basis

  vector<double> _vinit;
  vector<double> _a;
  vector<double> _b;
  vector<double> _e;
  vector<double> _ev;



  Lanczos(T* model, int nsize)
  {
    _model = model;
    _nstates = _model->get_nstates();
    _nsize = nsize;
    create_init_random_vector();
  }

  Lanczos(T* model, int nsize, const vector<double>& vinit)
  {
    _model = model;
    _nstates = _model->get_nstates();
    _nsize = nsize;
    assert(vinit.size() == _nstates);
    _vinit = vector<double>(_nstates,0.0);
    for (unsigned int i = 0; i < _nstates; i++)
    {
      _vinit[i] = vinit[i];
    }
  }

  void create_init_random_vector()
  {
    _vinit = vector<double>(_nstates,0.0);
    for (unsigned int i = 0; i < _nstates; i++)
    {
      int i1 = rand();
      double t1 = (i1 % 100000)/100000.0;
      _vinit[i] = t1;
    }
    normalize(_vinit);
  }

  // run() computes the _a and _b parameters (matrix diagonal
  // and off-diagonal) in the Lanczos basis from an initial vector
  // Initial vector can be generated or provided
  void run()
  {
    vector<double> v(_nstates,0.0);
    for (unsigned int i = 0; i < _nstates; i++)
    {
      v[i] = _vinit[i];
    }
    vector<double> vold(_nstates,0.0);
    _a = vector<double>(_nsize,0.0);
    _b = vector<double>(_nsize-1,0.0);
    vector<double> v2(_nstates,0.0);
    for (unsigned int i = 0; i < _nsize; i++)
    {
      START_TIMER();
      if (i > 0)
      {
        v2 = wst_daxpy(1.0,_model->apply(v),-_b[i-1],vold);
      }
      else
      {
        v2 = _model->apply(v);
      }
      END_TIMER("lanczos run 1");
      START_TIMER();
      _a[i] = dotvec(v,v2);
      END_TIMER("lanczos run 2");
      if (i < (_nsize-1))
      START_TIMER();
      {
        wst_daxpy_inplace(1.0,v2,-_a[i],v);
        _b[i] = norm2(v2);
        for (unsigned int iv = 0; iv < _nstates; iv++)
        {
          vold[iv] = v[iv];
          v[iv] = v2[iv]/_b[i];
        }
      }
      END_TIMER("lanczos run 3");
    }

    vector<double> mat(_nsize*_nsize,0.0);
    for (unsigned int i = 0; i < _nsize-1; i++)
    {
      mat[i*_nsize+i] = _a[i];
      mat[i*_nsize+i+1] = _b[i];
      mat[(i+1)*_nsize+i] = _b[i];
    }
    mat[_nsize*_nsize-1] = _a[_nsize-1];

    _e = vector<double>(_nsize,0.0);
    _ev = vector<double>(_nsize*_nsize, 0.0);
    diag_matrix(mat,_nsize,_e,_ev);
    printf("Lanczos: lowest eigenvalue is %15.8f\n\n", _e[0]);
  }

  // run(vector) computes the lowest eigenstate using the _a and _b
  // parameters that were already generated from a previous run
  vector<double> run(const vector<double>& vin)
  {
    // create return vector
    vector<double> rv(_nstates,0.0);
    // copy initial vector
    vector<double> v(_nstates,0.0);
    for (unsigned int i = 0; i < _nstates; i++)
    {
      v[i] = vin[i];
    }
    vector<double> vold(_nstates,0.0);
    vector<double> v2(_nstates,0.0);
    for (unsigned int i = 0; i < _nsize; i++)
    {
      // build up return value
      //wst_daxpy_inplace(1.0,rv,_ev[i*_nsize],v);
      for (unsigned int iv = 0; iv < _nstates; iv++)
      {
        rv[iv] += _ev[i*_nsize]*v[iv];
      }
      if (i > 0)
      {
        v2 = wst_daxpy(1.0,_model->apply(v),-_b[i-1],vold);
      }
      else
      {
        v2 = _model->apply(v);
      }
      if (i < (_nsize-1))
      {
        wst_daxpy_inplace(1.0,v2,-_a[i],v);
        for (unsigned int iv = 0; iv < _nstates; iv++)
        {
          vold[iv] = v[iv];
          v[iv] = v2[iv]/_b[i];
        }
      }
    }
    normalize(rv);
    return rv;
  }

  vector<double> lowstate()
  {
    return run(_vinit);
  }

};
//*****************************************************************************

#endif
