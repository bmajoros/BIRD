/****************************************************************
 Replicates.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Replicates_H
#define INCL_Replicates_H
#include <iostream>
#include "BOOM/Vector.H"
#include "Replicate.H"
using namespace std;
using namespace BOOM;


class Replicates {
  Vector<Replicate> reps;
public:
  Replicates();
  void add(const Replicate &);
  int size();
  Replicate &operator[](int i);
  void collapse();
  void clear();
};

#endif

