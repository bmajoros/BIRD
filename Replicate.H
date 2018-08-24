/****************************************************************
 Replicate.H
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Replicate_H
#define INCL_Replicate_H
#include <iostream>
using namespace std;


class Replicate {
  int ref, alt; // Read counts for reference and alternate alleles
public:
  Replicate();
  Replicate(int ref,int alt);
  int getRef() const;
  int getAlt() const;
  void setRef(int r);
  void setAlt(int a);
  void add(const Replicate &);
  void addPseudocount(int);
};

#endif

