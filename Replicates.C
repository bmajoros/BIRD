/****************************************************************
 Replicates.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Replicates.H"
using namespace std;
using namespace BOOM;

Replicates::Replicates()
{
  // ctor
}



void Replicates::add(const Replicate &r)
{
  reps.push_back(r);
}



Replicate &Replicates::operator[](int i)
{
  return reps[i];
}



void Replicates::collapse()
{
  Replicate r;
  const int n=reps.size();
  for(int i=0 ; i<n ; ++i) r.add(reps[i]);
  reps.resize(1);
  reps[0]=r;
}



void Replicates::clear()
{
  reps.purge();
}



int Replicates::size()
{
  return reps.size();
}



