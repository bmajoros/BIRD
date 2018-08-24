/****************************************************************
 Replicate.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Replicate.H"
using namespace std;


Replicate::Replicate()
  : ref(0), alt(0)
{
  // ctor
}


Replicate::Replicate(int ref,int alt)
  : ref(ref), alt(alt)
{
}



int Replicate::getRef() const 
{ 
  return ref; 
}



int Replicate::getAlt() const 
{ 
  return alt; 
}



void Replicate::setRef(int r) 
{ 
  ref=r; 
}



void Replicate::setAlt(int a) 
{ 
  alt=a; 
}



void Replicate::add(const Replicate &other)
{
  ref+=other.ref;
  alt+=other.alt;
}


void Replicate::addPseudocount(int pseudo)
{
  ref+=pseudo;
  alt+=pseudo;
}





