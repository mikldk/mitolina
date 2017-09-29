#ifndef MITOLINE_TYPES_H
#define MITOLINE_TYPES_H

#define CHECK_ABORT_EVERY 100000


//#define RCPP_XPTR_2ND_ARG true // ensures that finaliser is called
#define RCPP_XPTR_2ND_ARG false // do not call finalisers, I guess we live with some memory leaks for now...!!! FIXME

class Individual;
class Pedigree;
class Population;

//class SimulateChooseMother;
//class WFRandomMother;
//class GammaVarianceRandomMother;

#include "helper_Individual.hpp"

#include "class_Individual.hpp"
#include "class_Pedigree.hpp"
#include "class_Population.hpp"
#include "class_SimulateChooseMother.hpp"

#endif
