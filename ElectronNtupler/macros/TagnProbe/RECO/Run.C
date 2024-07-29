#include "TString.h"
#include "TROOT.h"

void Run(){
gROOT->ProcessLine(".L LeptonScaleLookup1.cc");
gROOT->ProcessLine(".L FinalEffSingle.C");
gROOT->ProcessLine("FinalEffSingle()");
}
