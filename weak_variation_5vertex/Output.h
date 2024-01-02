#ifndef OUTPUT_H
#define OUTPUT_H

#include"World.h"
#include<vector>
#include<fstream>
#include"Species.h"
#include"Source.h"
#include"Crosssection.h"
#include "Potentialsolver.h"
namespace Output 
{
    void fields(World &world, std::vector<Species> &species);
    void screenOutput(World &world, std::vector<Species> &species);
    void paramOutput(World &world, std::vector<Species> &species, std::vector<ColdBeamSource> &coldbeamsource, CrossSection &cross_section, Potentialsolver &solver);
    void particles (World &world, std::vector<Species> &species, int num_parts);
    void particles2 (World &world, Species &sp);
    void diagOutput(World &world, std::vector<Species> &species, int langkah_awal); //kita pakai diagoutput juga
    void flagOutput(World &world);
}

#endif
