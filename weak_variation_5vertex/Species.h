#ifndef SPECIES_H
#define SPECIES_H

#include<vector>
#include"Field.h"
#include"World.h"
#include<math.h>
#include<fstream>
#include<sstream>
#include<iostream>
#include<iomanip> 
#include"Crosssection.h"
//container untuk partikel (3 posisi, 3 kecepatan, 1 bobot makropartikel)
//meskipun mesh 2D
struct Particle
{
    //constructornya di sini
    Particle(double3 pos, double3 vel, double mpw) : pos{pos}, vel{vel}, mpw{mpw} { }
    double3 pos;
    double3 vel;
    double mpw;
};

//sementara tanpa trail dulu (trail 2 dimensi, agak hati-hati)

class Species
{
    public:
        //constructor parametrik dari spesies
        Species(std::string name, double mass, double charge, double mpw0, World &world, const int langkah_energi, CrossSection &cross_section, double Bz) : 
            name{name}, mass{mass}, charge{charge}, mpw0{mpw0}, den(world.nz,world.nr), world{world}, energy{langkah_energi}, eNonAxial{langkah_energi}, cross_section{cross_section}, magnetz{Bz} {
            }
            
        DField den;

        std::vector<Particle> particles;

        double getRealCount();

        void advance(); //untuk benchmark, sementara tanpa magnet. Jadi cukup leapfrog biasa
		
		void advanceM(double dt);

        void advanceBoris2(double dt); //boris solver untuk ion positif

        void advanceBorisElektron(double dt);

        void borisPush(double3 &pos, double3 &vel, double dt, double t, double s);

        void negativeRotation(double3& pos, double3 &vel);

        void positiveRotation(double3& pos, double3 &vel);

        bool boundaryCheck(double2 posawal, double2 &posakhir, double &dt, double3 &vel, double &mpw, int& dinding);

        void nullCollision(double3 &vel, double& mpw, double3 pos, int &elastic_count);

        void advanceBoris3(double dt); //boris solver untuk ion negatif

        void advanceBorisElektron2(double dt);

        void advanceBorisIon(double dt);

        void advanceCyclotronicDrift(double3 &pos, double3 &vel, double dt);

        void advanceCyclotronicKick(double3 &pos, double3 &vel, double dt);

        void advanceElectronCyclotronic(double dt);

        void advanceIonCyclotronic(double dt);

        void advanceCyclotronicDriftCartesian(double omega, double3 &pos, double3 &vel, double dt);

        void advanceCyclotronicMod(double dt);
        
        void advanceCyclotronicModIon(double dt);

        //tiap spesies memiliki nilai repat bilangan masing-masing, rapat muatan nanti tinggal dikalikan dengan muatan masing spesies, kemudian dijumlah
        void computeNumberDensity(); 

        void addParticle(double3 pos, double3 vel, double mpw); //kecepatan mundur setengah langkah (leapfrog)

        void addParticle2(double3 pos, double3 vel, double mpw); //mundur setengah langkah untuk boris silinder

        void addParticle3(double3 pos, double3 vel, double mpw); //apa adanya, untuk cyclotronic integrator [151022]

        void addParticle4(double3 pos, double3 vel, double mpw); //untuk cyclotronic di koordinat silinder dengan rotasi koordinat [121122]

        void momTransfer(double3 &vel, double vv, double Ek); //tumbukan elastik transfer momentum

        void elasticRot(double3 &vel); //tumbukan elastik transfer momentum, ini sama aja kayak momTransfer mestinya ya?

        void excitation(double3 &vel, double Eexc, double vv, double Ek);

        void ionization(double3 &vel, double3 pos, double mpw, double Eion, double vv, double Ek); //ionisasi, sejauh ini untuk heatsink aja

        double calcThetaSamping(double3 vel);

        double calcThetaBawah(double3 vel);

        double calcDelta(double delta_max, double E, double Eth, double Emax, double theta);

        void secondEmissionSamping(double delta_e, double3 pos, double3 vel, double vv, double temp);

        void secondEmissionCathode(double delta_e, double3 pos, double3 vel, double vv, double temp);

        void secondEmissionBawah(double delta_e, double3 pos, double3 vel, double vv, double temp);

        //ini nanti loading lewat source.h, jadi nggak perlu pakai loadbox.
        
        //tetap definisikan fungsi gamma, supaya nantinya praktis kalau code dijadikan relativistik
        double gamma(double3 vel)
        {
            double v2 = (vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2])/(Const::SPEEDC*Const::SPEEDC);
            return 1.0/sqrt(1-v2);
        }
        
        //ini untuk plot emitansinya
        void emittance()
        {
            std::ofstream fout("results/emittance.txt");
            for (Particle &part:particles) //ini lebih gampang dipanggil kalo dibikin di classnya yang terkait langsung
            {
                double gamm = gamma(part.vel);
                double3 mom = (gamm/Const::SPEEDC)*part.vel;
                fout << part.pos[1] << "\t" << mom[1] << "\t"  << "\n";
            }   
            fout.close();
        }
        
        //ini untuk menghitung rms emittance
        void rms_emittance()
        {
            double sumx = 0; 
            double sumpx = 0; 
            double sumxpx = 0; 
            double emitansi_rmsx; 
            for (Particle &part:particles)
            {
                double gamm = gamma(part.vel);
                sumx += part.pos[1]*part.pos[1];
                
                sumpx += gamm*gamm*part.vel[1]*part.vel[1]/pow(Const::SPEEDC,2);
                
                sumxpx += part.pos[1]*gamm*part.vel[1]/Const::SPEEDC;
                
            }
            double avgx = sumx/getNp(); 
            double avgpx = sumpx/getNp(); 
            double avgxpx = pow(sumxpx/getNp(),2); 
            emitansi_rmsx = sqrt(avgx*avgpx - avgxpx);
            std::ofstream fout("results/emittance_rms.txt");
            fout << "emitansi arah x = " << emitansi_rmsx << "\n";
            fout.close();
        }
        const std::string name;
        const double mass;
        const double charge;
        const double mpw0;
        
        //untuk plotting energi partikel, kita butuh semacam XtoL
        double EtoL(double E) const
        {
            return E/stepEnergi;
        }

        Energy energy; //jadi ada 10000 langkah ya, sementara coba untuk sampai maksimal 1keV, mestinya masih aman

        Energy eNonAxial; //untuk kasus kecepatan non-aksial (akan kita bandingkan)

        size_t getNp()
        {
            return particles.size();
        }

        //untuk mengisi vektor energi dengan menggunakan scatter
        void scatter_energi()
        {
            for (Particle &part:particles)
            {
                double E = 0.5*mass*(part.vel[0]*part.vel[0] + part.vel[1]*part.vel[1] + part.vel[2]*part.vel[2])/abs(charge);
                double Enon = 0.5*mass*(part.vel[1]*part.vel[1]+part.vel[2]*part.vel[2])/abs(charge);
                //kalau melebihi batas, jangan discatter
                if (E >= batasAtasEnergi )
                {
                    continue;
                }
                else
                {
                    double le = EtoL(E);
                    double nle = EtoL(Enon);
                    energy.scatter(le,part.mpw); //pakai bobot partikel masing-masing, bisa jadi nanti pakai bobot beda-beda untuk satu spesies kan
                    eNonAxial.scatter(nle,part.mpw); //skema scatternya sama untuk kecepatan non-aksial
                } 
            }
        }

        //void scatter_energi untuk kecepatan non axial.

        //untuk plot energi kinetik partikel, sementara masih pakai interpolasi orde 1. Harusnya pakai cubic spline.
        void plot_energi()
        {
            //pertama set nilai dari Energi dulu dengan menggunakan scatter
            std::ofstream fout("plot_energi.dat");
            //std::cout << std::fixed;
            //std::cout << std::setprecision(4);
            double de = 0.1*stepEnergi;
            int lank = 10*(energy.ni-1);
            double deawal = 0;
            for (int i=0; i<lank; i++)
            {
                double elogic = EtoL(deawal);
                //std::cout << i << "\t" << deawal << "\t" <<elogic << "\t" << energy.gather(elogic) << "\t" << de << "\n";
                deawal += de;
                fout << deawal << "\t" << energy.gather(elogic) << "\n"; 
            }
            fout.close();
        }

        //plot posisi dan kecepatan partikel. 
        void ruang_fase()
        {
            //kita pakai stringstream supaya bisa langsung pakai nama partikel
            //perlu ditambahkan juga ini langkah yang keberapa, lebar waktunya berapa, bobot partikel berapa
            std::stringstream namafile;
            namafile << "ruang_fase/ruang_fase_" << name << "_" << world.getTs()  <<".dat";
            std::ofstream fout(namafile.str());
            for (Particle &part:particles)
            {
                //ini urutannya z, r, vz, vr, vtheta
                fout << part.pos[0] << "\t" << part.pos[1] << "\t" << part.vel[0] << "\t" << part.vel[1] << "\t" << part.vel[2] << "\n";
            }
            fout.close();
        }

        //untuk mengeset stepEnergi
        void setStepEnergi(double de)
        {
            stepEnergi = de;
        }
        //untuk mengambil nilai step energi
        double getStepEnergi()
        {
            return stepEnergi;
        }
        void setBatasEnergi()
        {
            batasAtasEnergi = stepEnergi*(energy.ni-1);
        }
        void cleanupEnergy()
        {

        }

        //yang diatas mending dijadikan satu aja ya [200922]
        void distribusiEnergi(double de)
        {
            //jangan lupa kosongin dulu energy nya(enonaxial nya juga) [181022]
            energy = 0;
            eNonAxial = 0;
            setStepEnergi(de);
            setBatasEnergi();
            scatter_energi();
            //bagian plot energinya kita bikin di sini lagi aja, pakai stringstream
            std::stringstream nama_plot;
            nama_plot << "results/Plot_energi_" << name <<"_"<< world.getTs() << ".dat";
            std::ofstream fout(nama_plot.str());
            double dee = 0.1*stepEnergi;
            int lank = 10*(energy.ni-1);
            double deawal = 0;

            for (int i=0; i<lank; i++)
            {
                double elogic = EtoL(deawal);
                fout << deawal << "\t" << energy.gather(elogic) << "\t" << eNonAxial.gather(elogic) << "\n";
                deawal += dee;     
            }
            fout.close();
        }

    protected:
    World &world;
    CrossSection &cross_section;
    double magnetz;
    //ukuran step energi yang akan diplot
    double stepEnergi=1.0; //dalam eV.
    //batas atas energi (dalam eV)
    double batasAtasEnergi=1.0;
};



#endif
