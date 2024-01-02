#include"Source.h"
#include<iostream>

//sementara, untuk sample ini, sumber di seluruh permukaan tutup, bisa diatur nanti misal mau lebih kecil
void ColdBeamSource::sample()
{
    double2 dh = world.getDh();
    double2 x0 = world.getX0();
    double2 xm = world.getXm();

    double Lr = xm[1] - x0[1];//ini jari-jari tabung (atau cincin, di kasus kita tabung sih)
    double A = Const::PI*(xm[1]*xm[1]-x0[1]*x0[1]);//luas area cincin (atau tabung kalau x0[1]=0)

    //selanjutnya, jumlah partikel yang disimulasikan berdasar arus, kecepatan, permukaan, dan bobot makropartikel
    double num_real = arus*world.getDt()/abs(spec.charge);
	
	std::cout << "num_real = " << num_real << std::endl;

    //untuk jumlah partikel yang disimulasikan, digunakan ini
    int num_sim = (int)(num_real/spec.mpw0 + rnd());
	
	std::cout << "num_sim = " << num_sim << std::endl;
	
    for (int i=0; i<num_sim; i++)
    {
        double3 pos {x0[0]+world.getDh()[0], sqrt(Lr*Lr*rnd() + x0[1]*x0[1]), 0};
        double3 vel {v_drift, 0, 0}; //urutannya z, r, theta
        spec.addParticle(pos,vel,spec.mpw0); //agak beda dengan Brieda
    }
}

void ColdBeamSource::sample2(int num_sim)
{
    double2 dh = world.getDh();
    double2 x0 = world.getX0();
    double2 xm = world.getXm();
    double r0 = x0[1] +dh[1] ;

    double Lr = (xm[1] - x0[1])-2*dh[1];//ini jari-jari tabung (atau cincin, di kasus kita tabung sih)
    double A = Const::PI*(xm[1]*xm[1]-x0[1]*x0[1]);//luas area cincin (atau tabung kalau x0[1]=0)

    /*
    //selanjutnya, jumlah partikel yang disimulasikan berdasar arus, kecepatan, permukaan, dan bobot makropartikel
    double num_real = arus*world.getDt()/abs(spec.charge);
	
	std::cout << "num_real = " << num_real << std::endl;

    //untuk jumlah partikel yang disimulasikan, digunakan ini
    int num_sim = (int)(num_real/spec.mpw0 + rnd());
	
	std::cout << "num_sim = " << num_sim << std::endl;
	*/
    for (int i=0; i<num_sim; i++)
    {
        double3 pos {(xm[0]-4*dh[0])*rnd() + 2*dh[0], sqrt(Lr*Lr*rnd() + r0*r0), 0};
        double3 vel {0, 0, 0}; //urutannya z, r, theta
        spec.addParticle2(pos,vel,spec.mpw0); //agak beda dengan Brieda
    }
}

void ColdBeamSource::sampleThermalBoris(int num_sim, double Ek)
{
    double2 dh = world.getDh();
    double2 x0 = world.getX0();
    double2 xm = world.getXm();
    double r0 = x0[1] + 0.1*dh[1] ;

    double Lr = (xm[1] - x0[1]-3*dh[1]);//ini jari-jari tabung (atau cincin, di kasus kita tabung sih)
    double A = Const::PI*(xm[1]*xm[1]-x0[1]*x0[1]);//luas area cincin (atau tabung kalau x0[1]=0)

    double vth = sqrt(2*Ek*abs(Const::QE)/spec.mass);
    for (int i=0; i<num_sim; i++)
    {
        double3 pos {(xm[0]-3*dh[0])*rnd() , sqrt(Lr*Lr*rnd() + r0*r0), 0};
        //kecepatannya diundi terlebih dahulu
        
        double3 vel {(vth*(rnd()+rnd()+rnd()-1.5)), (vth*(rnd()+rnd()+rnd()-1.5)), (vth*(rnd()+rnd()+rnd()-1.5))}; //ini x, y. z
        spec.addParticle2(pos,vel,spec.mpw0); //agak beda dengan Brieda
    }
}

void ColdBeamSource::sampleThermalCyclotronic(int num_sim, double Ek)
{
    double2 dh = world.getDh();
    double2 x0 = world.getX0();
    double2 xm = world.getXm();
    double r0 = x0[1] + 0.1*dh[1] ;

    double Lr = (xm[1] - x0[1]-3*dh[1]);//ini jari-jari tabung (atau cincin, di kasus kita tabung sih)
    double A = Const::PI*(xm[1]*xm[1]-x0[1]*x0[1]);//luas area cincin (atau tabung kalau x0[1]=0)

    double vth = sqrt(2*Ek*abs(Const::QE)/spec.mass);
    for (int i=0; i<num_sim; i++)
    {
        double3 pos {(xm[0]-3*dh[0])*rnd() , sqrt(Lr*Lr*rnd() + r0*r0), 0};
        //kecepatannya diundi terlebih dahulu
        
        double3 vel {(vth*(rnd()+rnd()+rnd()-1.5)), (vth*(rnd()+rnd()+rnd()-1.5)), (vth*(rnd()+rnd()+rnd()-1.5))}; //ini x, y. z
        //spec.addParticle2(pos,vel,spec.mpw0); //agak beda dengan Brieda
        spec.addParticle3(pos,vel,spec.mpw0); //agak beda dengan Brieda
    }
}

/* kalau kasusnya axissymmetric, rasanya kurang pas mendefinisikan langkah yang konstan
void ColdBeamSource::sample2(double lank)
{
    double2 x0 = world.getX0();
    double2 xm = world.getXm();

    int nz = (int)((xm[0]-x0[0])/lank);
    int nr = (int)((xm[1]-x0[1])/lank);

    //kecepatan awal partikel (semua terkait partikel dikerjakan pada koordinat kartesian)
    double3 vel = {0,0,0};
    //jangan lupa pembobotannya
    for(int k=0; k<nz; k++)
    {
        for (int r=0; r<nr; r++)
        {
            double r = x0[1]+r*lank;
            double2 pos1 = {x0[0]+k*lank, x0[1]+r*lank};
        }
    }
}
*/

void ColdBeamSource::sampleIonHPlus()
{
    double3 vel1 = {0,0,0}; //ion bergerak sangat lambat
    size_t np = world.ionsHPlus.size();
    {
        for (int p = 0; p<np; p++)
        {
            spec.addParticle2(world.ionsHPlus[p].pos,vel1,spec.mpw0);
            //spec.addParticle3(world.ionsHPlus[p].pos,vel1,spec.mpw0);
        }
    }
}

void ColdBeamSource::sampleIonH2Plus()
{
    double3 vel1 = {0,0,0}; //ion bergerak sangat lambat
    size_t np = world.ionsH2Plus.size();
    {
        for (int p = 0; p<np; p++)
        {
            spec.addParticle2(world.ionsH2Plus[p].pos,vel1,spec.mpw0);
            //spec.addParticle3(world.ionsH2Plus[p].pos,vel1,spec.mpw0);
        }
    }
}

void ColdBeamSource::sampleIonHMin()
{
    double3 vel1 = {0,0,0}; //ion bergerak sangat lambat
    size_t np = world.ionsHMin.size();
    {
        for (int p = 0; p<np; p++)
        {
            spec.addParticle2(world.ionsHMin[p].pos,vel1,spec.mpw0);
            //spec.addParticle3(world.ionsHMin[p].pos,vel1,spec.mpw0);
        }
    }
}

void ColdBeamSource::sampleElecIonCyc() 
{
    size_t np = world.elec_ionization.size();
    {
        for (int p=0; p<np; p++)
        {
            //spec.addParticle2(world.elec_ionization[p].pos, world.elec_ionization[p].vel, spec.mpw0);
            spec.addParticle3(world.elec_ionization[p].pos, world.elec_ionization[p].vel, spec.mpw0); //addparticle3 untuk integrator cyclotronic [151022]
        }
    }
}

void ColdBeamSource::sampleElecIonBoris() 
{
    size_t np = world.elec_ionization.size();
    {
        for (int p=0; p<np; p++)
        {
            spec.addParticle2(world.elec_ionization[p].pos, world.elec_ionization[p].vel, spec.mpw0);
            //spec.addParticle3(world.elec_ionization[p].pos, world.elec_ionization[p].vel, spec.mpw0); //addparticle3 untuk integrator cyclotronic [151022]
        }
    }
}

void ColdBeamSource::sampleElec(double Ek)
{
    double vth = sqrt(2*Ek*abs(Const::QE)/spec.mass);
    
    size_t np = world.elec_secondary.size();
    {
        for (int p = 0; p<np; p++)
        {
            //mestinya mengarah ke atas atau ke bawah
            double3 vel {-abs(vth*(rnd()+rnd()+rnd()-1.5)), (vth*(rnd()+rnd()+rnd()-1.5)), (vth*(rnd()+rnd()+rnd()-1.5))};
            spec.addParticle2(world.elec_secondary[p].pos,vel,spec.mpw0);
        }
    }
}

void ColdBeamSource::sampleElec2(double Ek)
{
    double vth = sqrt(2*Ek*abs(Const::QE)/spec.mass);
    
    size_t np = world.elec_secondary.size();
    {
        for (int p = 0; p<np; p++)
        {
            //mestinya mengarah ke atas atau ke bawah
            double3 vel {abs(vth*(rnd()+rnd()+rnd()-1.5)), (vth*(rnd()+rnd()+rnd()-1.5)), (vth*(rnd()+rnd()+rnd()-1.5))};
            spec.addParticle2(world.elec_secondary2[p].pos,vel,spec.mpw0);
        }
    }
}
//suhu dinding di world

void ColdBeamSource::sampleElecWallCyc()
{
    size_t np = world.elec_wall.size();
    for(int p=0; p<np; p++)
    {
        //spec.addParticle2(world.elec_wall[p].pos, world.elec_wall[p].vel,spec.mpw0);
        spec.addParticle3(world.elec_wall[p].pos, world.elec_wall[p].vel,spec.mpw0);
    }
}

void ColdBeamSource::sampleElecWallBoris()
{
    size_t np = world.elec_wall.size();
    for(int p=0; p<np; p++)
    {
        spec.addParticle2(world.elec_wall[p].pos, world.elec_wall[p].vel,spec.mpw0);
        //spec.addParticle3(world.elec_wall[p].pos, world.elec_wall[p].vel,spec.mpw0);
    }
}


void ColdBeamSource::sampleBacaBoris(std::string namafile)
{
    std::ifstream baca(namafile.c_str());
    if (!baca.is_open())
    {
        std::cerr << "tidak dapat membuka file ->" << namafile.c_str() << std::endl;
        return;
    }
    while (true)
    {
        double z, r, vz, vr, vtheta;
        baca >> z >> r >> vz >> vr >> vtheta;
        //ini hati-hati sama bobotnya. kalau misalnya beda, bisa salah.
        spec.addParticle2({z,r,0},{vz,vr,vtheta}, spec.mpw0);
        if (baca.eof())
        {
            break;
        }
    }
}

void ColdBeamSource::sampleBacaPartialBoris(std::string namafile, double fraction)
{
    std::ifstream baca(namafile.c_str());
    if (!baca.is_open())
    {
        std::cerr << "tidak dapat membuka file ->" << namafile.c_str() << std::endl;
        return;
    }
    while (true)
    {
        double z, r, vz, vr, vtheta;
        baca >> z >> r >> vz >> vr >> vtheta;
        //ini hati-hati sama bobotnya. kalau misalnya beda, bisa salah.
        double undi = rnd();
        if (undi < fraction)
        {
            spec.addParticle2({z,r,0},{vz,vr,vtheta}, spec.mpw0); //nggak perlu digeser kedepan setengah langkah. Soalnya ini udah kegeser.
        }
        if (baca.eof())
        {
            break;
        }
    }
}

void ColdBeamSource::sampleBacaCyclotronic(std::string namafile)
{
    std::ifstream baca(namafile.c_str());
    if (!baca.is_open())
    {
        std::cerr << "tidak dapat membuka file ->" << namafile.c_str() << std::endl;
        return;
    }
    while (true)
    {
        double z, r, vz, vr, vtheta;
        baca >> z >> r >> vz >> vr >> vtheta;
        //ini hati-hati sama bobotnya. kalau misalnya beda, bisa salah.
        spec.addParticle3({z,r,0},{vz,vr,vtheta}, spec.mpw0); //nggak perlu digeser kedepan setengah langkah. Soalnya ini udah kegeser.
        if (baca.eof())
        {
            break;
        }
    }
}

void ColdBeamSource::sampleBacaPartialCyclotronic(std::string namafile, double fraction)
{
    std::ifstream baca(namafile.c_str());
    if (!baca.is_open())
    {
        std::cerr << "tidak dapat membuka file ->" << namafile.c_str() << std::endl;
        return;
    }
    while (true)
    {
        double z, r, vz, vr, vtheta;
        baca >> z >> r >> vz >> vr >> vtheta;
        //ini hati-hati sama bobotnya. kalau misalnya beda, bisa salah.
        double undi = rnd();
        if (undi < fraction)
        {
            spec.addParticle3({z,r,0},{vz,vr,vtheta}, spec.mpw0); //nggak perlu digeser kedepan setengah langkah. Soalnya ini udah kegeser.
        }
        if (baca.eof())
        {
            break;
        }
    }
}