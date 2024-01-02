#ifndef WORLD_H
#define WORLD_H

#include "Field.h"
#include<random>
#include<vector>
#include<chrono>
#include<string>
#include<iostream>
class Species; //why do we need this?

//mula-mula definisikan class random number generator
//algoritma yang digunakan adalah Mersenne twister
//Seed dari mersenne twister ditentukan oleh random_device, yang tergantung pada kondisi komputer
class Rnd
{
public:
    //constructor, pakai default. Digunakan untuk menentukan seed awal dan batas distribusi (pakai initializer list)
    Rnd() : mt_gen{std::random_device()()}, rnd_dist{-1.0,1.0} {}
    //overload operator ()
    double operator() () 
    {
        return  rnd_dist(mt_gen);
    }
protected:
    std::mt19937 mt_gen; //jadi dari sini, jenis mt_gen itu mt19937 (bukan sekedar int atau double)
    std::uniform_real_distribution<double> rnd_dist; //ini yang nanti digunakan untuk menentukan posisi awal partikel
};

extern Rnd rnd; //memberitahu compiler bahwa objek berjenis Rnd dengan nama rnd didefinisikan di tempat lain.

namespace Const
{
    const double EPS_0 = 8.8541878e-12; //epsilon nol
    const double QE = 1.602176565e-19; //muatan elektron coulomb (tanda masih positif)
    const double AMU = 1.660538921e-27; //atomic mass unit kg
    const double ME = 9.10938215e-31; //massa elektron kg
    const double K = 1.380648e-23; //konstanta Boltzmann
    const double PI = 3.141592653; //pi
    const double EvToK = QE/K; //electron volt ke kelvin
    const double SPEEDC = 299792458; //kecepatan cahaya dalam meter per detik
}

struct Ion
{
    Ion(double3 x) : pos{x} {} //constructornya gini aja.
    double3 pos;
};

//bikin satu struct lagi kali ya, untuk elektron yang menumbuk dinding
struct Walls
{
    Walls(double3 x, double3 v) : pos{x}, vel{v} {}
    double3 pos;
    double3 vel;
};

//world 2 dimensi, sementara tanpa magnet dulu (tanpa silinder juga, kebetulan sudah silinder kan domainnya)
class World
{
    //bikin constructornya dulu (parametrik)
    public:
    World(int nz, int nr); //jangan sampai salah di implementasinya

    void setextents(const double2 x0, const double2 xm);

    double2 getX0() const
    {
        return double2(x0);
    }
    double2 getXm() const
    {
        return double2(xm);
    }
    double2 getXc() const
    {
        return double2(xc);
    }
    double2 getDh() const
    {
        return double2(dh);
    }

    const int nn[2];
    const int nz,nr; //ini nggak ditaruh di bagian protected ya, biar gampang akses
    
    //bagian perhitungan muatan pada katoda [071022] ----------------------------------------------------
    void calcTopArea(double r0, double rm)
    {
        topArea = Const::PI*(rm*rm-r0*r0);
    }
    double getTopArea()
    {
        return topArea;
    }
    void addChargeIncrement(double charge, double mpw)
    {
        cathodeCharge += charge*mpw; //ini cukup sampai sini ding []
    }
    void setCathodeCharge(double charge)//Setiap iterasi, set cathode charge = 0 [131022]
    {
        cathodeCharge = charge; 
    }
    double getChargeIncrement()
    {
        return cathodeCharge;
    }
    //----------------------------------------------------------------------------------------------------


    //ini yang menyebabkan kita include Field.h, jadi medan langsung diinisiasi di sini
    //mungkin sekalian juga kan soalnya double3, field3 juga diinisiasi di sini
    DField phi;
    DField rho;
    DField node_vol;
    DField3 ef;
	DFieldI object_id;
    
    //tambahin XtoL (dua dimensi)
    double2 XtoL(double2 x) const
    {
        double2 lc;
        lc[0] = (x[0]-x0[0])/(dh[0]);
        lc[1] = (x[1]-x0[1])/(dh[1]);
        return lc;
    }
	
    //ini untuk koefisien produksi elektron sekunder yang sederhana [131022]

	void setGamma(double gam)
    {
        secondElecConst = gam;
    }
    double getGamma()
    {
        return secondElecConst;
    }

    void addBound(int gap);//buat ngeflag seluruh dinding batas.
    void addBoundCharged(int gap);//diketahui nilai rapat muatan permukaan pada katoda [051022]

    //untuk menentukan apakah partikel ada di luar batas atau enggak
    bool zoneOneBound(double2 pos);
    bool zoneTwoBound(double2 pos);
    bool zoneThreeBound(double2 pos);
    bool zoneFourBound(double2 pos);
    bool zoneFiveBound(double2 pos);
    bool zoneSixBound(double2 pos);
    bool zoneSevenBound(double2 pos);
    bool zoneEightBound(double2 pos);

    //untuk menentukan titik tumbuk partikel dengan batas domain
    void reflectionSheathOne(double2 posawal, double2& posakhir, double& dt, double3 &vel, int& dinding);
    void reflectionSheathTwo(double2 posawal, double2& posakhir, double& dt, double3 &vel, int& dinding);
    void reflectionSheathThree(double2 posawal, double2& posakhir, double& dt, double3 &vel, int& dinding);
    void reflectionSheathFour(double2 posawal, double2& posakhir, double& dt, double3 &vel, int& dinding);

    void getObjectId();

    void setVolt(double volt);

	//ada fungsi yang digunakan untuk mengubah posisi grid dari koordinat logika menuju ke koordinat real, sebagai berikut
	double2 pos(double2 li)
	{
		double2 x;
		x[0] = x0[0] + dh[0]*li[0];
		x[1] = x0[1] + dh[1]*li[1];
		return x;
	}
	//dapat dilihat kalau fungsi di atas inputnya adalah double precision. misal diinginkan input berupa integer, maka digunakan
	
	double2 pos(int z, int r)
	{
		double2 x{(double)z, (double)r};
		return pos(x); //di sini pos merujuk ke yang argumennya double
	}		
	
	//ini buat nentukan kalo posisi yang dituliskan ada di dalam bola aatau enggak
	//bool inSphere(double3 x);
	
	//ini buat ngeflag di bidang z=0, sekaligus ngeset tegangannya (syarat batas dirichlet)
	//void addInlet(double voltage);
	
    //perhitungan rapat muatan ditaruh di sini (mengambil dari scatter di class Field)
    void computeChargeDensity(std::vector<Species> &species); //maksudnya species apa spesifik ke class species? atau ini sembarang aja?

    bool inBounds(double2 pos)
    {
        for(int i=0; i<2; i++)
        {
            if(pos[i]<=x0[i] || pos[i]>=xm[i])
            {
                return false;
            }
        }
        return true;
    }

    //pergerakan waktu ditentukan di sini
    //set ukuran langkah waktu dan lebar langkahnya
    void setTime(double dt, int num_ts)
    {
        this->dt=dt;
        this->num_ts=num_ts;
    }

    //fungsi untuk mengakses informasi langkah keberapa
    int getTs() const
    {
        return ts;
    }
    
    int getMaxiter() const
    {
    	return num_ts;
    }
    //fungsi untuk mengakes nilai waktu
    double getTime() const
    {
        return time;
    }
    //untuk mengakses lebar langkah
    double getDt() const
    {
        return dt;
    }
    //buat menghentikan iterasi
    bool isLastTimeStep() const
    {
        return ts==num_ts-1; //ini return nya true atau false kan
    }
    //majukan langkah
    bool advanceTime()
    {
        time = time+dt;
        ts = ts+1;
        return ts <=num_ts; //true atau false nya di sini ya
    }

    void setTs(int num)
    {
        ts = num;
    }

	//double getPE(); //untuk menghitung energi potensial.
	
	double getWallTime(); //untuk mengubah wall time ke satuan sekon
	
    //untuk menentukan tegangan katoda
    double calcVolt(double resistor, double teganganAwal, double charge, double mpw, int subcycling)
    {
        double V = teganganAwal - (1.5*(2*secondary_elec - elecCathodeCount) - 2*prevElec)*charge*mpw*resistor/(subcycling*dt);
        prevElec = (2*secondary_elec-elecCathodeCount);
        return V;
    }

    //bagian magnet (kalau sekedar magnet konstan mungkin tidak apa-apa ya disertakan).
    //fase magnet
    void setParam(double frek, double fase, double Mag)
    {
        this->fase = fase;
        this->frek = frek;
        this->Mag = Mag;
        this->periode = 1.0/frek;
        this->gradienwaktu = 4.0*frek;
    }
    
    double getMagStrength()
    {
    	return Mag;
    }
    
    double getMagFreq()
    {
    	return frek;
    }
    
    double getMagPhase()
    {
    	return fase;
    }

    void set_mag(double magstrength)
    {
        Mag = magstrength;
    }

    void bacamagnet(std::string nama); 

    double getR(double j)
    {
        return x0[1] + dh[1]*j;
    }
   
    //untuk menentukan sin(wt+fase)
    void sinWaktu()
    {
        sinwaktu = sin(2 * Const::PI * frek * time + fase);
    }
    
    void segitigaWaktu()
    {
    	if(std::fmod(time,periode) <= periode/2 )
    	{
    		segitigawaktu = gradienwaktu*(std::fmod(time,periode))-1;
    	}
    	else
    	{
    		segitigawaktu = -gradienwaktu*(std::fmod(time,periode) - 0.5*periode) + 1;
    	}
    }

    void setParamWall(double temp, double E_nol, double E_max_wall, double delta_max)
    {
        this->tempWall = temp;
        this->E0_wall = E_nol;
        this->Emax_wall = E_max_wall;
        this->deltaMax_wall = delta_max;
    }

    double getWallTemp()
    {
        return tempWall;
    }

    double getE0Wall()
    {
        return E0_wall;
    }

    double getEMaxWall()
    {
        return Emax_wall;
    }

    double getDeltaMax()
    {
        return deltaMax_wall;
    }

    double getTempWall()
    {
        return tempWall;
    }

    //menentukan kuat magnetnya
	
    double3 magnett(double3 pos);//uji magnet konstan tanpa interpolasi, sementara belum didefinisikan fungsinya

    //buat flattening index, kayaknya kok nggak perlu ngambil field ya.
    int index2to1(int z, int r)
    {
        return z + r*nz ;
    }
	
    //rapat gas netral
    void setParamGas(double suhu, double tekanan)
    {
        temp = suhu;
        pressure = tekanan;
        ngas = pressure/(Const::K*temp);
    }
    //untuk mengambil parameter gas
    double getTempGas()
    {
        return temp;
    }
    double getPresGas()
    {
        return pressure;
    }

    double getNgas()
    {
        return ngas;
    }

    //untuk container ion.
    std::vector<Ion> ionsHPlus;
    std::vector<Ion> ionsH2Plus;
    std::vector<Ion> ionsHMin;
    std::vector<Ion> elec_secondary;
    std::vector<Ion> elec_secondary2;
    std::vector<Walls> elec_wall;
    std::vector<Walls> elec_ionization;

    void addElecWall(double3 pos, double3 vel)
    {
        elec_wall.emplace_back(pos,vel);
        elecWallCount += 1;
    }

    void addElecCathode(double3 pos, double3 vel)
    {
        elec_wall.emplace_back(pos,vel);
        elecCathodeCount += 1;
    }

    void eraseElecWall()
    {
        elec_wall.erase(elec_wall.begin(), elec_wall.end());
    }
    
    void setElecWall(int num)
    {
        elecWallCount = num;
    }

    int getElecWallCount()
    {
        return elecWallCount;
    }

    void addIonHPlus(double3 pos)
    {
        ionsHPlus.emplace_back(pos);
        ionizationHPlus += 1;
    }

    void addIonH2Plus(double3 pos)
    {
        ionsH2Plus.emplace_back(pos);
        ionizationH2Plus += 1;
    }

    void addIonHMin(double3 pos)
    {
        ionsHMin.emplace_back(pos);
        ionizationHMin += 1;
    }

    //untuk elektron yang terbentuk karena proses ionisasi
    void addElecIon(double3 pos, double3 vel)
    {
        elec_ionization.emplace_back(pos, vel);
        elecIonization += 1;
    }
    void setElecIonization(int num)
    {
        elecIonization = num;
    }
    void eraseElecIon()
    {
        elec_ionization.erase(elec_ionization.begin(), elec_ionization.end()); //mestinya ini bisa langsung ditulis di sample ya, biar nggak kepanjangan di world.h nya [151022]
    }

    void eraseIon()
    {
        ionsHPlus.erase(ionsHPlus.begin(), ionsHPlus.end());
        ionsH2Plus.erase(ionsH2Plus.begin(), ionsH2Plus.end());
        ionsHMin.erase(ionsHMin.begin(), ionsHMin.end());
    }

    int getIonHPlus()
    {
        return ionizationHPlus;
    }

    int getIonH2Plus()
    {
        return ionizationH2Plus;
    }

    int getIonHMin()
    {
        return ionizationHMin;
    }

    void setIon(int num) //perlu untuk balik ke nol, jadi dirapel aja ketiganya
    {
        ionizationHPlus = num;
        ionizationH2Plus = num;
        ionizationHMin = num;
    }

    void addElectron(double3 pos)
    {
        elec_secondary.emplace_back(pos);
        secondary_elec += 1;
    }

    void addElectron2(double3 pos)
    {
        elec_secondary2.emplace_back(pos);
        secondary_elec +=1;
    }

    void eraseElectron()
    {
        elec_secondary.erase(elec_secondary.begin(), elec_secondary.end());
        elec_secondary2.erase(elec_secondary2.begin(), elec_secondary2.end());
    }

    int getElec()
    {
        return secondary_elec;
    }

    int getElecCathode()
    {
        return elecCathodeCount;
    }

    void setElecCathode(int num)
    {
        elecCathodeCount = num;
    }

    void setElec(int num)
    {
        secondary_elec = num;
    }
    
    double getVoltage()
    {
        return voltage;
    }

    void addElCounter()
    {
        elecCounter += 1;
    }

    void addIonCounter()
    {
        ionCounter += 1;
    }

    int getElCounter()
    {
        return elecCounter;
    }

    int getIonCounter()
    {
        return ionCounter;
    }

    int getGap()
    {
        return gapPosition;
    }

    int getEjectedElectron() //untuk melihat seberapa banyak elektron yang muncul dari proses ionisasi [210922]
    {
        return ionizationElectron;
    }
    void addCountIonizeElectron()
    {
        ionizationElectron += 1;
    }
    void resetIonizationElectron() //untuk reset counter [210922]
    {
        ionizationElectron = 0;
    }
    void deleteKanan()
    {
        jumlahKanan += 1;
    }
    void deleteAtas()
    {
        jumlahAtas += 1;
    }
    int getDeleteKanan()
    {
        return jumlahKanan;
    }
    int getDeleteAtas()
    {
        return jumlahAtas;
    }
    void addTooHigh()
    {
        tooHigh += 1;
    }
    void getTooHigh()
    {
        std::cout << "\n Jumlah partikel dengan energi di atas ambang " << tooHigh << "\n";
        tooHigh = 0;
    }

    protected:
    double x0[2];
    double xm[2];
    double xc[2];
    double dh[2];
    double hd[2];

	double2 sphere_x0 {0,0};
	double sphere_rad2 = 0;
    //lebar langkah dan ukuran langkah ditaruh di sini
    double dt = 0;
	double time = 0; //waktu awal
    int num_ts = 0; //ini maksimal iterasinya ya
	int ts =-1;

    /* tanpa silinder
    double2 base_cyl{0,0};
    double rad_cyl = 0;
    double height_cyl = 0;
    */
	std::chrono::time_point<std::chrono::high_resolution_clock> time_start; //waktu awal untuk menentukan lama simulasi
	
    void computeNodeVolumes(); //ini ditaruh di protected (jadi yang bisa akses cuma member class ini, dalam hal ini node_vol)

    //buat magnet
    double fase = 0; //fase magnet
    double frek = 0; //frekuensi magnet (sementara lewat sini dulu sekedar untuk demonstrasi)
    double periode = 0; //periode osilasi magnet, T = 1/frekuensi
    double sinwaktu = 0; //menggambarkan sin(wt+fase)
    double segitigawaktu = 0; //perubahan kuat magnet dengan sinyal segitiga.
	double gradienwaktu = 0; //mempercepat perhitungan segitigawaktu
    //dengan asumsi bahwa wilayah yang ada magnetnya berbentuk balok
    
    double Mag = 0; //amplitudo medan magnet dalam satuan tesla

    //parameter terkait gas netral
    double temp = 1;
    double pressure = 1;
    double ngas = 1;

    int tooHigh = 0;

    //untuk ionisasi, sementara pakai trik ini
    int ionizationHPlus = 0;
    int ionizationH2Plus = 0;
    int ionizationHMin = 0; 
    int elecIonization = 0;
    int secondary_elec = 0;
    int elecWallCount = 0;
    int elecCathodeCount = 0;
    double voltage;
    int ionizationElectron = 0;

    int ionCounter = 0;
    int elecCounter = 0;

    //besaran terkait dinding
    double tempWall = 0;
    double E0_wall = 0;
    double Emax_wall = 0;
    double deltaMax_wall = 0;

    int prevElec = 0;

    int gapPosition;
    int jumlahKanan = 0;
    int jumlahAtas = 0;

    double topArea = 1; //luas tutup tabung default (dalam meter) [051022]

    double cathodeCharge = 0; //awalnya mulai dari nol coulomb

    double secondElecConst = 0;
};

//kalau pakai dua macam world kok error. Mungkin kita perlu bikin class yang berbeda, tapi isinya sama persis dengan world

class Low_res
{
    public:
        Low_res(int nz, int nr) : nz{nz}, nr{nr}, nn{nz,nr}, pixel(nz,nr)  {}

        DField pixel;
        double2 XtoL(double2 x) const
        {
            double2 lc;
            lc[0] = (x[0]-x0[0])/(dh[0]);
            lc[1] = (x[1]-x0[1])/(dh[1]);
            return lc;
        }
        void setextents(const double2 x0, const double2 xm);

        double2 getX0() const
        {
            return double2(x0);
        }
        double2 getXm() const
        {
            return double2(xm);
        }
        double2 getXc() const
        {
            return double2(xc);
        }
        double2 getDh() const
        {
            return double2(dh);
        }

        const int nn[2];
        const int nz,nr; //ini nggak ditaruh di bagian protected ya, biar gampang akses

    protected:
        double x0[2];
        double xm[2];
        double xc[2];
        double dh[2];
        double hd[2];

};


#endif
