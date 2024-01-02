#include "World.h"
#include "Species.h"
#include<random>
#include<math.h>
#include "Field.h"
#include<string>
#include<fstream>
//implementasinya, dengan initializer list
using namespace std;

Rnd rnd;

//world 2 dimensi taruh di sini, nama medannya bisa sama, soalnya classnya beda
World::World(int nz, int nr) : nz{nz}, nr{nr}, nn{nz,nr}, phi(nz,nr), rho(nz,nr), node_vol(nz,nr),
	ef(nz,nr), object_id(nz,nr)
	{
		//ketika world2 diinisiasi, langsung masukkan time startnya
		time_start = chrono::high_resolution_clock::now();
		
	}

void World::setextents(double2 _x0, double2 _xm)
{
	for (int i=0; i<2; i++)
	{
		x0[i] = _x0[i];
		xm[i] = _xm[i];
		dh[i] = (xm[i]-x0[i])/(nn[i] - 1);
		xc[i] = 0.5*(xm[i] + x0[i]);
		hd[i] = 1.0/dh[i];
	}
	computeNodeVolumes();
	calcTopArea(x0[1],xm[1]);
}

//simetri Z-R, berarti perhitungan volume node perlu sedikit koreksi
//oh ya, jadi indeks i itu Z, sementara indeks j itu R
void World::computeNodeVolumes()
{
	for (int z=0; z<nz; z++)
	{
		for (int r=0; r<nr; r++)
		{
			double r2 = getR(min(r+0.5, nr-1.0)); //nggak perlu pakai if, karena nggak sekedar dikalikan setengah
			double r1 = getR(max(r-0.5, 0.0));
			double dz = dh[0];
			if (z==0 || z==nz-1)
			{
				dz *= 0.5; //untuk pojok bawah 
			}
			node_vol[z][r] = dz*Const::PI*(r2*r2-r1*r1);
		}
	}
}

void World::computeChargeDensity(std::vector<Species> &species)
{
	rho = 0;
	for (Species &sp:species)
	{
		if (sp.charge==0)
		{
			continue;
		}
		rho += sp.charge*sp.den;
	}
}

/* kayaknya nggak butuh bola, axissymmetry kan
//kalo misalnya butuh bola, bisa buat sugarcubing ini (sugarcubing 2)
void World2::addSphere(double2 x0, double radius, double phi_sphere)
{
	sphere_x0 = x0;
	sphere_rad2 = radius*radius;

	for (int i=0; i<ni; i++)
	{
		for (int j=0; j<nj; j++)
		{
			double2 x = pos(i,j);
		}
	}
}

bool World2::inCyl(double2 x)
{

}
*/

//kita bakal butuh addBound, untuk menentukan nilai tegangan di syarat batas (ada celah sedikit antara katoda dan anoda)
//i = z, j = r, ni = nz, nj = nr.
//r == 0, nggak dikasih konduktor. 
void World::addBound(int gap) //ini gap minimal 0, maksimal nz-1
{
	for (int z=0; z<nz; z++)
	{
		for (int r=0; r<nr; r++)
		{
			gapPosition = gap;
			object_id[z][r] = 0;
			if( (z<=nz-1-gapPosition) && (r==nr-1) )
			{
				object_id[z][r] = 1;
			}
			if ( z==nz-1 )
			{
				object_id[z][r] = 1;
			}
		}
	}
}

void World::addBoundCharged(int gap) //bedanya untuk katoda flagnya == 2 [051022]
{
	for (int z=0; z<nz; z++)
	{
		for (int r=0; r<nr; r++)
		{
			gapPosition = gap;
			object_id[z][r] = 0;
			if( (z<=nz-1-gapPosition) && (r==nr-1) )
			{
				object_id[z][r] = 1;
			}
			/* //ini sementara kita hapus, mau lihat efeknya
			if ( (z<=nz-1-gapPosition) && (r==nr-2) && (z>=nz-8-gapPosition) )
			{
				object_id[z][r] = 1;
			}
			*/
			if ( z==nz-1 ) //ini untuk katoda, kita flag = 2 (jadi nggak syarat batas dirichlet) [051022]
			{
				object_id[z][r] = 3;
			}
		}
	}
}

//jadi dengan bound sederhana di atas, algoritma pantulan elektronnya bisa dibuat dengan jauh lebih mudah.
bool World::zoneOneBound(double2 pos) //ini atas ya
{
	if((pos[0]>x0[0]) && (pos[0]<xm[0]) && (pos[1]>xm[1]))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool World::zoneTwoBound(double2 pos) //ini kanan
{
	if((pos[1]>x0[1]) && (pos[1]<xm[1]) && (pos[0]>xm[0]))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool World::zoneThreeBound(double2 pos) //ini bawah
{
	if((pos[0]>x0[0]) && (pos[0]<xm[0]) && (pos[1]<x0[1]))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool World::zoneFourBound(double2 pos) //kiri
{
	if((pos[1]>x0[1]) && (pos[1]<xm[1]) && (pos[0]<x0[0]))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool World::zoneFiveBound(double2 pos) //kiri atas
{
	if((pos[0]<x0[0]) && (pos[1]>xm[1]))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool World::zoneSixBound(double2 pos) //kanan atas
{
	if((pos[0]>xm[0]) && (pos[1]>xm[1]))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool World::zoneSevenBound(double2 pos) //kanan bawah
{
	if( (pos[0]>xm[0]) && (pos[1]<x0[1]))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool World::zoneEightBound(double2 pos) //kiri bawah
{
	if((pos[0]<x0[0]) && (pos[1]<x0[1]))
	{
		return true;
	}
	else
	{
		return false;
	}
}

//terkait pembalikan kecepatan, jadi perlu diingat bahwa sebenernya metode yang digunakan ini menggerakkan partikel dengan kecepatan di koordinat silinder, tetapi posisi di koordinat kartesian
//artinya kalau menumbuk dinding yang arahnya ke r, kecepatan vr (vel[1]) aja yang dibalik
//untuk syarat batas terkait simetri, aturan yang sama kayaknya nggak selalu berlaku
//selanjutnya kita cari titik potong antara posawal dengan posakhir
//nilai dinding diubah tergantung di mana partikel menumbuk dinding, 1 vertikal, 2 horizontal (aslinya 0)
void World::reflectionSheathOne(double2 posawal, double2 &posakhir, double& dt, double3 &vel, int& dinding)//dipakai untuk wilayah 5-1-6
{
	//tentukan jarak yang sudah dilintasi
	double total_segment = sqrt((posakhir[0]-posawal[0])*(posakhir[0]-posawal[0]) + (posakhir[1]-posawal[1])*(posakhir[1]-posawal[1]));
	//tentukan gradiennya
	double grad = (posakhir[1]-posawal[1])/(posakhir[0]-posawal[0]);
	//kalau grad positif, wilayah 1-6, kalau negatif wilayah 5-1
	if (grad > 0)
	{
		double za = (xm[1]-posawal[1])/grad + posawal[0];
		if (za > xm[0]) // berarti numbuk katoda (perbatasan wilayah 2, vertikal)
		{
			posakhir = {xm[0], grad*(xm[0] - posawal[0]) + posawal[1]};
			vel = {-vel[0], vel[1], vel[2]};
			dinding = 1;
		}
		else //berarti numbuk perbatasan wilayah 1. (horizontal)
		{
			posakhir = {za, xm[1]};
			vel = {vel[0], -vel[1], vel[2]};
			dinding = 2;
		}
	}
	else if (grad < 0)
	{
		double za = (xm[1] - posawal[1])/grad + posawal[0]; //dari sini langsung kelihatan kalau za lebih kecil dibanding posawal[0]
		if (za > x0[0]) //berarti numbuk selubung (horizontal)
		{
			posakhir = {za,xm[1]};
			vel = {vel[0], -vel[1], vel[2]};
			dinding = 2;
		}
		else //menumbuk permukaan 4
		{
			posakhir = {x0[0], grad*(x0[0]-posawal[0]) + posawal[1] };
			vel = {-vel[0], vel[1], vel[2]};
			dinding = 1;
		}
	}
	else //grad = 0 artinya dia mendatar, mestinya nggak numbuk.
	{
		posakhir = posakhir;
	}
	//setelah dapat nilai posakhir, tentukan dt yang tersisa.
	double first_segment = sqrt((posakhir[0]-posawal[0])*(posakhir[0]-posawal[0]) + (posakhir[1]-posawal[1])*(posakhir[1]-posawal[1]));

	dt *= 1.0-first_segment/total_segment;
}
//beuh ada empat wilayah nih, harus dikerjakan satu-satu ya.
void World::reflectionSheathTwo(double2 posawal, double2 &posakhir, double& dt, double3 &vel, int& dinding)
{
	//tentukan jarak yang sudah dilintasi
	double total_segment = sqrt((posakhir[0]-posawal[0])*(posakhir[0]-posawal[0]) + (posakhir[1]-posawal[1])*(posakhir[1]-posawal[1]));
	//tentukan gradiennya
	double grad = (posakhir[1]-posawal[1])/(posakhir[0]-posawal[0]);
	//kalau positif, berarti wilayah 6-2, kalau negatif berarti wilayah 2-7
	if (grad > 0)
	{
		double ra = grad*(xm[0] - posawal[0]) + posawal[1];
		if (ra > xm[1]) //berarti menumbuk selubung 1
		{
			posakhir = {(xm[1]-posawal[1])/grad + posawal[0] ,xm[1]};
			vel = {vel[0], -vel[1], vel[2]};
			dinding = 2;
		}
		else //berarti menumbuk katoda
		{
			posakhir = {xm[0],ra};
			vel = {-vel[0], vel[1], vel[2]};
			dinding = 1;
		}
	}
	else if(grad < 0)
	{
		double ra = grad*(xm[0] - posawal[0]) + posawal[1];
		if (ra >= x0[1]) //menumbuk katoda
		{
			posakhir = {xm[0],ra};
			vel = {-vel[0], vel[1], vel[2]}; 
			dinding = 1;
		}
		else //melewati sumbu r = 0 
		{
			//posakhir = {(x0[1]-posawal[1])/grad + posawal[0] , x0[1]};
			//vel = {vel[0], -vel[1], vel[2]};
			dinding = 2;
		}
	}
	else //grad = 0, bergerak mendatar, otomatis nabrak katoda ya
	{
		posakhir = {xm[0],posakhir[1]};
		vel = {-vel[0], vel[1], vel[2]};
		dinding = 1;
	}
	double first_segment = sqrt((posakhir[0]-posawal[0])*(posakhir[0]-posawal[0]) + (posakhir[1]-posawal[1])*(posakhir[1]-posawal[1]));

	dt *= 1.0-first_segment/total_segment;
}

void World::reflectionSheathThree(double2 posawal, double2 &posakhir, double& dt, double3 &vel, int& dinding) //sebenernya nggak sheath juga sih, ini untuk mempermudah refleksi pada sumbu simetri r = 0
{
	//tentukan jarak yang sudah dilintasi
	double total_segment = sqrt((posakhir[0]-posawal[0])*(posakhir[0]-posawal[0]) + (posakhir[1]-posawal[1])*(posakhir[1]-posawal[1]));
	//tentukan gradien terlebih dahulu
	double grad = (posakhir[1]-posawal[1])/(posakhir[0]-posawal[0]); 
	//gradien positif berarti wilayah 8-3, sementara negatif berarti wilayah 3-7
	if (grad > 0)
	{
		double ra = grad*(x0[0]-posawal[0]) + posawal[1];
		if (ra >= x0[1]) //berarti menumbuk permukaan 4.
		{
			posakhir = {x0[0], ra}; //ribet juga kalo aturannya mau konsisten sama yang sebelumnya ya
			vel = {-vel[0], vel[1], vel[2]};
			dinding = 1;
		}
		else //berarti permukaan 3. kayaknya kalo menumbuk permukaan 3 nggak perlu diapa-apakan. 
		{
			//posakhir = {(x0[1]-posawal[1])/grad + posawal[0] , x0[1]};
			//vel = {vel[0], -vel[1], vel[2]};
			dinding = 2;
		}
	}
	else if (grad < 0) //wilayah 3-7
	{
		double ra = grad*(xm[0]-posawal[0])+posawal[1];
		if (ra > x0[1]) //berarti menumbuk permukaan 2
		{
			posakhir = {xm[0], ra};
			vel = {-vel[0], vel[1], vel[2]};
			dinding = 1;
		}
		else //berarti menumbuk permukaan 3
		{
			//posakhir = {(x0[1]-posawal[1])/grad + posawal[0] , x0[1]};
			//vel = {vel[0], -vel[1], vel[2]};
			dinding = 2;
		}
	}
	else //ini gerak mendatar, mestinya nggak masuk ke dalam domain dari fungsi ini
	{
		posakhir = posakhir;
	}
	double first_segment = sqrt((posakhir[0]-posawal[0])*(posakhir[0]-posawal[0]) + (posakhir[1]-posawal[1])*(posakhir[1]-posawal[1]));

	dt *= 1.0-first_segment/total_segment;
}
void World::reflectionSheathFour(double2 posawal, double2 &posakhir, double &dt, double3 &vel, int& dinding)
{
	//tentukan jarak yang sudah dilintasi
	double total_segment = sqrt((posakhir[0]-posawal[0])*(posakhir[0]-posawal[0]) + (posakhir[1]-posawal[1])*(posakhir[1]-posawal[1]));
	//tentukan gradien
	double grad = (posakhir[1]-posawal[1])/(posakhir[0]-posawal[0]);
	if (grad > 0) //wilayah 8-4
	{
		double ra = grad*(x0[0]-posawal[0]) + posawal[1];
		if(ra >= x0[1]) //menumbuk permukaan 4
		{
			posakhir = {x0[0],ra};
			vel = {-vel[0], vel[1], vel[2]};
			dinding = 1;
		}
		else //menumbuk permukaan 3
		{
			//posakhir = {(x0[1]-posawal[1])/grad + posawal[0] ,x0[1]};
			//vel = {vel[0], -vel[1], vel[2]};
			dinding = 2;
		}
	}
	else if (grad < 0) //wilayah 4-5
	{
		double ra = grad*(x0[0]-posawal[0]) + posawal[1];
		if (ra > xm[1]) //menumbuk permukaan 1
		{
			posakhir = {(xm[1]-posawal[1])/grad + posawal[0] , xm[1]};
			vel = {vel[0], -vel[1], vel[2]};
			dinding = 2;
		}
		else //menumbuk permukaan 4
		{
			posakhir = {x0[0], ra};
			vel = {-vel[0], vel[1], vel[2]};
			dinding = 1;
		}
	}
	else //bergerak mendatar
	{
		posakhir = {x0[0], posakhir[1]};
		vel = {-vel[0], vel[1], vel[2]};
		dinding = 1;
	}
	double first_segment = sqrt((posakhir[0]-posawal[0])*(posakhir[0]-posawal[0]) + (posakhir[1]-posawal[1])*(posakhir[1]-posawal[1]));

	dt *= 1.0-first_segment/total_segment;
}
void World::getObjectId()
{
	ofstream fout("tes_object_id.dat");
	for (int z=0; z<nz; z++)
	{
		for (int r=0; r<nr; r++)
		{
			fout << z << "\t" << r << "\t" << object_id[z][r] << "\n";
		}
	}
	fout.close();
}

//berhubung belum pakai simetri cermin di bidang z=0, maka diset tegangan di tutup atas sama tutup bawah
void World::setVolt(double tegangan)
{
	for(int r=0; r<nr; r++)
	{
		phi[nz-1][r] = tegangan;
		//phi[nz-2][r] = tegangan; //lebih tebal kali ya
		voltage = tegangan;
	}
}

double World::getWallTime()
{
	auto time_now = chrono::high_resolution_clock::now();
	chrono::duration<double> time_delta = time_now - time_start;
	return time_delta.count(); //perlu dipelajari lagi library chrono ini																			
}

//sementara, asumsi magnet homogen, merata di seluruh domain. Medan tetep vektor 3 komponen.
//ini nanti cross productnya harus hati-hati, urutannya (Z,R,phi)

void Low_res::setextents(double2 _x0, double2 _xm)
{
	for (int i=0; i<2; i++)
	{
		x0[i] = _x0[i];
		xm[i] = _xm[i];
		dh[i] = (xm[i]-x0[i])/(nn[i] - 1);
		xc[i] = 0.5*(xm[i] + x0[i]);
		hd[i] = 1.0/dh[i];
	}
}
