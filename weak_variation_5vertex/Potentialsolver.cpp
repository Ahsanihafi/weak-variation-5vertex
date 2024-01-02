#include<math.h>
#include<iostream>
#include"World.h"
#include"Potentialsolver.h"
#include"Field.h"
#include<stdlib.h>
#include<string.h>
#include <fstream>
#include<iomanip>
using namespace std;
using namespace Const;

using dvector = vector<double>;
/*
//pertama, overload operator perkalian matriks dengan vektor
dvector Matrix::operator* (dvector &v)
{
    dvector r(nu);
    for (int u=0; u<nu; u++)
    {
        auto &row = rows[u]; //sekedar mempermudah
        r[u] = 0;
        for (int i=0; i<nvals; i++)
        {
            if (row.col[i] >= 0)
            {
                r[u] += row.a[i]*v[row.col[i]];
            }
            else
            {
                break;
            }
        }
    }
    return r;
}
*/
dvector Matrix::operator*(dvector &v) {
	dvector r(nu);
	for (int u=0;u<nu;u++) {
		auto &row = rows[u];
		r[u] = 0;
		for (int i=0;i<nvals;i++){
			if (row.col[i]>=0) r[u]+=row.a[i]*v[row.col[i]];
			else break;	//end at the first -1
		}
	}
	return r;
}

double& Matrix::operator() (int r, int c)
{
    auto &row = rows[r];
    int i; //i nya dideclare di sini, nggak di dalam 
    for (i=0; i<nvals; i++)
    {
        //cout << "row.col[" << i << "]= \t" << row.col[i] << "\t row.a[" << i << "]= \t" << row.a[i] <<"\n" ;
        if (row.col[i] == c)
        {
            break; //mencari i berapa supaya row.col[i] == c
        }
        //ini if kedua ya, bukan else if.
        if (row.col[i]<0)
        {
            row.col[i] = c; //ini kalau ternyata belum diassign, sekalian assign. Jadi operator akses bisa bisa buat akses
            break;
        }
    }
    //cout << i << "\n";
    assert(i!=nvals);
    return row.a[i];
}

double Matrix::summy()
{
    double summ = 0;
    for (int i=0; i<nu; i++)
    {
        for (int j = 0; j<nvals; j++)
        {
            if (rows[i].col[j] >= 0)
            {
                summ += rows[i].a[j];
            }
            else
            {
                continue;
            }
        }
    }
    return summ;
}

void Matrix::printMatrix()
{
    ofstream fout("matrix.csv");
    assert(fout.is_open());
    for (int i=0; i<nu; i++)
    {
        fout << i ;
        for (int j=0; j<nvals; j++)
        {
            if (rows[i].col[j] >= 0)
            {
                fout << ","  << rows[i].a[j] << "col" << rows[i].col[j] ;
            }
            else
            {
                continue;
            }
        }
        fout << "\n";
    }
    fout.close();
}

//diagInverse untuk preconditioning (outputnya matriks, bukan vektor. Kalau mau vektor bisa juga sih)
Matrix Matrix::diagInverse()
{
    Matrix M(nu); //buat matriks baru kosongan
    for (int i=0; i<nu; i++)
    {
        M(i,i) = 1.0/(*this)(i,i); //this itu terkait dengan class yang sudah dipanggil
    }
    return M;
}

//mengurangkan komponen diagonal matriks dengan suatu vektor
Matrix Matrix::diagSubstract(dvector &diagonal)
{
    Matrix M(nu);
    for (int i=0; i<nu; i++)
    {
        M(i,i) = 1.0/(*this)(i,i) - diagonal[i];
    }
    return M;
}

//perkalian satu baris dari matriks dengan vektor
double Matrix::multRow(int r, dvector &v)
{
    double summ = 0;
    auto &row = rows[r];
    //cout << "indeks ke = " << r << "\t" ;
    for (int i=0; i<nvals; i++)
    {
        //cout << "i = " << i << "\t" << "kolom =" << row.col[i] << "\t";
        if (row.col[i]>=0)
        {
            summ += row.a[i]*v[row.col[i]];
        }
        else
        {
            break;
        }
    }
    //cout << "\n";
    return summ;
}

//overload beberapa operator (ini di luar matrix ya?)
//overload operator di luar class itu gimana ya (mungkin general ya. yang jelas tergantung sama argumennya apa)
dvector operator+ (const dvector&a , const dvector &b)
{
    size_t n = a.size();
    dvector v(n);
    for (size_t i=0; i<n; i++)
    {
        v[i] = a[i] + b[i];
    }
    return v;
}
dvector operator- (const dvector&a , const dvector &b)
{
    size_t n = a.size();
    dvector v(n);
    for (size_t i=0; i<n; i++)
    {
        v[i] = a[i] - b[i];
    }
    return v;
}
/*
dvector operator* (const dvector &b, const dvector &a) //lumrahnya skalar dikalikan vektor ya, bukan sebaliknya
{
    size_t n = a.size();
    dvector v(n);
    for (int i=0; i<n; i++)
    {
        v[i] = a[i] + b[i];
    }
    return v;
}
*/
dvector operator* (const dvector &a, const double s)
{
	size_t n=a.size();
	dvector v(n);
	for (size_t i=0; i<n; i++)
	{
		v[i] = a[i]*s;
	}
	return v;
}

dvector operator* (const double s, const dvector &a)
{
	size_t n=a.size();
	dvector v(n);
	for (size_t i=0; i<n; i++)
	{
		v[i] = a[i]*s;
	}
	return v;
}

namespace vec
{
    //kumpulan fungsi, manggilnya pakai vec::fungsi
    double dot (dvector a, dvector b)
    {
        size_t n = a.size();
        double summ = 0;
        for (size_t i=0; i<n; i++)
        {
            summ += a[i]*b[i];
        }
        return summ;
    }

    double norm(dvector a)
    {
        size_t n = a.size();
        double summ = 0;
        for (int i=0; i<n; i++)
        {
            summ += a[i]*a[i];
        }
        return sqrt(summ);
    }

    //mengubah matriks rank 2 ke rank 1
    dvector deflate(DField &f)
    {
        dvector v(f.nz*f.nr);
        for(int z=0; z<f.nz; z++)
        {
            for (int r=0; r<f.nr; r++)                
            {
                v[f.index2to1(z,r)] = f[z][r]; //ini ngisinya nggak berurutan, tapi hasil akhirnya urut
            }
        }
        return v;
    }
    
    void inflate(DField &f, dvector &v) //ubah field yang udah ada, nggak perlu bikin yang baru
    {
        for (int z = 0; z<f.nz; z++)
        {
            for (int r = 0; r<f.nr; r++)
            {
                f[z][r] = v[f.index2to1(z,r)];
            }
        }
    }
    void tampilkan(dvector &v)
	{
		size_t nmax = v.size();
		for (int n=0; n<nmax; n++)
		{
			cout << v[n] << "\t";
		}
		cout << "\n";
	}
    void printkan(dvector &v)
    {
        size_t nmax = v.size();
        ofstream fout("printkanx.csv");
        assert(fout.is_open());
        for (int i=0; i<nmax; i++)
        {
            fout << i << "," << v[i] << "\n";
        }
        fout.close();
    }
};

//selanjutnya masuk ke bagian class potential solver. 
//Pertama bikin matriks pentadiagonal. 
//Untuk simetri Z-R, ketika r=0, harus diterapkan syarat batas neumann, kalau enggak bisa blow up.

void Potentialsolver::buildMatrix()
{
    double2 dh = world.getDh();
    double idz = 1.0/dh[0];
    double idr = 1.0/dh[1];
    double idz2 = idz*idz;
    double idr2 = idr*idr;
    int nz = world.nz;
    int nr = world.nr;
    int nf = nz*nr; //flattened index
    node_type.reserve(nf);
    for (int r = 0; r < nr; r++)
    {
        for (int z = 0; z < nz; z++)
        {
            int id1 = world.index2to1(z,r); //z dulu lho, pastikan yang di atas ikut juga.
            A.clearRow(id1); 
            if (world.object_id[z][r] == 1) //artinya logam
            {
                //cout << id1 << "Dirichlet \n";
                A(id1,id1) = 1;
                node_type[id1] = DIRICHLET;
                continue; //lanjut iterasi berikutnya, yang di bawah diskip semua
            }
            
            if (world.object_id[z][r] == 3) //ini ketika ada rapat muatan permukaan pada logam, sementara cuma di pojokan aja [041022]
            {
                if (z==0)
                {
                    //selain untuk z==nz-1 didefinisikan belakangan [051022]
                }
                else if (z==nz-1) //yang penting di sini, lokasi katoda [051022], kita ngikuti skema verboncoeur [121022]
                {
                    A(id1,id1) = idz + (world.getDt()/world.getTopArea()*resistance*EPS_0);
                    A(id1,id1-1) = -idz;
                }
                else if(r==0)
                {
                    //selain untuk z==nz-1 didefinisikan belakangan [051022]
                }
                else if (r==nr-1)
                {
                    //selain untuk z==nz-1 didefinisikan belakangan [051022]
                }
                node_type[id1] = CHARGED;
                continue; //lanjut iterasi berikutnya, yang di bawah diskip semua
            }
            
            node_type[id1] = NEUMANN; //default kalau nggak dikasih addBound (mungkin di ujung z dibuka aja)

            if(z == 0)
            {
                //cout << id1 << "Neumann1 \n";
                A(id1,id1) = idz;
                A(id1,id1+1) = -idz;
            }
            else if(z == nz-1)
            {
                //cout << id1 << "Neumann2 \n";
                A(id1,id1) = idz;
                A(id1,id1-1) = -idz;
            }
            else if(r == 0) //bagian ini masih agak meragukan, tapi coba dulu, dibenchmark.
            {
                /*
                //cout << id1 << "Neumann3a \n";
                A(id1,id1) = idr;
                //cout << id1 << "Neumann3b \n";
                A(id1,id1+nz) = -idr;
                */
                A(id1,id1+nz) = 4*idr2;
                A(id1,id1) = -2.0*(2*idr2+idz2);
                A(id1,id1+1) = idz2;
                A(id1,id1-1) = idz2;

                node_type[id1] = AXIS;
            }
            else if(r==nr-1)
            {
               // cout << id1 << "Neumann4 \n";
                A(id1,id1) = idr;
                A(id1,id1-nz) = -idr;
            }

            else //ini untuk node yang nggak di ujung
            {
                //cout << id1 << "Regular \n";
                A(id1,id1-nz) = idr2*(1.0-0.5/r); //ini r mulai dari 1. r=0 dianggap syarat batas neumann. 
                A(id1,id1-1) = idz2;
                A(id1,id1) = -2.0*(idr2+idz2);
                A(id1,id1+1) = idz2;
                A(id1,id1+nz) = idr2*(1.0+0.5/r);
                node_type[id1] = REGULAR; //jangan lupa flag jenis node nya.
            }
        }
    }
    //untuk memastikan beneran kebuild apa enggak, 
    cout << "Matrix is built successfully " << endl;
    //cout << "summy =" << A.summy() << endl;
    //A.printMatrix();
    dvector ax = vec::deflate(world.phi);
    dvector ab = vec::deflate(world.rho);
    linearSolveGS(A,ax,ab);
    vec::inflate(world.phi,ax);
}

//solver PCG linear, pengalaman kemarin, PCG linear nggak bisa dipakai kalau di domain nggak ada muatan sama sekali
bool Potentialsolver::solvePCGLinear()
{
    bool converged = false; //by default
    int nu = A.nu; //ambil ukuran matriksnya

    //matriks kolom
    dvector x = vec::deflate(world.phi);
    dvector b = vec::deflate(world.rho); 
    //misal ada konduktor dalam domain, set matriks kolom ruas kanan (rho) sama dengan matriks ruas kiri (tegangan)
    for (int u = 0; u<nu; u++)
    {
        if (node_type[u] == DIRICHLET)
        {
            b[u] = x[u];
        }
        else if(node_type[u] == NEUMANN)
        {
            b[u] = 0; //muatan terlokalisir, muatan di pusat dianggap nol juga waktu perhitungan nih (semoga gak ngaruh ya?)
        }
        else if(node_type[u] == CHARGED)
        {
            //kita pakai skema verbonce, modifikasi syarat batas. Tapi menggunakan parameter arus, bukan muatan [131022]
            b[u] = (chargeSurfaceDensity + world.getChargeIncrement()/world.getTopArea() + sourcePotential*world.getDt()/(world.getTopArea()*resistance) + 0.5*b[u]*world.getDh()[0])/EPS_0; //ini tandanya positif, bagian matriksnya sudah disesuaikan [071022]
            //ini berarti rho nya bisa dianggap berubah gitu ya. Ini skemanya dikopi buat solver yang lain juga.
        }
        else
        {
            b[u] = -b[u]/EPS_0;
        }
    }

    //untuk menentukan konvergensi dari PCG 
    double l2 = 0;
    Matrix M = A.diagInverse(); //untuk preconditioning
    //inisiasi iterasi
    dvector g = A*x - b; //residu awal
    //cout << g.size() << "\t";
    dvector s = M*g; //preconditioning
    dvector d = -1.0*s; //arah gradien (negatif)
    for (unsigned it = 0; it<max_solver_it; it++)
    {
        //vec::tampilkan(g);
        dvector z = A*d; //belum bisa di luar kepala nih, nanti dipelajari lagi
        //di bawah ini mekanisme untuk mencari konjugat (kalau nggak salah gram schimdt orthogonalization)
        double alpha = vec::dot(g,s); 
        double beta = vec::dot(d,z);
        x = x + (alpha/beta)*d;
        g = g + (alpha/beta)*z; //ini harus dipasikan kalau beta nggak nol

        s = M*g;
		
		//cout << "pos 1" << endl;
        
		beta = alpha;
        alpha = vec::dot(g,s);

        d = (alpha/beta)*d-s; //ini arah baru ya
		
		//cout << "pos 2" << endl;
		
        l2 = vec::norm(g); //norma residu yang baru
        if (l2<tolerance)
        {
            converged = true;
            break;
        }
    }
    if (!converged) //ini artinya converged == false?
    {
        cerr << "PCG gagal konvergen dengan norma(g)= " << l2 << ", kerjakan ulang dengan gauss seidel \n";
        //perlu definisikan ulang kali ya
        
        dvector xx = vec::deflate(world.phi);
        dvector bb = vec::deflate(world.rho);
        linearSolveGS(A,xx,bb);
        dvector gg = A*xx-bb;
        l2 = vec::norm(gg);
        if (l2 <tolerance)
        {
            converged = true;
        }
        vec::inflate(world.phi,xx);
        
    }
    
    vec::inflate(world.phi,x);
    return converged;
}

bool Potentialsolver::linearSolveGS(Matrix &A, dvector &x, dvector &b)
{
    double l2 = 0; //nilai awal untuk menentukan norma
    bool converged = false;

    int nu = A.nu;
    for (int u = 0; u<nu; u++)
    {
        if (node_type[u] == DIRICHLET)
        {
            b[u] = x[u];
        }
        else if(node_type[u] == NEUMANN)
        {
            b[u] = 0; //muatan terlokalisir, muatan di pusat dianggap nol juga waktu perhitungan nih (semoga gak ngaruh ya?)
        }
        else if(node_type[u] == CHARGED)
        {
            b[u] = (chargeSurfaceDensity + world.getChargeIncrement()/world.getTopArea() + sourcePotential*world.getDt()/(world.getTopArea()*resistance) + 0.5*b[u]*world.getDh()[0])/EPS_0; //ini tandanya positif, bagian matriksnya sudah disesuaikan [071022]
        }
        else
        {
            b[u] = -b[u]/EPS_0;
        }
    }
    for (unsigned it = 0; it<max_solver_it; it++)
    {
        for (int u=0; u<nu; u++)
        {
            //cout << "indeks ke = " << u << "\n";
            double S = A.multRow(u,x) - A(u,u)*x[u]; //jadi dikurangi bagian diagonal
            double phibaru = (b[u] - S)/A(u,u); //bagian gauss seidel

            x[u] = x[u] +1.4*(phibaru - x[u]);
        }
        //cek konvergensi
        if(it%25==0)
        {
            dvector r = A*x - b; //kalo konvergen, ini nilainya mendekati nol (atau tergantung kita mau toleransi berapa)
            l2 = vec::norm(r);
            if( l2<tolerance )
            {
                converged = true;
                break;
            }
        }
    }
    if(!converged)
    {
        cerr << "GS untuk koreksi PCG juga gagal konvergen, l2 = " << l2 << "\n"; 
    }
    return converged;
}

bool Potentialsolver::solveGSLinear()
{
    double l2 = 0; //nilai awal untuk menentukan norma
    bool converged = false;

    int nu = A.nu;
    dvector x = vec::deflate(world.phi);
    dvector b = vec::deflate(world.rho);

    for (int u = 0; u<nu; u++)
    {
        if (node_type[u] == DIRICHLET)
        {
            b[u] = x[u];
        }
        else if(node_type[u] == NEUMANN)
        {
            b[u] = 0; //muatan terlokalisir, muatan di pusat dianggap nol juga waktu perhitungan nih (semoga gak ngaruh ya?)
        }
        else if(node_type[u] == CHARGED)
        {
            b[u] = (chargeSurfaceDensity + world.getChargeIncrement()/world.getTopArea() + sourcePotential*world.getDt()/(world.getTopArea()*resistance) + 0.5*b[u]*world.getDh()[0])/EPS_0; //ini tandanya positif, bagian matriksnya sudah disesuaikan [071022]
            //b[u] = (chargeSurfaceDensity + sourcePotential*world.getDt()/(world.getTopArea()*resistance))/EPS_0;
        }
        else if(node_type[u] == AXIS)
        {
            b[u] = -b[u]/EPS_0;
        }
        else
        {
            b[u] = -b[u]/EPS_0;
        }
    }
    for (unsigned it = 0; it<max_solver_it; it++)
    {
        for (int u=0; u<nu; u++)
        {
            //cout << "indeks ke = " << u << "\n";
            double S = A.multRow(u,x) - A(u,u)*x[u]; //jadi dikurangi bagian diagonal
            //double phibaru = (b[u] - S)/A(u,u); //bagian gauss seidel
            x[u] = (b[u] - S)/A(u,u);
            //x[u] = x[u] +1.4*(phibaru - x[u]);
        }
        //cek konvergensi
        if(it%25==0)
        {
            dvector r = A*x - b; //kalo konvergen, ini nilainya mendekati nol (atau tergantung kita mau toleransi berapa)
            l2 = vec::norm(r);
            if( l2<tolerance )
            {
                converged = true;
                break;
            }
        }
    }
    if(!converged)
    {
        cerr << "GS untuk koreksi PCG juga gagal konvergen, l2 = " << l2 << "\n"; 
    }
    vec::inflate(world.phi,x);
    return converged;
}

bool Potentialsolver::solveNR()
{
	/*main NR iteration loop*/
	const int NR_MAX_IT=20;		/*maximum number of NR iterations*/
	const double NR_TOL = 1e-3;
	int nu = A.nu;

	Matrix J(nu);
	dvector P(nu);
	dvector y(nu);
	dvector x = vec::deflate(world.phi);
	dvector b = vec::deflate(world.rho);

	/*set RHS to zero on boundary nodes (zero electric field)
      and to existing potential on fixed nodes */
    for (int u=0;u<nu;u++)
    {
		if (node_type[u]==NEUMANN) b[u] = 0;			/*neumann boundary*/
        else if (node_type[u]==DIRICHLET) b[u] = x[u];	/*dirichlet boundary*/
        else b[u] = -b[u]/EPS_0;            /*regular node*/
    }

	double norm;
	bool converged=false;
	for(int it=0;it<NR_MAX_IT;it++)
	{
		/*compute F by first subtracting the linear term */
		dvector F = A*x-b;

		/*Compute P, diagonal of d(bx)/dphi*/
		for (int n=0;n<nu;n++)
		{
			if (node_type[n]==REGULAR)
				P[n] = 0;
		}

		/*Compute J = A-diag(P)*/
		Matrix J = A.diagSubstract(P);

		/*solve Jy=F*/
        bool lin_converged = false;
        lin_converged = linearSolvePCG(J,y,F);
	
		/*clear any numerical noise on Dirichlet nodes*/
		for (int u=0;u<nu;u++)
			if (node_type[u]==DIRICHLET) y[u]=0;

		/*x=x-y*/
		x = x-y;

		norm=vec::norm(y);
		//cout<<"NR norm: "<<norm<<endl;

		if (norm<NR_TOL)
		{
			converged=true;
			break;
		}
	}

	if (!converged)
		cout<<"NR+PCG failed to converge, norm = "<<norm<<endl;

	/*convert to 3d data*/
	vec::inflate(world.phi,x);
	return converged;
}

/*PCG solver for a linear system Ax=b*/
bool Potentialsolver::linearSolvePCG(Matrix &A, dvector &x, dvector &b)
{
	bool converged= false;

	double l2 = 0;
	Matrix M = A.diagInverse(); //inverse of Jacobi preconditioner

	/*initialization*/
	dvector g = A*x-b;
	dvector s = M*g;
	dvector d = -1*s;

	for (unsigned it=0;it<max_solver_it;it++)
	{
		dvector z = A*d;
		double alpha = vec::dot(g,s);
		double beta = vec::dot(d,z);

		x = x+(alpha/beta)*d;
		g = g+(alpha/beta)*z;
		s = M*g;

		beta = alpha;
		alpha = vec::dot(g,s);

		d = (alpha/beta)*d-s;
		l2 = vec::norm(g);
		if (l2<tolerance) {converged=true;break;}
	}

	if (!converged)	cerr<<"PCG failed to converge, norm(g) = "<<l2<<endl;
    return converged;
}

//ini korrdinat silinder kan, harus cek dulu kalau ada bedanya sama koordinat kartesian
void Potentialsolver::computeEF()
{
    DField &phi = world.phi;
    DField3 &ef = world.ef; //komponen theta selalu nol di mana-mana ya

    int nz = world.nz;
    int nr = world.nr;

    double2 dh = world.getDh();
    double dz = dh[0];
    double dr = dh[1];
	
	//int counter = 0;
    for (int z=0; z<nz; z++)
    {
        for (int r=0; r<nr; r++)
        {
            //jadi untuk arah medan yang berbeda, posisi finite differencenya juga beda
            if (z==0)
            {
				//cout << "pos a" << endl;
                ef[z][r][0] = -(-3*phi[z][r] + 4*phi[z+1][r] - phi[z+2][r])/(2*dz);  	
            }
            else if (z==nz-1)
            {
				//cout << "pos b" << endl;
                ef[z][r][0] = (-3*phi[z][r] + 4*phi[z-1][r] - phi[z-2][r])/(2*dz);
			}
			else 
			{
				ef[z][r][0] = -(phi[z+1][r] - phi[z-1][r])/(2*dz);
			}
			//ini untuk r, jadi kelihatannya if baru itu beda dengan yang lama (kena override?)
            if (r==0)
            {
				//cout << "pos c" << endl;
                //ef[z][r][1] = -(-3*phi[z][r] + 4*phi[z][r+1] - phi[z][r+2])/(2*dr);
                ef[z][r][1] = 0;
			}
            else if (r==nr-1)
            {
				//cout << "pos d" << endl;
                ef[z][r][1] = (-3*phi[z][r] + 4*phi[z][r-1] - phi[z][r-2])/(2*dr);
			}
            else
            {
				//cout << "pos f" << endl;
				//salahnya di sini
				ef[z][r][1] = -(phi[z][r+1] - phi[z][r-1])/(2*dr);
            }
			//cout << r << "\t" << z << endl;
			//cout << counter << endl;
			//counter += 1;
            ef[z][r][2] = 0; //komponen phi
        }
    }
}

//perhitungan rapat muatan pada katoda ini dilakukan setelah potensial (dan medan listrik sekalian ding) selesai disolve [131022]
//ini sementara untuk skemanya verbonce kan, skema asalnya masih satu dimensi. Kita perlu pelajari skema dua dimensi oleh Vahedi juga [131022]
void Potentialsolver::calcSigma()
{
    //ini kayaknya nilai tegangan di katoda untuk r beda nilainya beda, sementara mungkin kita rerata [131022]
    double phiCathode = 0;
    int nz = world.nz;
    int nr = world.nr;
    DField &phi = world.phi;
    for (int r=0; r<nr; r++)
    {
        phiCathode += phi[nz-1][r];
    }
    phiCathode = phiCathode/nr; //skema reratanya kayak gini kali [131022]
    chargeSurfaceDensity += world.getChargeIncrement()/world.getTopArea() + ((sourcePotential-phiCathode)*world.getDt())/(world.getTopArea()*resistance);
    //kita sekalian itung arusnya ya
    externalCurrent = (sourcePotential-phiCathode)/resistance;
    //cout << setprecision(12);
    //cout << "\t sigma=" << chargeSurfaceDensity << "\n" ;
}

//GSStandard itu maksudnya solver gauss seidel tanpa menggunakan matriks tridiagonal. Lebih mudah buat dicek. [03012023]
bool Potentialsolver::solveGSStandard()
{
    DField &phi = world.phi; //nggak perlu bikin kopi field baru, cukup cek alamatnya aja.
    DField &rho = world.rho;

    double2 dh = world.getDh();
    double idz = 1.0/dh[0];
    double idr = 1.0/dh[1];
    double idz2 = idz*idz;
    double idr2 = idr*idr;

    double L2 = 0;
    bool converged = false; 

    double crz = 0.5/(idz2+idr2); //untuk selain bagian syarat batas

    //iterasi di sini
    for (unsigned it = 0; it<max_solver_it; it++)
    {
        for (int i=0; i<world.nz; i++)
        {
            for (int j=0; j<world.nr; j++)
            {
                //by default, syarat batas di boundary pakai syarat batas Neumann. Kecuali kalau ada logam atau di r=0
                if (world.object_id[i][j] == 1) //logam anoda
                {
                    continue; //nggak diapa-apain
                }
                else if (world.object_id[i][j] == 3) //logam katoda
                {
                    //kita pakai algoritma Vahedi
                    double denom = world.nr*EPS_0*idz + (world.getDt()/(world.getTopArea()*resistance));
                    double sumasi1 = 0; double sumasi2 = 0;
                    for(int k=0; k<world.nr; k++)
                    {
                        sumasi1 += phi[world.nz-2][k];
                        sumasi2 += rho[world.nz-1][k];
                    }
                    double sumasi = sumasi1*EPS_0*idz + 0.5*sumasi2*dh[0];
                    phi[i][j] = (chargeSurfaceDensity + sumasi + sourcePotential*world.getDt()/(resistance*world.getTopArea()) - world.getChargeIncrement()/world.getTopArea() )/denom;
                    continue; //nggak perlu pertimbangin yang di bawah, langsung lanjut iterasi
                }
                
                if (i==0) //titik simetri z=0;
                {
                    phi[i][j] = phi[i+1][j];
                }
                else if (i==world.nz-1) //ini nggak kepakai sih, posisi katoda
                {
                    //phi[i][j] = phi[i-1][j];
                    continue;
                }
                else if (j==0) //kasus r=0, pakai teorema gauss, bisa cek di Birdsall
                {
                    phi[i][j] = ((rho[i][j]/EPS_0) + 4*phi[i][j+1]*idr2 + (phi[i+1][j]+phi[i-1][j])*idz2)/(4*idr2+2*idz2);
                }
                else if (j==world.nr-1) //ini sebagian kepakai, sebagian kena anoda
                {
                    //phi[i][j] = phi[i][j-1]; //ini salah nih, mungkin kita perlu ngikuti skema potential drop nya vahedi
                    //coba kita pakai linear voltage drop
                    phi[i][j] = phi[world.nz-1][j]*(i-(world.nz-1-world.getGap()))/world.getGap();
                }
                else //bagian dalam
                {
                    //ini kebetulan r0=0, jadi kita pakai rj = j*dr
                    double phi_baru = crz*( (rho[i][j]/EPS_0) + idz2*(phi[i+1][j] + phi[i-1][j]) + phi[i][j+1]*(idr2 + 0.5*idr2/j) + phi[i][j-1]*(idr2 - 0.5*idr2/j) ); 
                    //coba pakai SOR di koordinat silinder, mana tau bisa
                    phi[i][j] = phi[i][j] + 1.4*(phi_baru-phi[i][j]);
                    //phi[i][j] = phi_baru;
                }
            }
        }
        //periksa konvergensi solver di sini. mungkin setiap 50 iterasi aja ya. biar nggak bingung sementara tanpa syarat batas
        if (it%50 == 0)
        {
            double sum = 0;
            for (int i=1; i<world.nz-1; i++)
            {
                for (int j=1; j<world.nr-2; j++)
                {
                    double R = -2*phi[i][j]*(idr2+idz2)
                    +((rho[i][j]/EPS_0)
                    +idz2*(phi[i+1][j]+phi[i-1][j])
                    +phi[i][j+1]*(idr2+0.5*idr2/j)
                    +phi[i][j-1]*(idr2-0.5*idr2/j));
                    sum += R*R;
                }
            }
            L2 = sqrt(sum/(world.nz*world.nr)); 
            if (L2 < tolerance)
            {
                converged=true;
                break;
            }
        }
    }
    if(!converged)
    {
        cerr <<"Gauss seidel standar gagal konvergen, L2="<<L2<<endl;
    }
    //cout << setprecision(12);
    //cout << "\t sigma=" << chargeSurfaceDensity << "\t" ;
    //mungkin perhitungan rapat muatan permukaan yang baru bisa ditaruh di sini ya.
    double sumasi0 = 0; double sumasi1 = 0; double sumasi2 = 0;
    for(int k=0; k<world.nr; k++)
    {
        sumasi0 += phi[world.nz-1][k];
        sumasi1 += phi[world.nz-2][k];
        sumasi2 += rho[world.nz-1][k];
    }
    chargeSurfaceDensity = EPS_0*idz*(sumasi0-sumasi1) - 0.5*dh[0]*sumasi2;
    externalCurrent = (sourcePotential-phi[world.nz-1][3])/resistance;
    return converged;
}