#ifndef SOLVER_H
#define SOLVER_H

#include"World.h"
#include<assert.h>
#include<iostream>
template<int S>
struct Row
{
    //default constructor dulu, ada 2xS entri setiap Row
    Row()
    {
        for (int i=0; i<S; i++)
        {
            a[i] = 0;
            col[i] = -1; //-1, buat nandain kalo belum disertakan ke dalam matriks (mana ada indeks kolom -1)
        }
    }

    void operator= (const Row &O)
    {
        for (int i=0; i<S; i++)
        {
            a[i] = O.a[i];
            col[i] = O.col[i];
        }
    }

    double a[S];
    int col[S]; //nggak dynamically allocated ya
};

class Matrix
{
    public:
    //Constructor harus parametrik, jumlah barisnya langsung dispesifikasi. Dengan initializer list
        Matrix(int nr) : nu{nr}
        {
            rows = new Row<nvals>[nr]; //udah diinisiasi lewat initializer list di atas
            //std::cout << "rows[1].col[1] = " << rows[1].col[1] << "\n";
        }
        //Copy constructor juga
        Matrix(const Matrix &O) : Matrix(O.nu)
        {
            for (int i=0; i<nu; i++)
            {
                rows[i] = O.rows[i];
            }
        }
        //Destructornya juga
        ~Matrix()
        {
            if(rows) //rows nya bool? otomatis iterasi ya.
            {
                delete[] rows;
            }
        }

        //operator akses, definisi di Potentialsolver.cpp
        double& operator() (int r, int c);

        double summy();

        //operator perkalian matriks dengan vektor
        dvector operator* (dvector &v);

        //renisialisasi row
        void clearRow(int r) //r row yang dituju
        {
            rows[r] = Row<nvals>(); //ini default constructor ya, jadi a[i] = 0 dan col[i] = -1
        }

        //ada beberapa fungsi yang praktis
        Matrix diagSubstract(dvector &diagonal); //mengurangi bagian diagonal matriks dengan vektor
        Matrix diagInverse();  //inverse dari bagian diagonal matriks
        double multRow(int r, dvector &v); //perkalian vektor dengan salah satu baris matriks (otomatis hasilnya double)
        void printMatrix();
        const int nu;
        static constexpr int nvals = 5;
    protected:
        Row<nvals> *rows;
};
enum SolverType {GS, PCG, STANDARD};
//solvernya sementara PCG dulu, kalau ternyata ada kasus yang nggak konvergen, nanti GS ditambahkan
class Potentialsolver
{
    public:
        //langsung parametric constructor
        Potentialsolver(World &world, SolverType type, int max_it, double tol) : 
        world(world), solver_type{type}, max_solver_it(max_it), tolerance(tol), A(world.nz*world.nr)
        {
            buildMatrix();
        }
        //fungsi buildmatrix heptadiagonal (tapi ini dua dimensi ding, berarti pentadiagonal ya)
        void buildMatrix();
        
        bool solve()
        {
    	    switch(solver_type)
    	    {
    		    case GS: return solveGSLinear();
    		    case PCG: return solvePCGLinear();
                case STANDARD: return solveGSStandard();
    		    default: return false;
    	    }
        }
       
        //untuk menghitung medan listrik
        void computeEF();
        //nanti coba cek mekanisme set tegangan boundary yang lewat addBound

        //di sini kita letakkan parameter terkait dengan rangkaian eksternal
        void setCircuitParameter(double v, double r, double sigma)
        {
            sourcePotential = v;
            resistance = r;
            chargeSurfaceDensity = sigma; //ini syarat awal dari rapat muatan. Kalau nggak diketahui, nilai nol masih masuk akal [131022]
        }
        //untuk menghitung rapat muatan katoda setelah nilai tegangan katoda diketahui (untuk digunakan lagi pada iterasi berikutnya) []
        void calcSigma(); //sekalian buat ngitung arusnya juga [181022]
        double getCurrent()
        {
            return externalCurrent;
        }
        double getVoltage()
        {
            return sourcePotential;
        }
        double getResistance()
        {
            return resistance;
        }
    protected:
        World &world; //world 2 dimensi, harus hati-hati
        Matrix A; //ini bikin matriks baru
        SolverType solver_type;
        bool solveNR();
        bool linearSolvePCG(Matrix &A, dvector &x, dvector &b);
        bool linearSolveGS(Matrix &A, dvector &x, dvector &b);
        //fungsi untuk menyelesaikan PCG linear
        bool solvePCGLinear();
        bool solveGSLinear();
        bool solveGSStandard();
        //beberapa jenis node
        enum NodeType {REGULAR, NEUMANN, DIRICHLET, CHARGED, AXIS}; //tipe charged ditambahkan ketika rapat muatan permukaan pada suatu konduktor diketahui [041022]
        std::vector<NodeType> node_type; //vektor, karena dimensinya sudah dideflate
        unsigned max_solver_it;
        double tolerance;

        //di bawah ini adalah parameter sirkuit [131022]
        double sourcePotential = 0; //ini positif, nanti tegangan katodanya menyesuaikan [131022]
        double resistance = 1;
        double chargeSurfaceDensity = 0; //rapat muatan pada katoda [131022]
        double externalCurrent = 0; //ini untuk ngitung arus yang lewat di resistor [181022]

};

#endif