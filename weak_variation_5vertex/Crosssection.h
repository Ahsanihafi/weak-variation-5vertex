#ifndef CROSS_H
#define CROSS_H

#include<fstream>
#include<vector>
#include<string>
#include"Field.h"
#include"World.h"
#include<algorithm>
//#include"Species.h"
//jadi kita tinggal baca data yang sudah diolah sebelumnya. kita perlu interpolasi, jadi kayaknya nggak bisa pakai std::vector ya, kita pakai energy.
//ukuran energy dan maxval kita tentukan menggunakan file output pada program pengolah data. 
//CrossSection di atas hanya untuk satu spesies partikel dengan gas netral. Kalau ada beberapa spesies, maka memerlukan beberapa class.
//CrossSection ada di atas species, jadi seharusnya di dalam argumen CrossSection nggak ada spesies. Nanti di main harus dideclare duluan sebelum Species.
//kita enumerasikan tipe tumbukan
enum CollisionType {Elastic, Exc_1, Exc_2, Exc_3, Null_Collision};

class CrossSection
{
    public:
        CrossSection(const int ni, World &world, std::string nama, double mass, double stepEnergi, int jumlah_proses) : ne{ni}, world{world}, nama{nama}, mass{mass}, stepEnergi{stepEnergi}, jumlahProses{jumlah_proses}
        { 
            addProcess(jumlahProses, ne);
            // misal kita inisialisasi di main aja gimana? di sini bermasalah mulu
            readTable(nama); //langsung baca file
            frequency();
            set_max_val();
            normalize();
            prob_collision();
            setEmax();
            //plotTable(); []   
        }

        //jadi kita declare beberapa Energy di sini, pada initializer list dimasukkan jumlah langkahnya. kayaknya nanti mending dijadikan satu sama parameter.txt yang lengkap
        //kalau dibikin rapi pake std::vector<Energy>, jadinya malah nggak tau urutan ke berapa itu proses apa, jadi mending dituliskan pelan-pelan aja, nggak perlu terlalu ringkas & rapi
        
        /*
        Energy energy_axis; //untuk sumbu x
        Energy v1; //tumbukan elastis
        Energy v2; //ionisasi
        */
        std::vector<Energy> proses; //mulai dari nol ya
        //inisiasinya harus dari luar nih.
        void addProcess(int jum_pros, int ni)
        {
            for (int j = 0 ; j < jum_pros; j++)
            {
                proses.emplace_back(ni);
            }
        }
        //coba panggil otomatis
        void readTable(std::string namanya);

        void plotTable()
        {
            std::ofstream fout("null_collision.dat");
            for (int i=0; i<ne; i++)
            {
                for (int j=0; j<jumlahProses; j++)
                {
                    fout << proses[j].get(i) << "\t";
                }
                fout << "\n";
            }
            fout.close();
        }

        //hitung frekuensi maksimal, kita perlu bikin container baru ding
        //Energy vmax;
        //hitung frekuensinya
        void frequency();
        //kemudian kita tentukan nilai maksimum untuk null collision
        void set_max_val()
        {
            for (int i=0; i<ne; i++)
            {
                maxi.emplace_back(proses[jumlahProses-1][i]);
            }
            max_val = *max_element(maxi.begin(), maxi.end());
        }
        //kemudian kita normalisasi
        void normalize();
        //hitung peluang terjadinya tumbukan
        void prob_collision();
        
        double getProbability()
        {
            return probability;
        }
        double EtoL(double E)
        {
            return E/stepEnergi; //ini karena energi selalu mulai dari nol
        }

        double max_val = 1;
        const int ne;
        std::string nama;
        const double toSI = 1e-20; //nggak perlu diset di constructor lah, ribet.
        void setEmax()
        {
            Emax = stepEnergi*(ne-1);
        }
        double getEmax()
        {
            return Emax; //ini dalam eV.
        }
        //alat bantu untuk mencari frekuensi maksimal
        std::vector<double> maxi;

        double get_probability() //untuk mengambil nilai peluang maksimal [200922]
        {
            return probability;
        }

        int get_jumlahProses() //untuk mengambil jumlah proses yang ditinjau [200922]
        {
            return jumlahProses;
        }

    protected:
        World &world;
        CollisionType collision_type;
        double mass;
        double stepEnergi;
        double Emax=0;
        double probability=0;

        int jumlahProses;
        
};

#endif