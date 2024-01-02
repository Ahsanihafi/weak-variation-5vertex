#include<fstream>
#include<vector>
#include<string>
#include"Field.h"
#include"Crosssection.h"
#include<math.h>
#include"World.h"
#include<iostream>
void CrossSection::readTable(std::string namanya)
{
    std::ifstream baca(namanya.c_str());
    for (int i=0; i<ne; i++)
    {
        double *a;
        a = new double [jumlahProses];
        //double a,b,c;
        //baca >> energy_axis[i] >> v1[i] >> v2[i] >> v3[i] >> v4[i];
        for (int j=0; j < jumlahProses; j++)
        {
            baca >> a[j];
        }
        //baca >> a >> b >> c ;
        //energy_axis.setValue(i,a);
        //v1.setValue(i,b);
        //v2.setValue(i,c);

        for (int j=0; j < jumlahProses; j++)
        {
            proses[j].setValue(i,a[j]);
        }

        delete[] a; //bener gini nggak ya? diisi, dibuang, diisi, dibuang
    }
}

//JANGAN LUPA KALAU PROSES[0] ITU SUMBU ENERGI, BUKAN TAMPANG LINTANG

//untuk frekuensi, nilainya diubah ke SI dulu, supaya konsisten. 
void CrossSection::frequency()
{
    double n = world.getNgas(); //ini n acuan dari tekanan. n yang dipakai di bawah nggak harus sama, diset di sini sekalian juga. 
    double nn[jumlahProses];
    /*
    for (int j=0; j<jumlahProses; j++)
    {
        if (j<4)
        {
            nn[j] = n;
        }
        else
        {
            nn[j] = 0.1*n*exp(-j+4); //ini baru estimasi, nanti perlu diganti yang lebih tepat lagi
        }
    }
    */
    for (int i=0; i<ne; i++)
    {
        //ini kalau prosesnya lebih banyak, v4 harus diganti lho, hati-hati
        //kerjakan dalam SI semua, jadi ubah ke joule dulu
        //ini untuk kasus Dissociative electron attachment, nilai dari n harus divariasikan terhadap mode vibrasi (datanya ada di paper)
        //v2.multValue(i,n*toSI*sqrt(2*energy_axis.get(i)*abs(Const::QE)/mass));
        //v1.multValue(i,n*toSI*sqrt(2*energy_axis.get(i)*abs(Const::QE)/mass));

        for (int j = 1; j<jumlahProses; j++) //kalo energy axis nggak diapa-apain, pakai ini
        {
            proses[j].multValue(i,n*toSI*sqrt(2*proses[0].get(i)*abs(Const::QE)/mass));
        }
    }
}

void CrossSection::normalize()
{
    for (int i = 0; i<ne; i++)
    {
        double a = 1.0/max_val;
        //std::cout << max_val << "\t" << a << "\n";
        //v2.multValue(i,a);
        //v1.multValue(i,a);
        for (int j=1; j<jumlahProses; j++)
        {
            proses[j].multValue(i,a);
        }
    }
}

void CrossSection::prob_collision()
{
    probability = 1.0 - exp(-max_val*world.getDt());
}