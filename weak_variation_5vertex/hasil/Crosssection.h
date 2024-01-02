#ifndef CROSS_H
#define CROSS_H

#include<fstream>
#include<vector>
#include<string>
#include"Field.h"
#include<bits/stdc++.h>
//container untuk tabel
struct CrossTable
{
    //constructornya dulu
    CrossTable(double E, double sigma) : energVal{E}, crossSection{sigma} {}
    double energVal; //kolom pertama
    double crossSection; //kolom kedua
};

//ini cross section untuk tiap proses yang mungkin terjadi, jadi nanti ada sendiri-sendir untuk tiap proses (nggak dijadikan satu dulu)
//setelah dibaca dan disimpan dalam bentuk crosstable, mungkin kita perlu menyeragamkan interval. 
class C_Section
{
    public:
        //pertama terkait constructornya
        C_Section(const int langkah, double lebar_langkah) : ne{langkah}, Cross_S(langkah), stepEnergi{lebar_langkah} 
        {
            batasEnergi = (langkah-1)*lebar_langkah;
        }
        //data yang sudah dibaca disimpan di sini
        std::vector<CrossTable> crosstable;
        void clearTable()
        {
            crosstable.erase(crosstable.begin(),crosstable.end());
        }
        //kita butuh juga algoritma baca data, pastikan data berupa tabel dua kolom
        void ifTable(std::string nama);
        //untuk interpolasi data dengan argumen nilai energi
        double interpolate(double E);
        //Vektor di bawah merupakan luaran yang diharapkan dari class ini
        Energy Cross_S;
        //seperti biasa, kita butuh converter dari koordinat real ke koor
        double EtoL(double E) const
        {
            return E/stepEnergi;            
        }
        //selanjutnya kita butuh fungsi untuk menentukan nilai dari Cross_S, stepsize dalam satuan eV
        void setVal();

        //jadi nantinya gini, misal data mulai dari E0=10, semetara Cross_S[0] itu terkait E=0, maka dari interpolate di atas, otomatis yang diassign adalah nilai nol. Jadi di scatter(double E, double val), val yang dimasukkan adalah nol, nggak masalah ya. 
        const int ne;
        double getStep()
        {
            return stepEnergi;
        }

        std::string unitOmission(std::string s); //string dipecah di bagian space.
        std::string firstOmission(std::string s);
        std::string nthOmission(std::string s, int n);
        void stringwrite(int itermax);
        void nfinder(int &n, std::string namafile, std::string namabatas);
    protected:
        double stepEnergi; //lebar langkahnya
        double batasEnergi;

};

//Setelah seluruh file dibaca, file kemudian disatukan dalam bentuk sigma_T. 

/*
class Sigma_T
{
    public:
        Sigma_T(const int langkah, double lebar_langkah) : stepEnergi{lebar_langkah}, ne{langkah}, argon_elastis(ne,stepEnergi), argon_exc1(ne,stepEnergi), argon_exc2(ne,stepEnergi), v1(ne,stepEnergi), v2(ne,stepEnergi)
        {
            batasEnergi = (langkah-1)*lebar_langkah;
            
            //baca data
            argon_elastis.ifTable("cross.txt");
            //set nilai Cross_S
            argon_elastis.setVal();
            
            //baca data
            argon_exc1.ifTable("1_exc.txt");
            //set nilai Cross_S
            argon_exc1.setVal();
            
            //baca data
            argon_exc2.ifTable("2_exc.txt");
            //set nilai Cross_S
            argon_exc2.setVal();

            v1(argon_elastis.Cross_S);
                        
        }
        //pertama untuk tumbukan elastis
        C_Section argon_elastis;
        //eksitasi pertama
        C_Section argon_exc1;
        //eksitasi kedua
        C_Section argon_exc2;

        Energy v1;

        Energy v2;

        const int ne;
    protected:
        double stepEnergi;
        double batasEnergi;
};
*/
#endif