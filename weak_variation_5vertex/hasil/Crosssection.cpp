#include<fstream>
#include<vector>
#include<string>
#include"Field.h"
#include"Crosssection.h"
#include<iostream>
#include<iomanip>
using namespace std;

void C_Section::ifTable(string nama)
{
    ifstream baca(nama.c_str());
    while(!baca.eof())
    {
        double E;
        double sigmaE;
        baca >> E >> sigmaE;
        crosstable.emplace_back(E,sigmaE);
    }
}

//dari interpolasi linear ini, harapannya kita bisa assign nilai pada vektor dengan step yang seragam (di bawah ini stepnya nggak seragam, makanya nggak bisa pakai EtoL)
double C_Section::interpolate(double E)
{
    size_t Ne = crosstable.size();
    if (E<crosstable[0].energVal)
    {
        return 0;
    }        
    else if (E>crosstable[Ne-1].energVal)
    {
        return 0;
    }
    else
    {
        double sigmaE = 0;
        //sementara scan semua, berhubung nggak banyak datanya, kayaknya nggak terlalu masalah
        for (int i=0; i<Ne; i++)
        {   
            if( (E>=crosstable[i].energVal) && (E<=crosstable[i+1].energVal) )
            {
                double sigmae = ((crosstable[i+1].crossSection - crosstable[i].crossSection)*(E-crosstable[i].energVal)/(crosstable[i+1].energVal-crosstable[i].energVal)) + crosstable[i].crossSection;
                sigmaE = sigmae;
                break;
            }
            else
            {
                continue;
                if (i==Ne-1)
                {
                    cout << "tidak terjadi break \n"; 
                }
            }
        }
        return sigmaE;
    }
}

void C_Section::setVal()
{
    double E0 = 0;
    for (int i=0; i<ne; i++)
    {
        double elogic = EtoL(E0); //ini berhubung energi selalu kelipatan stepEnergi, berarti elogic selalu bulat ya?
        Cross_S.scatter(elogic,interpolate(E0)); //scatternya bilangan bulat terus, jadi nggak ada titik yang ndobel.
        E0 += stepEnergi;
    }
}

string C_Section::unitOmission(string s)
{
    stringstream ss(s);
    string word;
    while (ss >> word)
    {
        //cout << word << endl;
        break; //ini cuma ambil string pertama doang.
    }
    return word;
}

string C_Section::firstOmission(string s)
{
    stringstream ss(s);
    string word;
    int i=0;
    while (ss >> word)
    {
        if (i==1)
        {
            //cout << word << "\n";
        }
        i++;
    }
    return word;
}

string C_Section::nthOmission(string s, int n)
{
    stringstream ss(s);
    string word;
    int i=0; 
    while(ss >> word)
    {
        if (i==n)
        {
            break;
        }
        i++;
    }
    return word;
}

void C_Section::stringwrite(int itermax){
    vector<double> medan; vector<double> posisi; ofstream gout("parameter_cavity.dat"); gout << setprecision(12); double scaler = 0; double freq = 0;
    for(int iter=0; iter<itermax; iter++){
        if(iter%100==0){
            cout << "iterasi ke=" << iter<< "\n";}
        string kalimat; //ditaruh di sini atau di dalam loop ya?
        string beforeEqual = ""; string afterEqual = ""; stringstream nama; //stringstream ini buat baca. file output cuma ada satu aja
        stringstream nama_output; nama << "input_sfo/MODEL_" << iter << ".SFO" ; //disesuaikan nanti isinya apa
        nama_output << "output_dat/medanout_" << iter << ".dat"; ofstream fout(nama_output.str()); fout << setprecision(12); ifstream baca(nama.str());
        int n=1; int nfile = 0; int nakhir = 0; int nfreq = 0;
        nfinder(nfile, nama.str(), "Ez(V/m)"); nfinder(nakhir, nama.str(), "Total cavity stored energy (from program Fish)"); nfinder(nfreq, nama.str(), "RF cavity resonant frequency");
        while (getline(baca, kalimat)) {
            if(n==nfreq){
                double frekuensi = stod(nthOmission(kalimat, 1)); gout << frekuensi << "\t";}
            if((n>nfile) && (n<(nakhir))){
                double medanz = stod(firstOmission(kalimat)); double posz = stod(unitOmission(kalimat)); medan.emplace_back(medanz); posisi.emplace_back(posz);
                if (n==nakhir-1) {
                    getline(baca, beforeEqual, ':'); getline(baca, afterEqual, '\n'); scaler = 1/sqrt(stod(unitOmission(afterEqual)));
                    gout << scaler << "\t";}}
            n++;}
        baca.close();
        double medanmax = *max_element(medan.begin(), medan.end()); gout << medanmax << "\n";
        for (int i=0; i<medan.size(); i++){
            if (i<medan.size()-1){
                fout << posisi[i] << "\t" << medan[i]*scaler << "\n";}
            else{
                fout << posisi[i] << "\t" << medan[i]*scaler ;}}
        medan.erase(medan.begin(), medan.end()); posisi.erase(posisi.begin(),posisi.end()); fout.close();} gout.close();}

void C_Section::nfinder(int &n, string namafile, string namabatas)
{
    int i=1;
    ifstream file(namafile);
    string line = "";
    while (getline(file, line))
    {
        if(line.find(namabatas) != string::npos)
        {
            //cout << "ditemukan di line=" << i << "\n";
            n=i;
            break;
        }
        else
        {
            i++;
        }
    }
    file.close();
}
