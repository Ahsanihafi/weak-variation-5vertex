#include"Field.h"
#include"Crosssection.h"
#include<fstream>
#include<vector>
#include<algorithm>
#include<iostream>
#include<math.h>
#include<iomanip>
#include<sstream>

using namespace std;

int main()
{
    //panggil classnya
    //tumbukan elastis
    int nn = 25;
    double dn = 0.3;
    //elastis antara elektron dengan H2
    C_Section data_final(nn,dn); //sekali aja
    
    int itermax = 100;
    data_final.stringwrite(itermax);

    //buka output
    ofstream fout("output_csv/outputmedan.csv");     
    ifstream baca("parameter_cavity.dat");
    for (int i=1; i<itermax; i++)
    {
        double frekuensi, skala, max_medan;
        baca >> frekuensi >> skala >> max_medan;

        stringstream nama_input;
        nama_input << "output_dat/medanout_"<<i<<".dat";
        cout << nama_input.str() << "\n";
        data_final.ifTable(nama_input.str());
        data_final.setVal();
        for(int j=0; j<nn; j++)
        {
            if(j==nn-1) //paling akhir (kanan)
            {
                fout << data_final.Cross_S[j] << "," << frekuensi*1000 << "," << skala*100000; //<< "," << max_medan << "\n" ;
            }
            else
            {
                fout << data_final.Cross_S[j] << ",";
            }
        }
        if(i==itermax-1)
        {
            continue;
        }
        else
        {
            fout << "\n";
        }
        //clear Cross_S
        data_final.clearTable();
        data_final.Cross_S.clear();
    }
    baca.close();
    fout.close();
    return 0;
}
