#include"World.h" //buat random generator
#include"Output.h"
#include"cav_generator.h"
using namespace std;



int main()
{
    World dunia(1001,1001);
    Low_res miniatur(21,21);
    dunia.setextents({0,0},{9.0,9.0});
    miniatur.setextents({0,0},{10.0,10.0});
    //kita coba flagnya pakai objectId aja ya. 
    int jumlah_model = 10;
    int nomor_awal_model = 0;
    cav_generator pillbox(7,dunia, miniatur); //nggak termasuk (0,0)
    pillbox.batWrite(nomor_awal_model, jumlah_model);
    /*
    pillbox.addVertex({0.0,5.0}); //1 - pojok fix x
    pillbox.addVertex({0.6,5.0}); //2
    pillbox.addVertex({1.2,5.0}); //3
    pillbox.addVertex({1.8,5.0}); //4
    pillbox.addVertex({2.5,5.0}); //5 - pojok
    pillbox.addVertex({2.5,4.4}); //6
    pillbox.addVertex({2.5,3.8}); //7
    pillbox.addVertex({2.5,3.2}); //8
    pillbox.addVertex({2.5,2.6}); //9
    pillbox.addVertex({2.5,2.0}); //10 - pojok
    pillbox.addVertex({7.5,2.0}); //11 - pojok fix x
    pillbox.addVertex({7.5,0.0}); //12 - pojok fix
    */
    //di bawah ini untuk kasus 7 titik
    pillbox.addVertex({0.0,5.0}); //1 pojok
    pillbox.addVertex({1.25,5.0}); //2
    pillbox.addVertex({2.5,5.0}); //3 pojok
    pillbox.addVertex({2.5,3.5}); //4
    pillbox.addVertex({2.5,2.0}); //5 pojok
    pillbox.addVertex({7.5,2.0}); //6 
    pillbox.addVertex({7.5,0.0}); //7 pojok

    pillbox.print_titik();
    pillbox.generate_out(nomor_awal_model, jumlah_model); // di sini udah langsung ada check inside dan flag outputnya juga
    
    return 0;
}
