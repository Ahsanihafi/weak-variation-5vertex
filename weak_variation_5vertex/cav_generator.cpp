#include"cav_generator.h"
#include<sstream>
#include<fstream>
#include<iostream>
#include"World.h"
#include<bits/stdc++.h>
void cav_generator::setTitikAwal(int i, double2 posisi)
{
    if((i>n_titik-1) || (i<0))
    {
        std::cout << "nilai i=" << i << " di luar indeks" << std::endl;
    }
    else
    {
        titik[0][i] = posisi[0];
        titik[1][i] = posisi[1];
        //std::cout << i << "\t" << titik[0][i] << "\t" << titik[1][i] << "\t" << posisi[0] << "\t" << posisi[1] << "\n";
    }
    print_titik();
    //std::cout << "tanda \t" << titik[0][0] << "\t" << titik[1][2] << "\n";
}

void cav_generator::addVertex(double2 posisi)
{
    vertex.emplace_back(posisi);
}

//jadi di sini, vertex itu fix, vertices itu divariasikan
void cav_generator::coarseVariation() //pergeseran besar lima vertex awal
{
    //ada tiga variasi, radius cavity, panjang cavity, radius drift tube
    double rad_cav_var = rnd();
    if (rad_cav_var>0)
    {
        rad_cav_var *= 2;
    }
    double rad_cav = vertex[0][1] + rad_cav_var;
    double len_cav_var = rnd();
    if((len_cav_var>0))
    {
        len_cav_var *= 2;
    }
    double len_cav = vertex[2][0] + len_cav_var;
    double rad_tube_var = rnd();
    double rad_tube = vertex[4][1] + rad_tube_var;
    //asumsi ada 7 vertex. jadi nilainya vertex[n][i] itu fix seterusnya, sementara vertices[n] berubah-ubah. ini sementara gapapa manual lah.
    //1,2,3
    for(int i=1; i<4; i++){
        vertices[i].x = (i-1)*len_cav/2; //kelipatan 0.6
        vertices[i].y = rad_cav;//ini nggak perlu ada resize nya
    }  
    //                                                                                              
    for(int i=4; i<6; i++){
        vertices[i].x = len_cav ; //x nya sama semua
        vertices[i].y = rad_cav - (i-3)*(rad_cav-rad_tube)/2;
    }
    vertices[6].x = vertex[5][0];
    vertices[6].y = rad_tube;
    vertices[7].x = vertex[6][0];
    vertices[7].y = vertex[6][1];
}
/*
void cav_generator::fineVariation() //variasi kecil dari s, mungkin baiknya kita pakai vektor aja biar bagus
{
    //fine variation ini selalu dilakukan setelah coarse variation. perlu diketahui dulu jarak vertex-vertex yang baru. kemudian diupdate
    double rad_cav = vertices[1].y;
    double len_cav = vertices[6].x;
    double rad_tube = vertices[11].y;
    //tiga besaran di atas digunakan untuk menentukan seberapa besar variasi dilakukan. Ini untuk kasus 12 buah vertex
    //nilai di bawah ini bisa dianggap jadi panjang vektornya. 
    double var_rad_mag = (rad_cav-rad_tube)/5;
    double var_len_mag = (len_cav)/4;

    //variasikan nilai vertices.
    //vertices[1] cuma bisa ke atas atau bawah
    vertices[1].y += var_rad_mag*rnd();
    double sintheta ; //bisa positif bisa negatif
    double costheta ; //selalu positif
    Point vektordatar; 
    //sementara, dari dua sampai sepuluh, dibikin bener-bener acak
    for(int i=2; i<5; i++)
    {
        sintheta = rnd();
        costheta = sqrt(1-sintheta*sintheta);
        vertices[i].x = vertices[i-1].x+var_len_mag*costheta; 
        vertices[i].y = vertices[i-1].y+var_len_mag*sintheta;
    }
    //titik pojok kanan atas variasnya sendiri, sementara vertices 5 fix
    //vertices[5].x += var_len_mag*rnd();
    //vertices[5].y += var_rad_mag*rnd();
    //var rad mag harusnya dihitung setelah posisi vertices[5] sudah ditentukan. kita tes pakai 1000 data dulu
    //tapi gimana misalnya vertices[5].y < rad_tube? mungkin posisi vertices 5 nggak ditentukan lewat vektor ya.
    var_rad_mag = (vertices[5].y-rad_tube)/5;
    for(int i=6; i<10; i++)
    {
        sintheta = rnd();
        costheta = sqrt(1-sintheta*sintheta);
        vertices[i].x = vertices[i-1].x + var_rad_mag*sintheta; //dikecilin aja bagian ini
        vertices[i].y = vertices[i-1].y - var_rad_mag*costheta;
    }
}
*/
void cav_generator::weakVariation()
{
    //jadi mestinya variasi x dari vertices[6-9] tergantung juga sama radius cavity nya. Jadi harus diketahui radius defaultnya
    double rad_cav_var_default = vertex[0][1] - vertex[4][1];
    double len_cav_var_default = vertex[2][0];
    double rad_cav = vertices[1].y;
    double len_cav = vertices[3].x;
    double rad_tube = vertices[5].y;
    //kasih kriteria dulu,
    double var_rad_mag;
    double var_len_mag;
    if(len_cav>len_cav_var_default)
    {
        var_rad_mag = 0.2*(rad_cav-rad_tube);
    }
    else
    {
        var_rad_mag = 0.2*(rad_cav-rad_tube)*len_cav/len_cav_var_default;
    }
    if(rad_cav > vertex[0][1])
    {
        var_len_mag = 0.2*len_cav;
    }
    else
    {
        var_len_mag = 0.2*len_cav*rad_cav/vertex[0][1];
    }
    vertices[2].x += var_len_mag*rnd();
    vertices[2].y += var_rad_mag*rnd();

    vertices[4].x += var_len_mag*rnd();
    vertices[4].y += var_rad_mag*rnd();

}

void cav_generator::nama_output(int i)
{
    //variasikan langsung di sini
    coarseVariation();
    //fineVariation();
    weakVariation();
    std::stringstream nama;
    nama << "hasil/input_sfo/model_" << i << ".af";
    std::ofstream fout(nama.str());
    //ini awalan namanya
    fout << "2.4-GHz TM010 Modified Pillbox Cavity \n";
    fout << "Pillbox cavity with a 2-cm-radius bore tube \n";
    fout << "\n";
    //ini sementara untuk kasus simetri aksial
    fout << "&reg kprob=1, \n";
    fout << "icylin=1 \n";
    fout << "dx=.125, \n";
    fout << "freq=2000, \n"; //tebakan frekuensi ya
    //double yvardrive = vertex[0][1] + rnd();
    fout << "xdri=0,ydri="<< vertices[1].y <<", \n"; //cari gunanya drive point buat apa, kayaknya kita selalu asumsikan nilainya di titik pojok kiri. 
    fout << "kmethod=1 \n";
    fout << "beta=0.95 & \n";
    fout << "\n";
    //titik selalu diawali dengan nol dan diakhiri dengan nol juga. 
    fout << "$po x=0.0,y=0.0 $ \n";
    //titik pertama cuma divariasi posisi y nya, sementara titik terakhir divariasi x nya
    for (int n=0; n<n_titik; n++)
    {
        if (n==0)
        {
            //std::cout << n << "\t" << vertex[n][0] << "\t" << vertex[n][1] << "\n";
            //yang divariasikan cuma y, x nya konstan
            //double yvar = titik[1][n] + rnd(); //skema variasinya sementara ngawur dulu, dicek dulu programnya bisa jalan atau enggak [20072023]
            fout << "$po x=0.0,y="<< vertices[n+1].y  << "$ \n";
            //vertices[n+1].x=0; 
            //vertices[n+1].y=yvardrive;
        }
        else if(n==n_titik-2)
        {
            //std::cout << n << "\t" << vertex[n][0] << "\t" << vertex[n][1] << "\n";
            //yang divariasikan cuma x, y nya konstan
            //double xvar = vertex[n][0]; //+ rnd();
            //double yvar = vertex[n][1];
            //std::cout << n << "\t" << titik[0][n] << "\t" << xvar <<  "\n";
            fout << "$po x="<< vertices[n+1].x <<",y=" << vertices[n+1].y << " $ \n";
            //vertices[n+1].x=xvar;
            //vertices[n+1].y=yvar;
        }
        else if(n==n_titik-1)
        {
            //std::cout << n << "\t" << vertex[n][0] << "\t" << vertex[n][1] << "\n";
            //yang divariasikan cuma x, y nya konstan
            //double xvar = vertex[n][0]; //+ rnd();
            //std::cout << n << "\t" << titik[0][n] << "\t" << xvar <<  "\n";
            fout << "$po x="<< vertices[n+1].x <<",y=0.0 $ \n";
            //vertices[n+1].x=xvar;
            //vertices[n+1].y=0;
        }
        else
        {
            //std::cout << n << "\t" << vertex[n][0] << "\t" << vertex[n][1] << "\n";
            //double xvar = vertex[n][0] + 0.3*rnd();
            //double yvar = vertex[n][1] + 0.3*rnd();
            fout << "$po x="<< vertices[n+1].x <<",y="<< vertices[n+1].y <<" $ \n";
            //vertices[n+1].x=xvar;
            //vertices[n+1].y=yvar;
        }
    }
    fout << "$po x=0.0,y=0.0 $ \n";
}

void cav_generator::generate_out(int iterawal, int jumlah_output)
{
    std::ofstream fout("model.csv");
    for (int iter=iterawal; iter<jumlah_output+iterawal; iter++)
    {
        if(iter%100 == 0)
        {
            std::cout << "model ke=" << iter << "\n";
        }
        nama_output(iter); //file .af muncul di sini
        //flagging node juga dituliskan di sini.
        nodeFlag(world);
        flagOutput(world, iter);
        reduce_resolution(world, low_res);
        pixelOutput(low_res, iter);
        //kita perlu juga bikin output csv di sini
        for(int i=0; i<low_res.nz; i++)
        {
            for(int j=0; j<low_res.nr; j++)
            {
                if((i==low_res.nz-1)&&(j==low_res.nr-1))
                {
                    fout << low_res.pixel[i][j];
                }
                else
                {
                    fout << low_res.pixel[i][j] << ",";
                }
            }
        }
        if(iter==jumlah_output-1)
        {
            continue;
        }
        else
        {
            fout<<"\n";
        }
    }
    fout.close();
}

void cav_generator::print_titik()
{
    std::cout << "print titik \n";
    for (int i=0; i<vertex.size(); i++)
    {
        std::cout << i << "\t" << vertex[i][0] << "\t" << vertex[i][1] << "\n";
    }
    std::cout << "================================= \n";
}

//fungsi buat ngecek kalau suatu titik ada di zona persegi panjang garis (garis tsb adalah salah satu diagonal dari persegi)
bool cav_generator::zona(line l1, Point p)
{
    if((p.x <= std::max(l1.p1.x,l1.p2.x)) && (p.x >= std::min(l1.p1.x,l1.p2.x)) 
        && (p.y <= std::max(l1.p1.y,l1.p2.y)) && (p.y >= std::min(l1.p1.y,l1.p2.y)) )
    {
        return true;
    }
    else
    {
        return false;
    }
} //kalau di luar zona nggak mungkin bisa titik 

//untuk menunjukkan kalau titik p ada pada garis l
bool cav_generator::onLine(line l, Point p)
{
    if(zona(l,p) == true) 
    {
        double val = (l.p2.x-l.p1.x)*(p.y-l.p1.y) - (l.p2.y-l.p1.y)*(p.x-l.p1.x);
        if(val == 0)// cross productnya nilainya nol
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

//perhitungan cross product dari dua vektor yang tersusun dari tiga titik (misal line 1 terdiri atas titik a dan b)
int cav_generator::cross_prod(Point a, Point b, Point c) //kita sekedar mau tau (a,b,c) itu colinear, mutar kanan, atau mutar kiri
{
    //ini vektornya b-a dan c-b, agak mbingungi ya
    double val = (b.y-a.y)*(c.x-b.x) - (c.y-b.y)*(b.x-a.x);
    if (val==0) //sejajar
    {
        return 0;
    }
    else if (val < 0) //berlawanan arah jarum jam
    {
        return 2;
    }
    else    //searah jarum jam
    {
        return 1;
    }
}

//kemudian andaikan ada dua garis, perlu dicek juga apakah saling bersilangan atau enggak. 
bool cav_generator::isIntersect(line l1, line l2)
{
    //paham maksudnya gimana. vektor yang terbentuk harus "tengok kanan kiri"
    int dir1 = cross_prod(l1.p1,l1.p2,l2.p1);
    int dir2 = cross_prod(l1.p1,l1.p2,l2.p2);
    int dir3 = cross_prod(l2.p1,l2.p2,l1.p1);
    int dir4 = cross_prod(l2.p1,l2.p2,l1.p2);
    if ((dir1 != dir2) && (dir3 != dir4)) 
    {
        return true;
    }
    else if((dir1 == 0) && (onLine(l1,l2.p1))) //ketika p1 dari line2 ada pada line1
    {
        return true;
    }
    else if((dir2 == 0) && (onLine(l1,l2.p2))) //ketika p2 dari line2 ada pada line1
    {
        return true;
    }
    else if((dir3 == 0) && (onLine(l2,l1.p1))) //ketika p1 dari line1 ada pada line2
    {
        return true;
    }
    else if((dir4 == 0) && (onLine(l2,l1.p2))) //ketika p2 dari line1 ada pada line2
    {
        return true;
    }
    else
    {
        return false;
    }
}

//ini perlu dicek juga 
bool cav_generator::checkInside(Point poly[], int n, Point p) //poly[] menunjukkan himpunan titik-titik vertex, n menunjukkan jumlah vertex
{
    if (n<3) //ini bukan polygon
    {
        return false ;
    }
    //triknya, buat garis horizontal dengan satu ujungnya di p, ujung lainnya di titik yang jauh
    line exline = {p, {9999, p.y+0.1}};
    int count = 0; //ini untuk menghitung berapa kali garis exline memotong garis sisi poligon
    int i=0; //indikator sisi poligon (jadi untuk suatu poligon, titik hanya terhubung ke titik terdekatnya saja)
    do
    {
        //bikin garis antara dua titik poligon yang berurutan
        line side = {poly[i],poly[(i+1)%n]}; //dikasih modulus supaya titik terakhir balik lagi ke titik pertama
        //jika garis berpotongan
        if (isIntersect(side, exline))
        {
            //misal sidenya juga mendatar, dan titik ada di side tsb, nggak perlu hitung yang lain
            if(cross_prod(side.p1,p,side.p2)==0) //colinear atau enggak. ini sebaiknya kita bikin hampir tidak mungkin terjadi. dengan memiringkan exline sedikit
            {
                return onLine(side,p); //bagian ini aku nggak paham
            }
            count ++;
        }
        i = (i+1)%n; //kalau i+1 = n, maka balik lagi nilainya ke nol.
    } while (i!=0); //berarti di iterasi pertama langsung digeser nilai i ya, sampai berulang (pakai modulus)
    return count & 1;
}

//untuk output kita bikin sendiri di sini, nggak mengikuti Output.h
void cav_generator::flagOutput(World &world, int i)
{
    //membuka file output, nama di atas digunakan di sini
    std::stringstream nama;
    nama << "hasil/hasil_vti/flag_" << i << ".vti";
    std::ofstream out(nama.str());
    if(!out.is_open())
    {
        std::cerr << "File flag tidak bisa dibuka" << std::endl;
        return;
    }

    //output dituliskan dalam format vtk dengan mesh kartesian terstruktur
    out<<"<VTKFile type=\"ImageData\">\n";
	double2 x0 = world.getX0();
	double2 dh = world.getDh();
	out<<"<ImageData Origin=\""<<x0[0]<<" "<<x0[1]<<" "<<0<<"\" "; //Karena  
	out<<"Spacing=\""<<dh[0]<<" "<<dh[1]<<" "<<1<<"\" "; //spaceing diset satu
	out<<"WholeExtent=\"0 "<<world.nz-1<<" 0 "<<world.nr-1<<" 0 "<<0<<"\">\n"; //di sini entri ketiga selalu nol (sudut azimuth), cuma ada satu entri

    //output disimpan dalam bentuk node (point data)
    out<<"<PointData>\n";

    //output medan pertama, object_id berupa medan skalar
    out<<"<DataArray Name=\"object_id\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
	out<<world.object_id;
	out<<"</DataArray>\n";

    //penutup
    out<<"</PointData>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();
}

void cav_generator::nodeFlag(World &world)
{
    for(int i=0; i<world.nz;i++)
    {
        for(int j=0; j<world.nr; j++)
        {
            double posx = i*world.getDh()[0];
            double posy = j*world.getDh()[1];
            Point p = {posx,posy};
            if(checkInside(vertices, n_titik+1, p))
            {
                //cout << "Titik ada di dalam dengan posisi x=" << p.x << " y=" << p.y << "\n";
                world.object_id[i][j] = 1;
            }
            else
            {
                //cout << "Titik ada di luar dengan posisi x=" << p.x << " y=" << p.y << "\n";
                world.object_id[i][j] = 0;
            }
        }
    }
}

void cav_generator::reduce_resolution(World &world, Low_res &low_res)
{
    //setelah selesai pakai nodeflag
    low_res.pixel.clear(); //kosongin dulu
    for(int i=0; i<world.nz; i++)
    {
        for(int j=0; j<world.nr; j++)
        {
            double nilai = world.object_id[i][j]; //ini nilainya sesudah diupdate yak
            double posx = i*world.getDh()[0];
            double posy = j*world.getDh()[1];
            //std::cout << "posx = " << posx << "\t posy = " << posy << "\n"; 
            //operasi XtoL di low_res.
            double2 lc = low_res.XtoL({posx,posy});
            //std::cout << "n=" << i*world.nz+j << "\t lc1=" << lc[0] << "\t lc2=" << lc[1] << "\n";
            low_res.pixel.scatterCartesian(lc, nilai);
        }
    }
}

void cav_generator::pixelOutput(Low_res &low_res, int i)
{
    //membuka file output, nama di atas digunakan di sini
    std::stringstream nama;
    nama << "hasil/hasil_vti/pixel_" << i << ".vti";
    std::ofstream out(nama.str());
    if(!out.is_open())
    {
        std::cerr << "File flag tidak bisa dibuka" << std::endl;
        return;
    }

    //output dituliskan dalam format vtk dengan mesh kartesian terstruktur
    out<<"<VTKFile type=\"ImageData\">\n";
	double2 x0 = low_res.getX0();
	double2 dh = low_res.getDh();
	out<<"<ImageData Origin=\""<<x0[0]<<" "<<x0[1]<<" "<<0<<"\" "; //Karena  
	out<<"Spacing=\""<<dh[0]<<" "<<dh[1]<<" "<<1<<"\" "; //spaceing diset satu
	out<<"WholeExtent=\"0 "<<low_res.nz-1<<" 0 "<<low_res.nr-1<<" 0 "<<0<<"\">\n"; //di sini entri ketiga selalu nol (sudut azimuth), cuma ada satu entri

    //output disimpan dalam bentuk node (point data)
    out<<"<PointData>\n";

    //output medan pertama, object_id berupa medan skalar
    out<<"<DataArray Name=\"object_id\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<low_res.pixel;
	out<<"</DataArray>\n";

    //penutup
    out<<"</PointData>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();
}

void cav_generator::batWrite(int iterawal, int jumlah_output)
{
    std::ofstream fout("run_model.bat"); 
    int itermax = iterawal + jumlah_output;
    for(int i=iterawal; i<itermax; i++)
    {
        fout << "start /W \" \"  \"%SFDIR%autofish\" model_" << i << "\n";
    }
    fout.close();
}