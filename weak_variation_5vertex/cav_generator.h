#ifndef CAV_GENERATOR_H
#define CAV_GENERATOR_H
#include<fstream>
#include<iostream>
#include<sstream>
#include<random>
#include"World.h"
#include"Field.h" //vec2 ada di sini 
#include<bits/stdc++.h>

//kita perlu juga instance untuk mengecilkan resolusinya
struct Point
{
    double x,y;
};

struct line
{
    Point p1, p2;
};

class cav_generator
{
    public:
        cav_generator(int n, World &world, Low_res &low_res) : n_titik{n}, world{world}, low_res{low_res} //(0,0) nggak dihitung 
        {
            titik = new double2[n];
            vertices = new Point[n+1]; //(0,0) dihitung
            clear();
        } 
        void clear()
        {
            for(int i=0; i<n_titik; i++)
            {
                titik[i] = 0;
                vertices[i+1].x = 0; vertices[i+1].y = 0; //vertex pertama selalu nol, aman lah
            }
        }
        void setTitikAwal(int i, double2 posisi); //ini mungkin diset beberapa kali ya
        void generate_out(int iterawal, int jumlah_output); //ini generate banyak output, file csv nya juga disatukan di sini
        void nama_output(int i); //ini untuk satu file aja
        void print_titik();
        void coarseVariation(); //variasi bentuk umum dari cavity
        void fineVariation(); //variasi bentuk khusus dari cavity
        void weakVariation();
        
        bool zona(line l1, Point p); 
        bool onLine(line l, Point p);
        int cross_prod(Point a, Point b, Point c);
        bool isIntersect(line l1, line l2);
        bool checkInside(Point poly[], int n, Point p);

        std::vector<double2> vertex;
        void addVertex(double2 posisi);
        int getn_titik()
        {
            return n_titik;
        }
        void nodeFlag(World &world); //untuk ngeflag node node yang ada di dalam poligon yang dibatasi titik-titik vertices
        void flagOutput(World &world, int i); //pengganti Output.h, bisa ganti ke file csv juga nanti kalau butuh
        
        void reduce_resolution(World &world, Low_res &low_res);
        void pixelOutput(Low_res &low_res, int i);

        //bikin file .bat untuk run banyak file sekaligus
        void batWrite(int iterawal, int jumlah_output); //sementara pakai model_i.af
        
        Point *vertices; //lebih nggak pusing pakai Point daripada pakai double2
    private:
        double2 *titik;
        int n_titik;
        World &world;
        Low_res &low_res;
};


#endif