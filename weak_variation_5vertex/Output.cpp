#include "Output.h"
#include "World.h"
#include "Species.h"
#include "Source.h"
#include "Potentialsolver.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

//bedanya dengan kasus 3 dimensi, untuk x[2] diset nol terus
void Output::fields(World &world, vector<Species> &species)
{
    //untuk memberikan nama file sesuai langkah waktu, tapi belum ditulis
    stringstream name;
    name << "results/fields_" << setfill('0') << setw(5) << world.getTs() << ".vti";

    //membuka file output, nama di atas digunakan di sini
    ofstream out(name.str());
    if(!out.is_open())
    {
        cerr << "File " << name.str() << " tidak bisa dibuka" << endl;
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

    //selanjutnya, medan volume node, skalar
    out<<"<DataArray Name=\"NodeVol\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.node_vol;
	out<<"</DataArray>\n";

    //potensial listrik, medan skalar
    out<<"<DataArray Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.phi;
	out<<"</DataArray>\n";

    //rapat muatan, medan skalar
    out<<"<DataArray Name=\"rho\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.rho;
	out<<"</DataArray>\n";

    //rapat bilangan untuk tiap-tiap spesies
    for (Species &sp:species)
	{
		out<<"<DataArray Name=\"nd."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den;
		out<<"</DataArray>\n";
	}

    /*ini kalau nanti fungsi pererata sudah dipakai
    for (Species &sp:species)
	{
		out<<"<DataArray Name=\"nd-ave."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den_ave;
		out<<"</DataArray>\n";
	}
    */

    //medan vektor, tiga komponen (tetap matriks rank 2 sih)
    out<<"<DataArray Name=\"ef\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.ef;
	out<<"</DataArray>\n";

    //penutup
    out<<"</PointData>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();
}

void Output::flagOutput(World &world)
{
	//membuka file output, nama di atas digunakan di sini
    ofstream out("flag.vti");
    if(!out.is_open())
    {
        cerr << "File flag tidak bisa dibuka" << endl;
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

//selanjutnya terkait output di terminal
void Output::screenOutput(World &world, vector<Species> &species)
{
	cout << "langkah: " << world.getTs();
	for (Species &sp:species)
	{
		cout << setprecision(3) << "\t " << sp.name << ":" <<sp.getNp() ;//<< "\t" << "Jumlah partikel:" << sp.getRealCount();
	}
	cout << endl;
}

//diagnostic yang ditaruh di file csv sementara nggak dipakai dulu. 

//untuk plot scatter partikel
/*
void Output::particles(World &world, vector<Species> &species, int num_parts)
{
    //jadi ada luaran file untuk masing-masing spesies dalam plasma
    for (Species &sp:species)
    {
        stringstream name; 
        name << "results/parts_" << sp.name << "_" << setfill('0') << setw(5) << world.getTs() << ".vtp"; //ekstensinya bukan vti ya]

        //buka file output
        ofstream out(name.str());
        if (!out.is_open())
        {
            cerr << "File " << name.str() << " tidak dapat dibuka" << endl;
            return;
        }

        //buat daftar partikel untuk output
        double dp = num_parts/(double)sp.getNp();
        double counter = 0;
        vector<Particle*> to_output;
        for (Particle &part:sp.particles)
        {
            counter +=dp;
            if (counter > 1)
            {
                to_output.emplace_back(&part);
                counter = 0;
            }
        }
        
        //header untuk filenya
        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out << "<PolyData>\n";
        out << "<Piece NumberOfPoints=\"" << to_output.size() << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
        out << "NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";

        //titik-titik yang diplot (scatter)
        out << "<Points>\n";
        out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (Particle *part:to_output)
        {
            out << part->pos <<"\n";
        }
        out << "</DataArray>\n";
        out << "</Points>\n";

        //selanjutnya, untuk kecepatan yang diplot
        out << "<PointData>\n";
        out << "<DataArray Name=\"vel." << sp.name << "\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (Particle *part:to_output)
        {
            out << part->vel << "\n";
        }
        out << "</DataArray>\n";
        out << "</PointData>\n";

		//kecepatan partikel
		out<<"<PointData>\n";
		out<<"<DataArray Name=\"vel."<<sp.name<<"\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		for (Particle *part: to_output)
		{
			out<<part->vel<<"\n";
		}
		out<<"</DataArray>\n";
		out<<"</PointData>\n";
		
        //penutup
        out << "</Piece>\n";
        out << "</PolyData>\n";
        out << "</VTKFile>\n";

        out.close();
    }
}
*/


void Output::particles(World &world, vector<Species> &species, int num_parts) {
	/*loop over all species*/

	for (Species &sp:species) {

		//open a phase_sp_it.vtp
		stringstream name;
		name<<"results/parts_"<<sp.name<<"_"<<setfill('0')<<setw(5)<<world.getTs()<<".vtp";

		/*open output file*/
		ofstream out(name.str());
		if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}

		/*build a list of particles to output*/
		/*
		double dp = num_parts/(double)sp.getNp();
		double counter = 0;
		vector<Particle*> to_output;
		for (Particle &part : sp.particles)	{
			counter+=dp;
			if (counter>1) //save particle
				{to_output.emplace_back(&part);counter=0;}
		}
		*/
		
		//mungkin kalau maunya satu partikel satu titik, bisa pakai algoritma berikut ini
		vector<Particle*> to_output;
		for (Particle &part:sp.particles)
		{
			to_output.emplace_back(&part);
		}

		/*header*/
		out<<"<?xml version=\"1.0\"?>\n";
		out<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
		out<<"<PolyData>\n";
		out<<"<Piece NumberOfPoints=\""<<to_output.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
		out<<"NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";

		/*points*/
		out<<"<Points>\n";
		out<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		for (Particle *part: to_output)
			out<<part->pos<<"\n";
		out<<"</DataArray>\n";
		out<<"</Points>\n";

		/*velocities*/
		out<<"<PointData>\n";
		out<<"<DataArray Name=\"vel."<<sp.name<<"\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		for (Particle *part: to_output)
			out<<part->vel<<"\n";
		out<<"</DataArray>\n";
		out<<"</PointData>\n";

		out<<"</Piece>\n";
		out<<"</PolyData>\n";
		out<<"</VTKFile>\n";

		out.close();
	}
}

void Output::particles2(World &world, Species &sp) {
	/*loop over all species*/
    //open a phase_sp_it.vtp
	stringstream name;
	name<<"results/parts_"<<sp.name<<"_"<<setfill('0')<<setw(5)<<world.getTs()<<".vtp";

	/*open output file*/
	ofstream out(name.str());
	if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}

		/*build a list of particles to output*/
		/*
		double dp = num_parts/(double)sp.getNp();
		double counter = 0;
		vector<Particle*> to_output;
		for (Particle &part : sp.particles)	{
			counter+=dp;
			if (counter>1) //save particle
				{to_output.emplace_back(&part);counter=0;}
		}
		*/
		
		//mungkin kalau maunya satu partikel satu titik, bisa pakai algoritma berikut ini
	vector<Particle*> to_output;
	for (Particle &part:sp.particles)
	{
		to_output.emplace_back(&part);
	}

	/*header*/
	out<<"<?xml version=\"1.0\"?>\n";
	out<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out<<"<PolyData>\n";
	out<<"<Piece NumberOfPoints=\""<<to_output.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
	out<<"NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";

	/*points*/
	out<<"<Points>\n";
	out<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (Particle *part: to_output)
		out<<part->pos<<"\n";
	out<<"</DataArray>\n";
	out<<"</Points>\n";

	/*velocities*/
	out<<"<PointData>\n";
	out<<"<DataArray Name=\"vel."<<sp.name<<"\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (Particle *part: to_output)
		out<<part->vel<<"\n";
	out<<"</DataArray>\n";
	out<<"</PointData>\n";

	out<<"</Piece>\n";
	out<<"</PolyData>\n";
	out<<"</VTKFile>\n";

	out.close();
	

		
	
}

void Output::paramOutput(World &world, vector<Species> &species, vector<ColdBeamSource> &coldbeamsource, CrossSection &cross_section, Potentialsolver &solver)
{
    ofstream fout("ruang_fase/parameter_simulasi.txt");
	fout <<setprecision(20);
    fout << "Domain Part (2D) -------------------------------------------------------- \n";
    fout << "zmin= " << "\t" << world.getX0()[0] << "\n";
    fout << "rmin= " << "\t" << world.getX0()[1] << "\n";
    fout << "zmax= " << "\t" << world.getXm()[0] << "\n";
    fout << "rmax= " << "\t" << world.getXm()[1] << "\n";
    fout << "nz= " << "\t" << world.nz << "\n";
    fout << "nr= " << "\t" << world.nr << "\n";
    fout << "Parameter simulasi-------------------------------------" << "\n";
    fout << "dt= " << "\t" <<  world.getDt() << "\n";
    fout << "maxiter= " << "\t" << world.getMaxiter() << "\n";
    fout << "tegangan= \t" << solver.getVoltage() << "\n";
	fout << "hambatan resistor= \t" << solver.getResistance() << "\n";
    fout << "magnet= \t" << world.getMagStrength() << "\n";
	fout << "tekanan gas= \t" << world.getPresGas() << "\n";
	fout << "suhu gas= \t" << world.getTempGas() << "\n";
	fout << "koefisien dinding= \t" << world.getGamma() << "\n";
    fout << "Parameter Spesies-------------------------------------" << "\n";
    for (Species &sp:species)
    {
        fout << sp.name <<  " +++++++++++++++++"  << "\n";
        fout << "Mass= " << "\t" << sp.mass << "\n";
        fout << "Charge= " << "\t" << sp.charge << "\n";
        fout << "Macroparticle Weight= " << "\t" << sp.mpw0 << "\n";
    }/*
    fout << "Parameter Source--------------------------------------" << "\n";
    for (ColdBeamSource &cold:coldbeamsource)
    {
        fout << cold.getSpec().name << "+++++++++++++++" << "\n"; 
        fout << "Drift Velocity= " << "\t" << cold.getVdrift()  << "\n";
        fout << "Current= " << "\t" << cold.getArus() << "\n";
		fout << "Beam radius= " << "\t" << cold.getRadius() << "\n";
    }
	*/
	fout << "Monte-Carlo part--------------------------------------------" << "\n"; //ini kayaknya buat elektron dengan H+H2 dulu, yang lain menyusul [200922]
	fout << "Jumlah proses = " << cross_section.get_jumlahProses() << "\n";
	fout << "Frekuensi maksimum = " << cross_section.max_val << "\n";
	fout << "Peluang maksimum = " << cross_section.get_probability() << "\n";

    fout.close();
}

//ini handle perlu supaya file tetep kebuka mungkin ya?
namespace Output
{
	std::ofstream f_diag;
}

void Output::diagOutput(World &world, vector<Species> &species, int langkah_awal)
{
	using namespace Output;
	stringstream nama;
	nama << "Particle_count_" << langkah_awal << ".csv";
	//pastikan dulu filenya kebuka, kalo kebuka nggak usah ditulis ini. karena ini row pertama.
	if (!f_diag.is_open())
	{
		f_diag.open(nama.str()); //sementara coba pakai csv dulu
		f_diag << "Langkah";
		for (Species &sp:species)
		{
			f_diag << "," << sp.name;
		}
		f_diag << endl;
		
	}
	f_diag << world.getTs();
	for (Species &sp:species)
	{
			f_diag << "," << sp.getNp() ; 
	}
	f_diag << "\n";
	if (world.getTs()%1000 == 0)
	{
		f_diag.flush();
	}
	
}