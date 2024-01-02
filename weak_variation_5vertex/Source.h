#ifndef SOURCE_H
#define SOURCE_H

#include"World.h"
#include"Species.h"
#include<string>
#include<fstream>

class ColdBeamSource
{
	public:
		//constructor beserta dengan initializer list
		ColdBeamSource(Species &species, World &world, double v_drift, double arus) : spec{species}, world{world}, v_drift{v_drift}, arus{arus}
		{ }

		void sample();
		//ini untuk sampling merata dengan kecepatan awal nol, jangan lupa pembobotannya juga, seperti pada sample()
		void sample2(int num_sim); //dipanggil secara acak, tidak terdistribusi secara uniform

		void sampleThermalCyclotronic(int num_sim, double Ek); //jadi kecepatan acak yang disampel akan tergantung pada energi yang diinginkan

		void sampleThermalBoris(int num_sim, double Ek); //jadi kecepatan acak yang disampel akan tergantung pada energi yang diinginkan

		double getArus()
		{
			return arus;
		}
		double getVdrift()
		{
			return v_drift;
		}
		double getRadius()
		{
			return radius;
		}
		Species getSpec() //diperlukan terutama nanti di paramOutput
		{
			return spec;
		}
		void sampleIonHPlus();

		void sampleIonH2Plus();

		void sampleIonHMin();

		void sampleElecIonCyc(); //untuk elektron hasil ionisasi [151022]

		void sampleElecIonBoris();

		void sampleElec(double Ek);

		void sampleElec2(double Ek);

		void sampleElecWallCyc(); //Ek terkait dengan suhu dari dinding

		void sampleElecWallBoris();

		void sampleBacaBoris(std::string namafile); //ini baca dari file, untuk penggerak boris

		void sampleBacaPartialBoris(std::string namafile, double fraction);

		void sampleBacaCyclotronic(std::string namafile); //ini baca dari file, untuk penggerak cyclotronic.

		void sampleBacaPartialCyclotronic(std::string namafile, double fraction);

	protected:
		World &world;
		Species &spec;
		double v_drift; //kecepatan drift rata-rata dari berkas partikel 
		double arus; //arus berkasnya berapa, jadi agak beda dengan yang digunakan sama Brieda (dia pakai density)
		double radius = 2e-2; 
		int jumlahmakro = 0;
};


#endif
