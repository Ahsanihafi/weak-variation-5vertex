//ini definisikan vektor dulu ya, pakai struct

#ifndef FIELD_H
#define FIELD_H

#include<vector>
#include<ostream>
#include<assert.h>
#include<iostream>

template<typename T>
struct vec3
{
    //pertama, constructor dengan initializer list
    vec3 (const T u, const T v, const T w) : d{u,v,w} { }
    vec3 (const T a[3]) : d{a[0],a[1],a[2]} { }
    vec3 () : d{0,0,0} {}
    //tiga constructor di atas sebenernya sama aja sih, cuma beda argumen
    //sama yang terakhir itu kalo nggak dispesifikasi dianggap vektor null

    //overload operator akses dan assignment. aksesnya ada dua
    T& operator[] (int i) {return d[i];}
    T operator() (int i) const {return d[i];}
    vec3<T>& operator=(double s) 
    {
        d[0] = s; d[1] = s; d[2] = s;
        return (*this);
    }
    //overload penjumlahan dan pengurangan
    vec3<T>& operator += (vec3<T> other)
    {
        for (int i=0; i<3; i++)
        {
            d[i] = d[i] + other[i];
        }
        return (*this);
    }
    vec3<T>& operator -= (vec3<T> other)
    {
        for (int i=0; i<3; i++)
        {
            d[i] = d[i] + other[i];
        }
        return (*this);
    }
    /*coba overload =, bisa enggak ya
    vec3<T>& operator = (vec3<T> other)
    {
        for (int i=0; i<3; i++)
        {
            d[i] = other[i];
        }
        return (*this); //ini itu perlu kalo objeknya semisal double& gitu ya?
    }
    */
protected:
    T d[3];
};

//operasi antara vektor dengan skalar
//perkalian dengan skalar (ngelengkapinya belakangan lah, toh belum dipake)
template<typename T>
vec3<T> operator* (const vec3<T> &a, T s)
{
    return vec3<T>(a(0)*s, a(1)*s, a(2)*s);
}
//perkalian dari kanan
template<typename T>
vec3<T> operator* (T s, const vec3<T> &a)
{
	return vec3<T>(a(0)*s, a(1)*s, a(2)*s);
}

//operasi antar dua vektor
//overload penjumlahan dua vektor
template<typename T>
vec3<T> operator+ (const vec3<T> &a, const vec3<T> &b)
{
    return vec3<T>(a(0)+b(0),a(1)+b(1),a(2)+b(2));
}
//overload pengurangan dua vektor
template<typename T>
vec3<T> operator- (const vec3<T> &a, const vec3<T> &b)
{
    return vec3<T>(a(0)-b(0),a(1)-b(1),a(2)-b(2));
}

//overload direct product
template<typename T>
vec3<T> operator* (const vec3<T> &a, const vec3<T> &b)
{
    return vec3<T>(a(0)*b(0),a(1)*b(1),a(2)*b(2));
}
//overload direct division
template<typename T>
vec3<T> operator/ (const vec3<T> &a, const vec3<T> &b)
{
    return vec3<T>(a(0)/b(0),a(1)/b(1),a(2)/b(2));
}

//overload buat print ke output
template<typename T>
std::ostream& operator<< (std::ostream &out, vec3<T> &v)
{
    out<<v[0]<<" "<<v[1]<<" "<<v[2];
    return out;
}

using double3 = vec3<double>;
using int3 = vec3<int>;
//masuk ke fieldnya

//di sini diletakkan container untuk vec2. (sebenernya pakai class juga bisa kan?)
template<typename T>
struct vec2
{
	vec2(const T u, const T v) : data{u,v} { }
	vec2(const T u[2] ) : data{u[0],u[1]} { }
	//perlu juga default constructor, mana tau lupa declare nilai awalnya
	vec2() : data{0,0} { }
	
	//beberapa operator sebaiknya dioverload juga
	//overload operator akses dan assignment. aksesnya ada dua
    T& operator[] (int i) {return data[i];}
    T operator() (int i) const {return data[i];}
    vec2<T>& operator=(double s) 
    {
        data[0] = s; data[1] = s; 
        return (*this);
    }
    //overload penjumlahan dan pengurangan
    vec2<T>& operator += (vec2<T> other)
    {
        for (int i=0; i<2; i++)
        {
            data[i] = data[i] + other[i];
        }
        return (*this);
    }
    vec2<T>& operator -= (vec2<T> other)
    {
        for (int i=0; i<2; i++)
        {
            data[i] = data[i] + other[i];
        }
        return (*this);
    }
protected:
	T data[2];
};

//meniru kasus vec3, operator lain dioverload juga di sini

template<typename T>
vec2<T> operator* (const vec2<T> &a, T s)
{
    return vec2<T>(a(0)*s, a(1)*s);
}
//perkalian dari kanan
template<typename T>
vec2<T> operator* (T s, const vec2<T> &a)
{
	return vec2<T>(a(0)*s, a(1)*s);
}

//operasi antar dua vektor
//overload penjumlahan dua vektor
template<typename T>
vec2<T> operator+ (const vec2<T> &a, const vec2<T> &b)
{
    return vec2<T>(a(0)+b(0),a(1)+b(1));
}
//overload pengurangan dua vektor
template<typename T>
vec2<T> operator- (const vec2<T> &a, const vec2<T> &b)
{
    return vec2<T>(a(0)-b(0),a(1)-b(1));
}

//overload direct product
template<typename T>
vec2<T> operator* (const vec2<T> &a, const vec2<T> &b)
{
    return vec2<T>(a(0)*b(0),a(1)*b(1));
}
//overload direct division
template<typename T>
vec2<T> operator/ (const vec2<T> &a, const vec2<T> &b)
{
    return vec2<T>(a(0)/b(0),a(1)/b(1));
}

//overload buat print ke output (vtk ngertinya 3 dimensi, jadi untuk entri ketiga, diset nol.
template<typename T>
std::ostream& operator<< (std::ostream &out, vec2<T> &v)
{
    out<<v[0]<<" "<<v[1] << " 0";
    return out;
}

using int2 = vec2<int>;
using double2 = vec2<double>;

//selanjutnya, dibuat class field untuk domain dua dimensi (berupa matriks rank 2)
template<typename T>
class DDField_
{
public:
    //constructor dulu, berupa inisiasi matriks rank 3 (untuk medan pada grid)
    //kita coba bedakan argumen dari constructor dengan instance pada initializer list
    DDField_(int nz, int nr) : nz{nz}, nr{nr}
    {
        data = new T*[nz];
        for (int z=0; z<nz; z++)
        {
            data[z] = new T[nr];
        }
        clear(); //ini nanti harus didefinisikan dulu
    }

    //konstruktor jenis lain
    DDField_(int2 nn) : DDField_(nn[0],nn[1]) {};

    //copy constructor, buat bikin medan yang sama 
    //ini pakai get, didefinisikan di bawah
    DDField_(const DDField_ &other) : DDField_{other.nz, other.nr}
    {
        for (int z=0; z<nz; z++)
        {
            for (int r=0; r<nr; r++)
            {
				data[z][r] = other(z,r);
            }
        }
    }

    //constructor move (kenapa argumennya pakai dua & ya)
    //initializer list buat mindahin ni, terus data dipindahin di dalam fungsinya (midahin itu sebenernya sekedar ganti nama ya?)
    DDField_(DDField_ &&other) : nz{other.nz}, nr{other.nr}
    {
        data = other.data;
        other.data=nullptr;
    }
    //operator assignment = (semacam copy constructor ya)
    DDField_& operator= (DDField_ &&f)
    {
        return (*this);
    }
    //definisikan destructornya
    ~DDField_()
    {
        if (data==nullptr) return;

        for (int z=0; z<nz; z++)
        {
            delete[] data[z];
        }
        delete[] data;
    }
    //operator akses langsung ke medannya
    T* operator[] (int i) 
    {
        return data[i];
    }
	//ini menggunakan const, menggambarkan tidak ada perubahan data
	T operator() (int z, int r) const 
	{
		return data[z][r];
	}
	
    //inisiasi nilai awal dengan skalar
    void operator= (double s)
    {
        for (int z=0; z<nz; z++)
        {
            for (int r=0; r<nr; r++)
            {
				data[z][r] = s;
            }
        }
    }
    //operator pembagian per anggota
    void operator /= (const DDField_& other)
    {
        for (int z=0; z<nz; z++)
        {
            for (int r=0; r<nr; r++)
            {
				if (other.data[z][r] != 0)
				{
					data[z][r] /= other.data[z][r];
				}
				else
				{
					data[z][r] = 0;
				}
            }
        }
    }
    //penjumlahan medan (dengan cara incremenet)
    DDField_& operator+= (const DDField_& other)
    {
        for (int z=0; z<nz; z++)
        {
            for (int r=0; r<nr; r++)
            {
				data[z][r] += other.get(z,r);
            }
        }
        return (*this);
    }
    //perkalian medan dengan skalar
    DDField_& operator*= (double s)
    {
        for (int z=0; z<nz; z++)
        {
            for (int r=0; r<nr; r++)
            {
				data[z][r] *= s;
            }
        }
        return (*this);
    }

    //operator perkalian medan dengan skalar menggunakan friend
    friend DDField_<T> operator*(double s, const DDField_<T>&f)
    {
        DDField_<T> r(f); //ini kayaknya copy constructor ya
        return r*=s;
    }

    //set seluruh entri sama dengan nol
    //kalo bentuknya kayak gini, apa artinya ini cuma bisa dipake di constructor ya?
    void clear()
    {
        (*this)=0;
    }
    //ambil data
    double get(int z, int r) const
    {
        return data[z][r];
    }
    //ini kayaknya yang bikin data bisa diakses dari luar
    //ini gak ada di buku seingetku
    template<typename S>
    friend std::ostream& operator<<(std::ostream &out, DDField_<S> &f); //kayak semacam prototype ya

    //scatter dan gather
    //scatter dilakukan untuk menentukan rapat bilangan di tiap node. Rapat bilangan ditentukan dengan metode orde pertama
    //berhubung ini untuk mesh axissymmetric, dengan volume sel yang secara umum berbeda-beda (tergantung r), maka gather dan scatternya agak sedikit berbeda
	void scatter(double2 lc, T value, double r0, double dr) //value itu bobotnya ya
    {
        //pertama tentukan indeks terkait dengan sel di mana partikel berada
        int i = (int)lc[0];
        int j = (int)lc[1];

        //pengaman
        if((i>=0) && (i<nz) && (j>=0) && (j<nr) )
        {
            //tentukan jarak ternormalisasi partikel relatif dari node i j k
            double di = lc[0] - i;
            double dj = lc[1] - j;
		
            //untuk mempermudah pembobotan volume
            double rj = r0+j*dr; //ini campur-campur antara koordinat real sama logical ya. hati-hati.
            double f1 = (rj + 0.5*dj*dr)/(rj+0.5*dr);
            double f2 = (rj + 0.5*(dj+1)*dr)/(rj+0.5*dr);
		    //update nilai rapat bilangan pada tiap node
		    data[i][j] += (T)value*(1-di)*(1-dj)*f2;
		    data[i+1][j] += (T)value*(di)*(1-dj)*f2;
		    data[i][j+1] += (T)value*(1-di)*(dj)*f1;
		    data[i+1][j+1] += (T)value*di*dj*f1;
        }
        else
        {
            outOfBound += 1;
        }
    }/**/
    void scatterCartesian(double2 lc, T value)
    {
        //pertama tentukan indeks terkait dengan sel di mana partikel berada
        int i = (int)lc[0];
        int j = (int)lc[1];
        //std::cout << "-3" << "\t";
        //std::cout << "i=" << i << "\t j=" << j << "\n";
        //std::cout << "-2" << "\t";
        //pengaman
        if((i>=0) && (i<nz) && (j>=0) && (j<nr) )
        {
            //std::cout << "-1" << "\t";
            //tentukan jarak ternormalisasi partikel relatif dari node i j k
            double di = lc[0] - i;
            double dj = lc[1] - j;
            //std::cout << "0" << "\t";
		    //update nilai rapat bilangan pada tiap node
            //std::cout << "a" << "\t";
		    data[i][j] += (T)value*(1-di)*(1-dj);
		    //std::cout << "b" << "\t";
            data[i+1][j] += (T)value*(di)*(1-dj);
		    //std::cout << "c" << "\t";
            data[i][j+1] += (T)value*(1-di)*(dj);
		    //std::cout << "d" << "\n";
            data[i+1][j+1] += (T)value*di*dj;
        }
        else
        {
            outOfBound += 1;
            std::cout << "Out of Bound \n";
        }
    }
    
    //selanjutnya gather, untuk menginterpolasikan medan sesuai dengan posisi partikel
    T gather(double2 lc, double r0, double dr)
    {
        //mirip dengan scatter, pertama tentukan indeks terkait dengan sel di mana partikel berada
        int i = (int)lc[0];
        int j = (int)lc[1];
        T val;
        if((i>=0) && (i<nz) && (j>=0) && (j<nr) )
        {
            //tentukan jarak ternormalisasi partikel relatif dari node i j k
            double di = lc[0] - i;
            double dj = lc[1] - j;
		
            //sama dengan kasus scatter
            double rj = r0+j*dr; //ini campur-campur antara koordinat real sama logical ya. hati-hati.
            double f1 = (rj + 0.5*dj*dr)/(rj+0.5*dr);
            double f2 = (rj + 0.5*(dj+1)*dr)/(rj+0.5*dr);
        
		    //interpolasi linear 2D, sebagai berikut (ada bobotnya untuk kasus axissymetric)
		    val = data[i][j]*(1-di)*(1-dj)*f2 +
			data[i+1][j]*(di)*(1-dj)*f2 +
			data[i][j+1]*(1-di)*(dj)*f1 +
			data[i+1][j+1]*di*dj*f1;
        }
        else
        {
            val = 0;
            outOfBound += 1;
        }
        return val;
    }

    //mungkin butuh satu instance lagi untuk menentukan partikel ada di dalam logam atau enggak menggunakan object_id
    T location(int i, int j)
    {
        //mirip dengan scatter, pertama tentukan indeks terkait dengan sel di mana partikel berada
        //int i = (int)lc[0];
        //int j = (int)lc[1];
        T val;
        if((i>=0) && (i<nz) && (j>=0) && (j<nr) )
        {
            val = data[i][j] + data[i+1][j] + data[i][j+1] + data[i+1][j+1];
        }
        else
        {
            val = 0;
            outOfBound +=1;
        }
        return val; //ini khusus untuk object_id
    }
	
	const int nz, nr; //jumlah grid untuk arah x, y, dan z
	
	//untuk index flattening, mulai dari i, terus j, terus k
	int index2to1(int z, int r)
	{
		return z + nz*r;
	}
	
protected:
    
    T **data;
    int outOfBound;
};

//mungkin ini setelah ada friend di class Field, bisa dicompile kenapa harus seperti ini aku belum paham
//sejelek-jeleknya pun, output nggak perlu juga dibikin code nya serapi ini.
//ini otomatis bisa langsung buat Field3 nggak ya?
template<typename T>
std::ostream& operator<<(std::ostream &out, DDField_<T> &f)
{
    for (int r=0; r<f.nr ; r++, out<<"\n")
    {
        for (int z=0; z<f.nz; z++)
        {
			out << f.data[z][r] << " ";
        }
    }
    return out;
} //harus hati-hati, di sini ni = nz, sementara nj=nr (kalau sulit dibetulkan, mungkin nanti perlu diganti dari awal)

//kayaknya nggak bisa ditulis 2DField, jadi pakai DField aja
using DField = DDField_<double>;
using DField3 = DDField_<double3>; //medannya tetep punya 3 komponen kan? atau enggak ya?
using DFieldI = DDField_<int>;
using dvector = std::vector<double>;


class Energy
{
    public:
        //constructornya dulu
        Energy(int n) : ni{n} 
        {
            data = new double[n];
            clear();
        }
        //mungkin perlu copy constructor
        Energy(const Energy &other) : Energy(other.ni)
        {
            for (int i=0; i<ni; i++)
            {
                data[i] = other(i);
            }
        }
        //ini harus ditulis, tapi kepakenya di mana ya?
        Energy(Energy &&other) : ni{other.ni}
        {
            data = other.data;
            other.data = nullptr;
        } 
        //destructornya
        ~Energy()
        {
            if (data==nullptr) 
            {
                return;
            }
            delete[] data;
        }
        void clear()
        {
            (*this)=0;
        }

        //operator assignment
        Energy& operator= (Energy&& f)
        {
            return (*this);
        }
        double operator[] (int i)
        {
            return data[i];
        }
        double operator() (int i) const
        {
            return data[i];
        }
        Energy& operator+= (Energy &f)
        {
            for (int i=0; i<ni; i++)
            {
                data[i] += f.get(i);
            }
            return (*this);
        }
        Energy& operator-= (Energy &f)
        {
            for (int i=0; i<ni; i++)
            {
                data[i] -= f.get(i);
            }
            return (*this);
        }
        //operator assignment dengan argumen satu skalar (seluruh komponen vektor memiliki nilai tersebut)
        void operator= (double s)
        {
            for (int i=0; i<ni; i++)
            {
                data[i] = s;
            }
        }
        //ambil data
        double get(int i) const
        {
            return data[i];
        }

        void setValue(int i, double a)
        {
            data[i] = a;
        }

        void multValue(int i, double a)
        {
            data[i] *= a;
        }

        //paling penting ini, butuh scatter (perlu diingat kalau E ada pada koordinat logika)
        void scatter(double E, double value)
        {
            //tentukan indeks di mana energi yang dipilih berada
            int i = (int)E;

            //tentukan jarak titik energi tersebut relatif terhadap node
            double di = E-i;
            
            //set nilai data dengan scatter orde 1
            data[i] += value*(1-di);
            data[i+1] += value*di;
        }

        //ketika plot, supaya mulus, kita gunakan interpolasi orde 1
        //sama seperti scatter, E ada pada koordinat logika
        double gather(double E)
        {
            //untuk menentukan di indeks mana interpolasi harus dilakukan
            
            int i = (int)E;
            double val;
            if(i < ni)
            {
                double di = E-i;

                //interpolasi linear
                val = data[i]*(1-di) + data[i+1]*di;
            }
            else
            {
                val = 0;
                outOfBound += 1;
            }
            return val;
        }

        const int ni;
    private:
        double *data;
        int outOfBound;
};


#endif
