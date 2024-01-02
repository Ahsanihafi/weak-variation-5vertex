#include"Species.h"
#include<math.h>
#include<iostream>
#include"Field.h"
#include<fstream>
#include<assert.h>
double Species::getRealCount()
{
    double summ = 0;
    for (Particle &part:particles)
    {
        summ += part.mpw;
    }
    return summ;
}

void Species::computeNumberDensity()
{
    den.clear(); //jangan lupa
    for (Particle &part:particles)
    {
        double2 lc = world.XtoL({part.pos[0],part.pos[1]});
        den.scatter(lc, part.mpw, world.getX0()[1], world.getDh()[1]);
    }
    //setelah diberi penyematan, harus dibagi dengan volume node
    den /= world.node_vol;
}

//ini tanpa ada "conditioning" di awal
void Species::addParticle(double3 pos, double3 vel, double mpw)
{
    if (!world.inBounds({pos[0],pos[1]}))
    {
        return; //kok bisa pakai return ya, padahal void
    }
    
    particles.emplace_back(pos,vel,mpw);
}

void Species::addParticle2(double3 pos, double3 vel, double mpw)
{
    if (!world.inBounds({pos[0],pos[1]}))
    {
        return; //kok bisa pakai return ya, padahal void
    }
    //by default, posisi azimuthal ada di theta = 0.
    //untuk addParticle, nggak harus subcycling
    double dt = -0.5*world.getDt();
    
    double t = charge*magnetz*dt/(2*mass);
    double s = 2*t/(1+t*t);

    double2 lc2 = world.XtoL({pos[0],pos[1]});
    double3 ef_part = (0.5*charge*dt/mass)*world.ef.gather(lc2,world.getX0()[1], world.getDh()[1]);

    double3 velmin = vel + ef_part;
    double3 v1 = {velmin[0], velmin[1] + velmin[2]*t, velmin[2]-velmin[1]*t};
    //akhirnya, diperoleh velplus
    double3 velplus = {velmin[0], velmin[1] + v1[2]*s, velmin[2] - v1[1]*s};
    //kecepatan final dan posisi final diberikan sebagai berikut, langsung update nilainya
    vel = velplus + ef_part; 
    
    particles.emplace_back(pos,vel,mpw);
}
//ini untuk cyclotronic versi Delzanno, kayaknya bermasalah untuk r kecil [121122]
void Species::addParticle3(double3 pos, double3 vel, double mpw)
{
    advanceCyclotronicDrift(pos,vel,world.getDt()*0.5);
    if (!world.inBounds({pos[0],pos[1]}))
    {
        return; //kok bisa pakai return ya, padahal void
    }
    particles.emplace_back(pos,vel,mpw);
}

void Species::addParticle4(double3 pos, double3 vel, double mpw)
{
    double omega0 = charge*magnetz/mass;
    advanceCyclotronicDriftCartesian(omega0,pos,vel,world.getDt()*0.5);
    negativeRotation(pos,vel);
    if (!world.inBounds({pos[0],pos[1]}))
    {
        return; //kok bisa pakai return ya, padahal void
    }
    particles.emplace_back(pos,vel,mpw);
}

//advance tanpa magnet dan nonrelativistik, menggunakan leapfrog
void Species::advance()
{
    double dt = world.getDt();
    for(Particle &part:particles)
    {
        double2 lc = world.XtoL({part.pos[0],part.pos[1]});
        double3 ef_part = world.ef.gather(lc, world.getX0()[1], world.getDh()[1]);
        //ini versi ringkasnya, kalau dibandingkan addParticle di atas
        part.vel += ef_part*(dt*charge/mass);
        part.pos += part.vel*dt; //ini bener gini ya? vel nya udah keupdate duluan dong
        //kalau keluar bound, langsung diubah bobotnya jadi nol
        if(!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue; //habiskan semua partikel dulu
        }
    }
    //cek mana aja partikel yang bobotnya nol, kemudian dihapus
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); //yang kosong (selisih yang terisi dengan awal) dihapus
}

//versi dengan magnet, dan dengan koordinat z-rho-phi
//dengan simetri domain, maka selalu diasumsikan bahwa phi=0

//bagian penggerak partikel (kayaknya interpolasinya langsung masuk sini aja ya)
void Species::borisPush(double3& pos, double3& vel, double dt, double t, double s)
{
    //set koordinat logika
    double2 lc2 = world.XtoL({pos[0],pos[1]});    
    //medan listrik tetap 3 dimensi (Er = Ex asalkan sudut phi = 0)
    double3 ef_part = (0.5*charge*dt/mass)*world.ef.gather(lc2,world.getX0()[1], world.getDh()[1]); //sekalian definisikan di sini
    //bagian algoritma boris, proses terjadi di koordinat kartesian
    //pertama pisahkan bagian listrik dengan magnet
    double3 velmin = vel + ef_part;
    //berikutnya, tentukan v1, t, dan s (besaran pembantu terkait rotasi) di sini
    double3 v1 = {velmin[0], velmin[1] + velmin[2]*t, velmin[2]-velmin[1]*t};
    //akhirnya, diperoleh velplus
    double3 velplus = {velmin[0], velmin[1] + v1[2]*s, velmin[2] - v1[1]*s};
    //kecepatan final dan posisi final diberikan sebagai berikut, langsung update nilainya
    vel = velplus + ef_part; 
    pos = pos + dt*vel; //secara umum pos[2] != 0 , artinya harus dirotasikan balik
}

//untuk partikel bermuatan negatif, arah rotasinya terbalik
void Species::negativeRotation(double3& pos, double3& vel)
{
    //rotasikan balik ke bidang simulasi. Kecepatan juga ikut dirotasikan (ingat, koordinatnya kartesian, bukan silindris)
    double r = sqrt(pos[1]*pos[1]+pos[2]*pos[2]);
    double costheta, sintheta;
    if (r!=0)
    {
        costheta = pos[1]/r;
        sintheta = pos[2]/r;
    }
    else
    {
        costheta = 1;
        sintheta = 0;
    }
    vel[1] = costheta*vel[1] + sintheta*vel[2];
    vel[2] = - sintheta*vel[1] + costheta*vel[2];

    pos = {pos[0],r,0}; //ini semacam putar balik
}

//rotasi untuk partikel bermuatan positif
void Species::positiveRotation(double3& pos, double3& vel)
{
    //rotasikan balik ke bidang simulasi. Kecepatan juga ikut dirotasikan (ingat, koordinatnya kartesian, bukan silindris)
    double r = sqrt(pos[1]*pos[1]+pos[2]*pos[2]);
    double costheta, sintheta;
    if (r!=0)
    {
        costheta = pos[1]/r;
        sintheta = pos[2]/r;
    }
    else
    {
        costheta = 1;
        sintheta = 0;
    }
    vel[1] = costheta*vel[1] - sintheta*vel[2];
    vel[2] = sintheta*vel[1] + costheta*vel[2];

    pos = {pos[0],r,0}; //ini semacam putar balik
}
//jadi int &dinding cuma penting untuk zona 5-6-7-8
//edit 18Januari2022, let's try the easiest way first. Jadi partikel langsung delete aja, spawn partikel baru
//partikel baru digerakkan di langkah berikutnya. Mestinya ini bisa dibetulkan lagi.
//ditemukan kesalahan terkait crossing z=0 [210922]
bool Species::boundaryCheck(double2 posawal, double2 &posakhir, double &dt, double3 &vel, double &mpw, int& dinding)
{
    //parameter dinding ditulis di sini

    if(world.zoneOneBound(posakhir)) //untuk zona 1 dan 6, perlu dipastikan lagi partikel menumbuk dinding atau bablas
    {
        world.reflectionSheathOne(posawal, posakhir, dt, vel, dinding);
        double celah = world.getXm()[0] - world.getGap()*world.getDh()[0];
        if (posakhir[0] >= celah) // nggak menumbuk dinding
        {
            mpw = 0; //langsung delete aja
            world.deleteAtas();
        }
        else //menumbuk permukaan 1
        {
            double vv = (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
            double Ek = 0.5*mass*vv/abs(charge); //spek dinding 
            //ini dt awal udah diubah nilainya kan. jadi partikel yang lama didelete
            double inc_theta = calcThetaBawah(vel);
            double delta_e = calcDelta(world.getDeltaMax(), Ek, world.getE0Wall(), world.getEMaxWall(), inc_theta);
            secondEmissionBawah(delta_e, {posakhir[0],posakhir[1],0}, vel, vv, world.getTempWall());
            mpw = 0;
            world.deleteAtas();
        }
        return true;
    }
    else if(world.zoneTwoBound(posakhir)) 
    {
        world.reflectionSheathTwo(posawal, posakhir, dt, vel, dinding);
        //pasti mantul dari kanan
        double vv = (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
        double Ek = 0.5*mass*vv/abs(charge); //spek dinding 
        //ini dt awal udah diubah nilainya kan. jadi partikel yang lama didelete
        double inc_theta = calcThetaSamping(vel);
        double delta_e = calcDelta(world.getDeltaMax(), Ek, world.getE0Wall(), world.getEMaxWall(), inc_theta);
        secondEmissionCathode(delta_e, {posakhir[0],posakhir[1],0}, vel, vv, world.getTempWall());
        mpw = 0;
        world.deleteKanan();
        return true;
    }
    else if(world.zoneThreeBound(posakhir)) //ini dipantulkan balik (berarti nggak diapa-apakan ya)
    {
        world.reflectionSheathThree(posawal, posakhir, dt, vel, dinding);
        return true;
    }
    else if(world.zoneFourBound(posakhir)) //ini juga dipantulkan balik
    {
        world.reflectionSheathFour(posawal, posakhir, dt, vel, dinding);
        return true;
    }
    else if(world.zoneFiveBound(posakhir)) //tergantung nilai dinding 
    {
        world.reflectionSheathOne(posawal, posakhir, dt, vel, dinding);
        //kalo numbuk atas, pake skema pantulan, kalo samping, nggak usah diapa-apakan
        if(dinding == 2)
        {
            double vv = (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
            double Ek = 0.5*mass*vv/abs(charge); //spek dinding 
            //ini dt awal udah diubah nilainya kan. jadi partikel yang lama didelete
            double inc_theta = calcThetaBawah(vel);
            double delta_e = calcDelta(world.getDeltaMax(), Ek, world.getE0Wall(), world.getEMaxWall(), inc_theta);
            secondEmissionBawah(delta_e, {posakhir[0],posakhir[1],0}, vel, vv, world.getTempWall());
            mpw = 0;
            world.deleteAtas();
        }
        return true;
    }
    else if(world.zoneSixBound(posakhir)) //tergantung nilai dinding, perlu tambahan seperti zona 1.
    {
        world.reflectionSheathOne(posawal, posakhir, dt, vel, dinding);
        if(dinding == 1) //nabrak samping
        {
            double vv = (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
            double Ek = 0.5*mass*vv/abs(charge); //spek dinding 
            //ini dt awal udah diubah nilainya kan. jadi partikel yang lama didelete
            double inc_theta = calcThetaSamping(vel);
            double delta_e = calcDelta(world.getDeltaMax(), Ek, world.getE0Wall(), world.getEMaxWall(), inc_theta);
            secondEmissionCathode(delta_e, {posakhir[0],posakhir[1],0}, vel, vv, world.getTempWall());
            mpw = 0;
            world.deleteKanan();
        }
        else if(dinding == 2) //nabrak atas
        {
            double celah = world.getXm()[0] - world.getGap()*world.getDh()[0];
            if (posakhir[0] >= celah) // nggak menumbuk dinding
            {
                mpw = 0; //langsung delete aja
                world.deleteAtas();
            }
            else //menumbuk permukaan 1
            {
                double vv = (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
                double Ek = 0.5*mass*vv/abs(charge); //spek dinding 
                //ini dt awal udah diubah nilainya kan. jadi partikel yang lama didelete
                double inc_theta = calcThetaBawah(vel);
                double delta_e = calcDelta(world.getDeltaMax(), Ek, world.getE0Wall(), world.getEMaxWall(), inc_theta);
                secondEmissionBawah(delta_e, {posakhir[0],posakhir[1],0}, vel, vv, world.getTempWall());
                mpw = 0;
                world.deleteAtas();
            }
        }
        return true;
    }
    else if(world.zoneSevenBound(posakhir)) //tergantung nilai dinding
    {
        world.reflectionSheathThree(posawal, posakhir, dt, vel, dinding);
        if(dinding == 1)
        {
            double vv = (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
            double Ek = 0.5*mass*vv/abs(charge); //spek dinding 
            //ini dt awal udah diubah nilainya kan. jadi partikel yang lama didelete
            double inc_theta = calcThetaSamping(vel);
            double delta_e = calcDelta(world.getDeltaMax(), Ek, world.getE0Wall(), world.getEMaxWall(), inc_theta);
            secondEmissionCathode(delta_e, {posakhir[0],posakhir[1],0}, vel, vv, world.getTempWall());
            mpw = 0;
            world.deleteKanan();
            //ini kan sebatas buat elektron ya
        }
        //kalo dinding == 2, biarkan apa adanya aja
        return true;
    }
    else if(world.zoneEightBound(posakhir)) //ini kedua permukaannya merupakan permukaan simetri, jadi nggak diapa-apakan
    {
        world.reflectionSheathThree(posawal, posakhir, dt, vel, dinding);
        return true;
    }
    else //ini dilanjutkan dengan algoritma null collision
    {
        return false;
    }
}

void Species::nullCollision(double3 &vel, double& mpw, double3 pos, int &elastic_count)
{
    double peluang = rnd(); //ini dibalik dari sebelumnya ya, kalo kurang dari peluang, maka algoritma tumbukan dijalankan. 
    double vv = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
    double Ek = 0.5*mass*vv/abs(charge);
    if (peluang < cross_section.getProbability())
    {
        if(Ek >= cross_section.getEmax())
        {
            //std::cout << "ada partikel dengan energi lebih besar dibandingkan dengan Emax \n";
            //nggak perlu continue kan, soalnya udah ada else.
            world.addTooHigh();
        }
        else
        {
            double le = cross_section.EtoL(Ek);
            double peluang2 = rnd();
            //katanya, meskipun ionisasi, sudut hamburannya tetap cosx, artinya hanya tergantung pada energi awal saja
            //diingat semua prosesnya: 
            /*
            proses[1] = elastis H2
            proses[2] = H+
            proses[3] = H2+
            proses[4] = Dissociative attachment
            proses[5] = elastis H
            
            */
           if (peluang2 <= cross_section.proses[1].gather(le)) //tumbukan elastis antara e dan H2
           {
               momTransfer(vel, vv, Ek);
               elastic_count += 1;
           }
           else if((peluang2 > cross_section.proses[1].gather(le)) && (peluang2 <= cross_section.proses[2].gather(le))) //ionisasi dengan disosiasi (e+H2 -> e + e + H + H+)
           {
               //energi ionisasinya
               double Eion = 13.6;
               //pastikan energinya memang lebih
               if(Ek <= Eion)
               {
                   std::cout << "\t Energi partikel lumayan rendah, tetapi terjadi ionisasi \n";
               }
               else //memang terjadi ionisasi
               {
                   ionization(vel, pos, mpw, Eion, vv, Ek);
                   world.addIonHPlus(pos);
               }
           }
           else if((peluang2 > cross_section.proses[2].gather(le)) && (peluang2 <= cross_section.proses[3].gather(le))) //ionisasi tanpa disosiasi (e + H2 -> e + e + H2+)
           {
               //energi ionisasi
               double Eion = 13.6;
               //pastikan energinya memang lebih
               if(Ek <= Eion)
               {
                   std::cout << "\t Energi partikel lumayan rendah, tetapi terjadi ionisasi \n";
               }
               else //memang terjadi ionisasi
               {
                   ionization(vel, pos, mpw, Eion, vv, Ek);
                   world.addIonH2Plus(pos);
               }
           }
           else if((peluang2 > cross_section.proses[3].gather(le)) && (peluang2 <= cross_section.proses[4].gather(le))) //seluruh dissociative electron attachment coba digabung menjadi satu (e + H2 -> H + H-)
           {
               world.addIonHMin(pos);
               mpw = 0; //delete elektron pada proses DA.
           }
           else if ((peluang2 > cross_section.proses[4].gather(le)) && (peluang2 <= cross_section.proses[5].gather(le))) //ionisasi e + H -> e + e + H+
           {
               //energi ionisasi
               double Eion = 13.6;
               //pastikan energinya memang lebih
               if(Ek <= Eion)
               {
                   std::cout << "\t Energi partikel lumayan rendah, tetapi terjadi ionisasi \n";
               }
               else //memang terjadi ionisasi
               {
                   ionization(vel, pos, mpw, Eion, vv, Ek);
                   world.addIonHPlus(pos);
               }
           }
           else if ((peluang2 > cross_section.proses[5].gather(le)) && (peluang2 <= cross_section.proses[6].gather(le))) //tumbukan elastis antara e dan H
           {
               momTransfer(vel,vv,Ek);
               elastic_count += 1;
           }
           //else null, nggak perlu ditulis sih. kalau nggak ada kondisi yang dipenuhi, nggak perlu ngelakukan apa-apa.
        }
    }
}

//jadi ini mending algoritma penentu jenis interaksi dengan dinding, dan algoritma null collision disendirikan di luar advance


void Species::advanceBorisElektron(double dt)
{
    int nelastic = 0;
    double t = charge*magnetz*dt/(2*mass);
    double s = 2*t/(1+t*t);
    for (Particle &part:particles)
    {
        int dinding_awal = 0;
        //jangan lupa, pos[0] = z, pos[1] = x, pos[2] = y
        double2 posawal = {part.pos[0],part.pos[1]};
        //untuk menentukan dt awal (mempermudah dalam menentukan sisa waktu setelah tumbukan)
        double dt_awal = dt;
        //pastikan partikel benar-benar di dalam domain
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        //dorong dengan boris push
        borisPush(part.pos, part.vel, dt_awal, t, s);
        //rotasikan balik
        negativeRotation(part.pos, part.vel);
        double2 posakhir = {part.pos[0],part.pos[1]};
        //ini versi yang lebih baru terkait skema elektron sekunder, untuk menentukan lokasi tumbukan dengan dinding
        //ini belum bener ding. belum mempertimbangkan waktu tumbukan, harusnya dipertimbangkan juga. (jadi partikel bergerak lumayan setelah tumbukan jauh gitu.)
        if(!boundaryCheck(posawal, posakhir, dt_awal, part.vel, part.mpw, dinding_awal))
        {
            //null collision
            nullCollision(part.vel, part.mpw, part.pos, nelastic);
            continue; //jadi di bawah ini, dilakukan lagi boris push menggunakan sisa waktu yang belum habis
            //kalau udah jelas di dalam domain, nggak perlu dikerjain, jadi tinggal kasih continue.
        }
        part.pos[0] = posakhir[0];
        part.pos[1] = posakhir[1];
        //sementara kasih itu dulu, nanti harus dilacak salahnya di mana (segmentation fault)
        /*
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        */
        //jadi misalkan udah numbuk dinding, set mpw=0. kalau melewati batas simetri, lanjutkan gerakan dengan dt_awal yang sudah diupdate
        if (part.mpw == 0)
        {
            continue;
        }
        else
        {
            borisPush(part.pos, part.vel, dt_awal, t, s);
            negativeRotation(part.pos,part.vel); //diputar lagi. mestinya sekali aja cukup sih ini.
        }
        
    }
    //hapus partikel yang berbobot nol
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); //yang kosong (selisih yang terisi dengan awal) dihapus
    //std::cout << "\t Jumlah tumbukan elastik = " << nelastic << "\t";
}

//ini advance boris elektron dengan skema produksi elektron sekunder yang sederhana [051022]
//bedanya dengan yang ion, nanti untuk ion bisa produksi elektron sekunder pada katoda, dan sementara ion nggak pakai null collision
void Species::advanceBorisElektron2(double dt)
{
    double t = charge*magnetz*dt/(2*mass);
    double s = 2*t/(1+t*t);
    for (Particle &part:particles)
    {
        int nelastic = 0;
        //kita pakai versi sederhana dulu [051022]
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        borisPush(part.pos, part.vel, dt, t, s);
        negativeRotation(part.pos,part.vel);

        //skema setelah gerakannya kayak gini [051022]
        if((part.pos[0] >= world.getXm()[0])) //partikel menumbuk katoda, rapat muatan permukaan berubah
        {
            world.addChargeIncrement(charge, part.mpw);
            world.deleteKanan();
            part.mpw = 0;
            continue; 
        }
        else if(part.pos[1] >= world.getXm()[1]) //partikel menumbuk anoda, langsung delete
        {
            part.mpw = 0;
            world.deleteAtas();
            continue;
        }
        else if ( part.pos[1] <= world.getX0()[1] )//partikel melewati garis r = 0
        {
            part.vel[1] = -part.vel[1]; //balikkan kecepatan arah r
            part.pos[1] =  2*world.getX0()[1] - part.pos[1];
            continue;
        }
        else if(part.pos[0] <= world.getX0()[0])
        {
            part.vel[0] = -part.vel[0];
            part.pos[0] = 2*world.getX0()[0] - part.pos[0];
            continue;
        }
        else //ini kalo masih di dalam, kenakan monte carlo collision
        {
            //cek dulu kalau partikel ada di wilayah dinding logam atau tidak [251022]
            double2 loc_part = world.XtoL({part.pos[0],part.pos[1]});
            int loc_val = world.object_id.location(loc_part[0],loc_part[1]);
            if (loc_val == 4)
            {
                part.mpw = 0;
                continue;
            }
            nullCollision(part.vel, part.mpw, part.pos, nelastic); 
        }
    }
    //hapus partikel yang berbobot nol
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); //yang kosong (selisih yang terisi dengan awal) dihapus
}

//ini untuk algoritma cyclotronic termodifikasi [101122]

void Species::advanceCyclotronicDriftCartesian(double omega, double3 &pos, double3 &vel, double dt)
{
    double3 x1 = pos;
    double3 v1 = vel;

    pos[0] = x1[0] + v1[0]*dt; //ini z
    pos[1] = x1[1] + (v1[2]-v1[2]*cos(omega*dt)+v1[1]*sin(omega*dt))/omega; //ini x
    pos[2] = x1[2] + (-v1[1] + v1[1]*cos(omega*dt) + v1[2]*sin(omega*dt) )/omega; //ini y

    vel[1] = v1[1]*cos(omega*dt) + v1[2]*sin(omega*dt);
    vel[2] = -v1[1]*sin(omega*dt) + v1[2]*cos(omega*dt);

}

void Species::advanceCyclotronicMod(double dt)
{
    int nelastic = 0;
    double omega0 = charge*magnetz/mass; //bisa negatif ya
    for (Particle &part:particles)
    {
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        //dengan rotasi koordinat untuk meniru gaya inersial
        advanceCyclotronicKick(part.pos, part.vel, dt);
        advanceCyclotronicDriftCartesian(omega0, part.pos, part.vel, dt);
        negativeRotation(part.pos, part.vel);
        //advanceCyclotronicDriftCartesian(omega0, part.pos, part.vel, 0.5*dt);
        //negativeRotation(part.pos, part.vel);
        if((part.pos[0] >= world.getXm()[0])) //partikel menumbuk katoda, rapat muatan permukaan berubah
        {
            world.addChargeIncrement(charge, part.mpw);
            world.deleteKanan();
            part.mpw = 0;
            continue; 
        }
        else if(part.pos[1] >= world.getXm()[1]) //partikel menumbuk anoda, langsung delete
        {
            part.mpw = 0;
            world.deleteAtas();
            continue;
        }
        else if ( part.pos[1] <= world.getX0()[1] )//partikel melewati garis r = 0
        {
            part.vel[1] = -part.vel[1]; //balikkan kecepatan arah r
            part.pos[1] =  2*world.getX0()[1] - part.pos[1];
            continue;
        }
        else if(part.pos[0] <= world.getX0()[0])
        {
            part.vel[0] = -part.vel[0];
            part.pos[0] = 2*world.getX0()[0] - part.pos[0];
            continue;
        }
        else //kalo masih dalam bound, hitung nullcollisionnya
        {
            double2 loc_part = world.XtoL({part.pos[0],part.pos[1]});
            int loc_val = world.object_id.location(loc_part[0],loc_part[1]);
            if (loc_val == 4)
            {
                part.mpw = 0;
                continue;
            }
            nullCollision(part.vel, part.mpw, part.pos, nelastic);
        }
    }
    //hapus partikel yang berbobot nol
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); //yang kosong (selisih yang terisi dengan awal) dihapus
    
}

void Species::advanceCyclotronicModIon(double dt)
{
    int nelastic = 0;
    double omega0 = charge*magnetz/mass; //bisa negatif ya
    for (Particle &part:particles)
    {
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        //dengan rotasi koordinat untuk meniru gaya inersial
        advanceCyclotronicKick(part.pos, part.vel, dt);
        advanceCyclotronicDriftCartesian(omega0, part.pos, part.vel, dt);
        negativeRotation(part.pos, part.vel);
        //advanceCyclotronicDriftCartesian(omega0, part.pos, part.vel, 0.5*dt);
        //negativeRotation(part.pos, part.vel);
        if((part.pos[0] >= world.getXm()[0])) //partikel menumbuk katoda, rapat muatan permukaan berubah
        {
            world.addChargeIncrement(charge, part.mpw);
            world.deleteKanan();
            part.mpw = 0;
            double emissionChance = rnd();
            if (emissionChance <= world.getGamma())
            {
                //fungsi kreasi elektron. tambah di suatu counter, terus dispawn pakai source. Jangan lupa dikurangi muatan (negatif) pada rapat muatan permukaan [061022]
                world.addChargeIncrement(charge, part.mpw); //satu elektron lepas sama dengan satu ion positif masuk (kecuali kalo yang numbuk ion negatif, tapi kayaknya nggak mungkin)
                //skemanya mungkin sementara elektronnya kecepatannya nol ya [061022]
                world.addElecCathode({world.getXm()[0] - 0.5*world.getDh()[0],part.pos[1],0}, {0,0,0});
            }
            continue; 
        }
        else if(part.pos[1] >= world.getXm()[1]) //partikel menumbuk anoda, langsung delete
        {
            part.mpw = 0;
            world.deleteAtas();
            continue;
        }
        else if ( part.pos[1] <= world.getX0()[1] )//partikel melewati garis r = 0
        {
            part.vel[1] = -part.vel[1]; //balikkan kecepatan arah r
            part.pos[1] =  2*world.getX0()[1] - part.pos[1];
            continue;
        }
        else if(part.pos[0] <= world.getX0()[0])
        {
            part.vel[0] = -part.vel[0];
            part.pos[0] = 2*world.getX0()[0] - part.pos[0];
            continue;
        }
        else //kalo masih dalam bound, hitung nullcollisionnya
        {
            double2 loc_part = world.XtoL({part.pos[0],part.pos[1]});
            int loc_val = world.object_id.location(loc_part[0],loc_part[1]);
            if (loc_val == 4)
            {
                part.mpw = 0;
                continue;
            }
        }
    }
    //hapus partikel yang berbobot nol
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); //yang kosong (selisih yang terisi dengan awal) dihapus
    
}

//integrator cyclotronic kita masukkan di sini [151022]
void Species::advanceCyclotronicDrift(double3 &pos, double3 &vel, double dt)
{
    //tetap menggunakan koordinat z-r-theta
    double3 x1 = pos;
    double3 v1 = vel;
    double omega = charge*magnetz/mass;
    double epsilon = 0.5*(v1[1]*v1[1]+v1[2]*v1[2]); //kecepatan non-axial ya
    double Ptheta = x1[1]*v1[2]+0.5*omega*x1[1]*x1[1];
    //update posisi
      
    pos[0] = x1[0] + v1[0]*dt;
    pos[1] = sqrt(x1[1]*x1[1]*cos(omega*dt) + (4*epsilon+2*omega*Ptheta)*(1-cos(omega*dt))/(omega*omega) + 2*x1[1]*v1[1]*sin(omega*dt)/omega);
    pos[2] = x1[2] ; //males ngitung bagian theta, belum kepakai toh axially symmetric
     
    //update kecepatan
    vel[0] = v1[0]; //pada bagian drift, nggak ada perubahan kecepatan arah z.
    vel[1] = (x1[1]*v1[1]*cos(omega*dt) - 0.5*omega*x1[1]*x1[1]*sin(omega*dt) + (2*epsilon + omega*Ptheta)*sin(omega*dt)/omega )/pos[1] ;
    vel[2] = Ptheta/pos[1] - 0.5*omega*pos[1]; 
}
void Species::advanceCyclotronicKick(double3 &pos, double3 &vel, double dt)
{
    if (world.inBounds({pos[0],pos[1]}))
    {
        double2 lc = world.XtoL({pos[0],pos[1]});
        double3 ef_part = world.ef.gather(lc, world.getX0()[1], world.getDh()[1]);
        vel += (charge*dt/mass)*ef_part;
    }
}

void Species::advanceElectronCyclotronic(double dt)
{
    int nelastic = 0;
    for (Particle &part:particles)
    {
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        advanceCyclotronicKick(part.pos,part.vel,dt);
        advanceCyclotronicDrift(part.pos,part.vel,dt);
        if(part.pos[0] >= world.getXm()[0])
        {
            world.addChargeIncrement(charge, part.mpw);
            part.mpw = 0;
            world.deleteKanan();
            continue; //habiskan semua partikel dulu (kalau di bawah ngapain pake continue ya)
        }
        else if(part.pos[1] >= world.getXm()[1])
        {
            part.mpw = 0;
            world.deleteAtas();
            continue; //habiskan semua partikel dulu (kalau di bawah ngapain pake continue ya)
        }
        else if ( part.pos[1] <= world.getX0()[1] )//partikel melewati garis r = 0
        {
            part.vel[1] = -part.vel[1]; //balikkan kecepatan arah r
            part.pos[1] =  2*world.getX0()[1] - part.pos[1];
            continue;
        }
        else if(part.pos[0] <= world.getX0()[0])
        {
            part.vel[0] = -part.vel[0];
            part.pos[0] = 2*world.getX0()[0] - part.pos[0];
        }
        else //kalo masih dalam bound, hitung nullcollisionnya
        {
            double2 loc_part = world.XtoL({part.pos[0],part.pos[1]});
            int loc_val = world.object_id.location(loc_part[0],loc_part[1]);
            if (loc_val == 4)
            {
                part.mpw = 0;
                continue;
            }
            nullCollision(part.vel, part.mpw, part.pos, nelastic);
        }
    }
    //hapus partikel yang berbobot nol
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); 
}

void Species::advanceIonCyclotronic(double dt)
{
    int nelastic = 0;
    for (Particle &part:particles)
    {
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        advanceCyclotronicDrift(part.pos,part.vel,0.5*dt);
        advanceCyclotronicKick(part.pos,part.vel,dt);
        advanceCyclotronicDrift(part.pos,part.vel,0.5*dt);
        if(part.pos[0] >= world.getXm()[0]) //partikel menumbuk dinding katoda
        {
            world.addChargeIncrement(charge, part.mpw);
            part.mpw = 0;
            //untuk produksi elektron sekunder, kita tiru skemanya.
            //kalau pakai data abadi, gamma = 0.1, artinya satu elektron sekunder untuk sepuluh ion yang menumbuk
            //mendingan dibikin acak aja ya, kalo rnd() < 0.1, maka spawn elektron sekunder. kalo enggak ya enggak. [061022]
            double emissionChance = rnd();
            if (emissionChance <= world.getGamma())
            {
                //fungsi kreasi elektron. tambah di suatu counter, terus dispawn pakai source. Jangan lupa dikurangi muatan (negatif) pada rapat muatan permukaan [061022]
                world.addChargeIncrement(charge, part.mpw); //satu elektron lepas sama dengan satu ion positif masuk (kecuali kalo yang numbuk ion negatif, tapi kayaknya nggak mungkin)
                //skemanya mungkin sementara elektronnya kecepatannya nol ya [061022]
                world.addElecCathode({world.getXm()[0] - 0.5*world.getDh()[0],part.pos[1],0}, {0,0,0});
            }
            continue; 
        }
        else if(part.pos[1] >= world.getXm()[1])
        {
            part.mpw = 0;
            continue;
        }
        else if ( part.pos[1] <= world.getX0()[1] )//partikel melewati garis r = 0
        {
            part.vel[1] = -part.vel[1]; //balikkan kecepatan arah r
            part.pos[1] =  2*world.getX0()[1] - part.pos[1];
            continue;
        }
        else if(part.pos[0] <= world.getX0()[0])
        {
            part.vel[0] = -part.vel[0];
            part.pos[0] = 2*world.getX0()[0] - part.pos[0];
        }
        else
        {
            double2 loc_part = world.XtoL({part.pos[0],part.pos[1]});
            int loc_val = world.object_id.location(loc_part[0],loc_part[1]);
            if (loc_val == 4)
            {
                part.mpw = 0;
                continue;
            }
        }
    }
    //hapus partikel yang berbobot nol
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); 
}

//kita bikin satu lagi untuk ion, yang kalo numbuk menghasilkan elektron sekunder, mungkin ada yang lebih efisien [061022]
void Species::advanceBorisIon(double dt)
{
    double t = charge*magnetz*dt/(2*mass);
    double s = 2*t/(1+t*t);
    for (Particle &part:particles)
    {
        int nelastic = 0;
        //kita pakai versi sederhana dulu [051022]
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        borisPush(part.pos, part.vel, dt, t, s);
        negativeRotation(part.pos,part.vel);

        //skema setelah gerakannya kayak gini [051022]
        if((part.pos[0] >= world.getXm()[0])) //partikel menumbuk katoda, rapat muatan permukaan berubah
        {
            world.addChargeIncrement(charge, part.mpw);
            part.mpw = 0;
            //untuk produksi elektron sekunder, kita tiru skemanya.
            //kalau pakai data abadi, gamma = 0.1, artinya satu elektron sekunder untuk sepuluh ion yang menumbuk
            //mendingan dibikin acak aja ya, kalo rnd() < 0.1, maka spawn elektron sekunder. kalo enggak ya enggak. [061022]
            double emissionChance = rnd();
            if (emissionChance <= world.getGamma())
            {
                //fungsi kreasi elektron. tambah di suatu counter, terus dispawn pakai source. Jangan lupa dikurangi muatan (negatif) pada rapat muatan permukaan [061022]
                world.addChargeIncrement(charge, part.mpw); //satu elektron lepas sama dengan satu ion positif masuk (kecuali kalo yang numbuk ion negatif, tapi kayaknya nggak mungkin)
                //skemanya mungkin sementara elektronnya kecepatannya nol ya [061022]
                world.addElecCathode({world.getXm()[0] - 0.5*world.getDh()[0],part.pos[1],0}, {0,0,0});
            }
            continue; 
        }
        else if(part.pos[1] >= world.getXm()[1]) //partikel menumbuk anoda, langsung delete
        {
            part.mpw = 0;
            continue;
        }
        else if ( part.pos[1] <= world.getX0()[1] )//partikel melewati garis r = 0
        {
            part.vel[1] = -part.vel[1]; //balikkan kecepatan arah r
            part.pos[1] =  2*world.getX0()[1] - part.pos[1];
            continue;
        }
        else if(part.pos[0] <= world.getX0()[0])
        {
            part.vel[0] = -part.vel[0];
            part.pos[0] = 2*world.getX0()[0] - part.pos[0];
            continue;
        }
        else
        {
            double2 loc_part = world.XtoL({part.pos[0],part.pos[1]});
            int loc_val = world.object_id.location(loc_part[0],loc_part[1]);
            if (loc_val == 4)
            {
                part.mpw = 0;
                continue;
            }
        }
        //sementara nggak pakai mcc, belum ada datanya [061022]
    }
    //hapus partikel yang berbobot nol
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); //yang kosong (selisih yang terisi dengan awal) dihapus
}


//untuk menentukan koefisien emisi sekunder (sudut datang ditentukan di luar kali ya)
//ini cari data parameternya dulu, sementara bisa ngikutin grafik yang udah ada (tapi jadinya ya agak ngawur sih)
double Species::calcDelta(double delta_max, double E, double Eth, double Emax, double theta)
{
    //ini asumsi seberapa halusnya permukaan
    double ksw = 1; //typical value
    double w = (E-Eth)/(Emax*(1 + ksw*theta*theta/(2*Const::PI) ) - Eth );
    double k; //konstanta pangkat
    if (w < 1)
    {
        k = 0.56;
    }
    else if(w>=1 && w<= 3.6)
    {
        k = 0.25;
    }
    else
    {
        k = 0;
    }
    double Ww; //container untuk menentukan ketergantungan koefisien delta terhadap energi
    if(w <= 3.6)
    {
        Ww = (w*exp(1-w))*k;
    }
    else
    {
        Ww = 1.125/pow(w,0.35);
    }
    //ini konstanta terkait smoothness juga
    double ksd = 1;

    return delta_max*(1+ksd*theta*theta/(2*Const::PI))*Ww;
}

//ini tergantung kecepatan arah azimuth juga
double Species::calcThetaSamping(double3 vel)
{
    double vv = sqrt(vel[1]*vel[1]+vel[2]*vel[2]);
    double tantheta = vv/vel[0];
    return atan(tantheta); //otomatis bisa buat kecepatan negatif juga ya
}

double Species::calcThetaBawah(double3 vel)
{
    double vv = sqrt(vel[1]*vel[1]+vel[2]*vel[2]);
    double tantheta = vel[0]/vv;
    return atan(tantheta);
}
//kecepatan acaknya langsung dari sini aja ya. nggak di source
void Species::secondEmissionSamping(double delta_e, double3 pos, double3 vel, double vv, double temp) 
{
    double vth = sqrt(2*temp*abs(Const::QE)/mass); //temperatur dalam electron volt

    int int_delta = (int)delta_e; double del_delta = delta_e - int_delta;
    double r = rnd();
    int jum_par;
    double3 tempVel; 
    if (r < del_delta)
    {
        jum_par = int_delta + 1;
    }
    else if (r >= del_delta)
    {
        jum_par = int_delta;
    }
    for (int i=0; i<jum_par; i++)
    {
        double cc = rnd(); //jadi mungkin emang didefinisikan di sini jenis emisinya
        if(cc<0.9) //elektron sekunder murni
        {
            //kecepatan untuk spawn di permukaan mestinya ya beda dengan untuk kasus spawn di volume ya. Pake hemisphere gitu.
            if(vel[0]>0) //permukaan di sebelah kiri nih
            {
                double sin_theta = rnd(); //nggak ada sintheta negatif
                double cos_theta = sqrt(1-cos_theta*cos_theta);
                double phii = 2*Const::PI*rnd();
                tempVel = {-vth*cos_theta, vth*sin_theta*cos(phii), vth*sin_theta*sin(phii)};
            }
            else
            {
                double sin_theta = rnd(); //nggak ada sintheta negatif
                double cos_theta = sqrt(1-cos_theta*cos_theta);
                double phii = 2*Const::PI*rnd();
                tempVel = {vth*cos_theta, vth*sin_theta*cos(phii), vth*sin_theta*sin(phii)};
            }
        }
        else if((cc >= 0.9) && (cc<0.97)) //backscattering
        {
            if (vel[0] > 0) //datang dari kiri
            {
                double multiplier = rnd();
                tempVel = {-multiplier*vel[0], multiplier*vel[1], multiplier*vel[2]}; //kecepatan arah z harus dibalik
            }
            else //datang dari kanan
            {
                double multiplier = rnd();
                tempVel = {multiplier*vel[0], multiplier*vel[1], multiplier*vel[2]};
            }
        }
        else //refleksi
        {
            if (vel[0] > 0) //datang dari kiri
            {
                tempVel = {-vel[0], vel[1], vel[2]}; //kecepatan arah z harus dibalik

            }
            else //datang dari kanan
            {
                tempVel = {vel[0], vel[1], vel[2]};
            }
        }
        world.addElecWall(pos, tempVel); 
    }
}

void Species::secondEmissionCathode(double delta_e, double3 pos, double3 vel, double vv, double temp) 
{
    double vth = sqrt(2*temp*abs(Const::QE)/mass); //temperatur dalam electron volt

    int int_delta = (int)delta_e; double del_delta = delta_e - int_delta;
    double r = rnd();
    int jum_par;
    double3 tempVel; 
    if (r < del_delta)
    {
        jum_par = int_delta + 1;
    }
    else if (r >= del_delta)
    {
        jum_par = int_delta;
    }
    for (int i=0; i<jum_par; i++)
    {
        double cc = rnd(); //jadi mungkin emang didefinisikan di sini jenis emisinya
        if(cc<0.9) //elektron sekunder murni
        {
            //kecepatan untuk spawn di permukaan mestinya ya beda dengan untuk kasus spawn di volume ya. Pake hemisphere gitu.
            if(vel[0]>0) //permukaan di sebelah kiri nih
            {
                double sin_theta = rnd(); //nggak ada sintheta negatif
                double cos_theta = sqrt(1-sin_theta*sin_theta);
                double phii = 2*Const::PI*rnd();
                tempVel = {-vth*cos_theta, vth*sin_theta*cos(phii), vth*sin_theta*sin(phii)};
            }
            else
            {
                double sin_theta = rnd(); //nggak ada sintheta negatif
                double cos_theta = sqrt(1-sin_theta*sin_theta);
                double phii = 2*Const::PI*rnd();
                tempVel = {vth*cos_theta, vth*sin_theta*cos(phii), vth*sin_theta*sin(phii)};
            }
        }
        else if((cc >= 0.9) && (cc<0.97)) //backscattering
        {
            if (vel[0] > 0) //datang dari kiri
            {
                double multiplier = rnd();
                tempVel = {-multiplier*vel[0], multiplier*vel[1], multiplier*vel[2]}; //kecepatan arah z harus dibalik
            }
            else //datang dari kanan
            {
                double multiplier = rnd();
                tempVel = {multiplier*vel[0], multiplier*vel[1], multiplier*vel[2]};
            }
        }
        else //refleksi
        {
            if (vel[0] > 0) //datang dari kiri
            {
                tempVel = {-vel[0], vel[1], vel[2]}; //kecepatan arah z harus dibalik

            }
            else //datang dari kanan
            {
                tempVel = {vel[0], vel[1], vel[2]};
            }
        }
        world.addElecCathode({world.getXm()[0] - 0.01*world.getDh()[0],pos[1],0}, tempVel); 
    }
}

//untuk kasus ini masih ragu-ragu, karena permukaannya kan mengarah ke r, sementara kecepatannya diarahkan ke vx sama vy.
void Species::secondEmissionBawah(double delta_e, double3 pos, double3 vel, double vv, double temp) 
{
    double vth = sqrt(2*temp*abs(Const::QE)/mass); //temperatur dalam electron volt
    int int_delta = (int)delta_e; double del_delta = delta_e - int_delta;
    double r = rnd();
    int jum_par;
    double3 tempVel;
    
    //untuk pengujian, kita gunakan cout
    //std::cout << "\t nilai r = " << r << "\t delta_e = " << delta_e << "\t"; 
    if (r < del_delta)
    {
        jum_par = int_delta + 1;
    }
    else if (r >= del_delta)
    {
        jum_par = int_delta;
    }
    for (int i=0; i<jum_par; i++)
    {
        //tentukan jenis emisinya terlebih dahulu, oya, ini logikanya pasti kecepatan partikel arah r positif
        double cc = rnd();
        if (cc < 0.9) //emisi sekunder beneran, selalu muncul ke bawah kan
        {
            double sin_theta = rnd(); //nggak ada sintheta negatif
            double cos_theta = sqrt(1-sin_theta*sin_theta);
            double phii = 2*Const::PI*rnd();
            tempVel = {vth*sin_theta*cos(phii),-vth*cos_theta,vth*sin_theta*sin(phii)};
        }
        else if((cc >= 0.9) && (cc<0.97)) //backscattering
        {
            double multiplier = rnd();
            tempVel = {multiplier*vel[0], -multiplier*vel[1], -multiplier*vel[2]};
        }
        else
        {
            tempVel = {vel[0], -vel[1], -vel[2]};
        }
        //std::cout << "jum_par = " << jum_par << "\n" ;
        world.addElecWall({pos[0],world.getXm()[1] - 0.01*world.getDh()[1],0}, tempVel);
    }
}
void Species::excitation(double3 &vel, double Eexc, double vv, double Ek) //excitation sekedar energy sink
{
    double normav = sqrt(vv);
    //normalisasi, untuk menentukan vektor arah
    double3 vi = {vel[0]/normav, vel[1]/normav, vel[2]/normav};
    double sinth = sqrt(1.0-vi[1]*vi[1]);
    double phi1 = 2*Const::PI*rnd();
    double R = rnd();
    double cosx = (2 + Ek - 2*pow(1+Ek,R))/Ek; //oya, ini untuk elektron yang numbuk gas ya. Kalo ion yang numbuk ya beda
    double sinx = sqrt(1.0-cosx*cosx);
    double Escat = Ek - Eexc;
    double vscat = sqrt(2*Escat/mass);
    double vvx = vi[1]*cosx - (vi[2]*vi[2]+vi[0]*vi[0])*sinx*cos(phi1)/sinth;
    double vvy = vi[2]*cosx + vi[0]*sinx*sin(phi1)/sinth + vi[2]*vi[1]*sinx*cos(phi1)/sinth;
    double vvz = vi[0]*cosx - vi[2]*sinx*sin(phi1)/sinth + vi[0]*vi[1]*sinx*cos(phi1)/sinth;

    vel[0] = vscat*vvz;
    vel[1] = vscat*vvx;
    vel[2] = vscat*vvy; //ini energinya udah berkurang
}

void Species::momTransfer(double3 &vel, double vv, double Ek)
{
    double phi = 2*Const::PI*rnd();
    double R = rnd();
    double cosx = (2 + Ek - 2*pow(1+Ek,R))/Ek; //sementara, untuk tumbukan elastik, tidak ada transfer energi menuju gas netral
    double sinx = sqrt(1.0-cosx*cosx);
    double normav = sqrt(vv);
    double3 vi = {vel[0]/normav, vel[1]/normav, vel[2]/normav};
    double sinth = sqrt(1.0-vi[1]*vi[1]);
    double vvx = vi[1]*cosx - (vi[2]*vi[2]+vi[0]*vi[0])*sinx*cos(phi)/sinth;
    double vvy = vi[2]*cosx + vi[0]*sinx*sin(phi)/sinth + vi[2]*vi[1]*sinx*cos(phi)/sinth;
    double vvz = vi[0]*cosx - vi[2]*sinx*sin(phi)/sinth + vi[0]*vi[1]*sinx*cos(phi)/sinth;
    vel[0] = normav*vvz;
    vel[1] = normav*vvx;
    vel[2] = normav*vvy;
}
void Species::ionization(double3 &vel, double3 pos, double mpw, double Eion, double vv, double Ek)
{
    double normav = sqrt(vv);
    //normalisasi, untuk menentukan vektor arah
    double3 vi = {vel[0]/normav, vel[1]/normav, vel[2]/normav};
    double sinth = sqrt(1.0-vi[1]*vi[1]);
    double phi1 = 2*Const::PI*rnd();
    double phi2 = 2*Const::PI*rnd();
    double R = rnd();
    double Eej;
    double cosx = (2 + Ek - 2*pow(1+Ek,R))/Ek; //sementara, untuk tumbukan elastik, tidak ada transfer energi menuju gas netral
    double sinx = sqrt(1.0-cosx*cosx);
    //kita pisahkan kasus di bawah 20 eV dengan di atasnya, karena untuk kasus di bawah 20 eV, ada pendekatan yang dapat dilakukan.
    if (Ek <= 20) //sekali lagi, ini dalam eV
    {
        Eej = R*0.5*(Ek-Eion);
    }
    else //ini nggak pakai pendekatan, lebih ribet juga. Tapi ini berlaku kurang lebih sampai 70 eV
    {
        //ini angka 10 di bawah cuma aproksimasi, yang valid sampai 70eV saja, buat atom jenis lain gimana ya?
        Eej = 10*tan(R*atan( (Ek-Eion)/20.0 ));
    }
    //kayaknya nggak harus menaati kelestarian momentum ya, jadi arah elektron yang terhambur dengan yang terejeksi nggak harus berlawanan
    double Escat = Ek - Eion - Eej;
    double vscat = sqrt(abs(2*Escat*charge/mass)); //amplitudo scattered
    double vej = sqrt(abs(2*Eej*charge/mass)); //amplitudo ejected

    //semua yang di bawah ini sudah ternormalisir
    double vvx = vi[1]*cosx - (vi[2]*vi[2]+vi[0]*vi[0])*sinx*cos(phi1)/sinth;
    double vvy = vi[2]*cosx + vi[0]*sinx*sin(phi1)/sinth + vi[2]*vi[1]*sinx*cos(phi1)/sinth;
    double vvz = vi[0]*cosx - vi[2]*sinx*sin(phi1)/sinth + vi[0]*vi[1]*sinx*cos(phi1)/sinth;
    double vex = vi[1]*cosx - (vi[2]*vi[2]+vi[0]*vi[0])*sinx*cos(phi2)/sinth;
    double vey = vi[2]*cosx + vi[0]*sinx*sin(phi2)/sinth + vi[2]*vi[1]*sinx*cos(phi2)/sinth;
    double vez = vi[0]*cosx - vi[2]*sinx*sin(phi2)/sinth + vi[0]*vi[1]*sinx*cos(phi2)/sinth;

    vel[0] = vscat*vez;
    vel[1] = vscat*vex;
    vel[2] = vscat*vey;
    world.addElecIon(pos,{vez*vej,vex*vej,vey*vej});
    world.addCountIonizeElectron();
}

void Species::advanceBoris2(double dt)
{
    int nelastic = 0;
    double t = charge*magnetz*dt/(2*mass);
    double s = 2*t/(1+t*t);
    for (Particle &part:particles)
    {
        int dinding_awal = 0;
        //jangan lupa, pos[0] = z, pos[1] = x, pos[2] = y
        double2 posawal = {part.pos[0],part.pos[1]};
        //untuk menentukan dt awal (mempermudah dalam menentukan sisa waktu setelah tumbukan)
        double dt_awal = dt;
        //pastikan partikel benar-benar di dalam domain
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        //dorong dengan boris push
        borisPush(part.pos, part.vel, dt_awal, t, s);
        //rotasikan balik
        negativeRotation(part.pos, part.vel);
        double2 posakhir = {part.pos[0],part.pos[1]};
        //ini versi yang lebih baru terkait skema elektron sekunder, untuk menentukan lokasi tumbukan dengan dinding
        //ini belum bener ding. belum mempertimbangkan waktu tumbukan, harusnya dipertimbangkan juga. (jadi partikel bergerak lumayan setelah tumbukan jauh gitu.)
        if(!boundaryCheck(posawal, posakhir, dt_awal, part.vel, part.mpw, dinding_awal))
        {
            //null collision
            //nullCollision(part.vel, part.mpw, part.pos, nelastic);
            continue; //jadi di bawah ini, dilakukan lagi boris push menggunakan sisa waktu yang belum habis
            //kalau udah jelas di dalam domain, nggak perlu dikerjain, jadi tinggal kasih continue.
        }
        part.pos[0] = posakhir[0];
        part.pos[1] = posakhir[1];
        
        //jadi misalkan udah numbuk dinding, set mpw=0. kalau melewati batas simetri, lanjutkan gerakan dengan dt_awal yang sudah diupdate
        if (part.mpw == 0)
        {
            continue;
        }
        else
        {
            borisPush(part.pos, part.vel, dt_awal, t, s);
            negativeRotation(part.pos,part.vel); //diputar lagi. mestinya sekali aja cukup sih ini.
        }

    }
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); //yang kosong (selisih yang terisi dengan awal) dihapus
    //std::cout << "\t Jumlah tumbukan elastik = " << nelastic << "\t";
}


void Species::advanceBoris3(double dt)
{
    int nelastic = 0;
    double t = charge*magnetz*dt/(2*mass);
    double s = 2*t/(1+t*t);
    for (Particle &part:particles)
    {
        int dinding_awal = 0;
        //jangan lupa, pos[0] = z, pos[1] = x, pos[2] = y
        double2 posawal = {part.pos[0],part.pos[1]};
        //untuk menentukan dt awal (mempermudah dalam menentukan sisa waktu setelah tumbukan)
        double dt_awal = dt;
        //pastikan partikel benar-benar di dalam domain
        if (!world.inBounds({part.pos[0],part.pos[1]}))
        {
            part.mpw = 0;
            continue;
        }
        //dorong dengan boris push
        borisPush(part.pos, part.vel, dt_awal, t, s);
        //rotasikan balik
        negativeRotation(part.pos, part.vel);
        double2 posakhir = {part.pos[0],part.pos[1]};
        //ini versi yang lebih baru terkait skema elektron sekunder, untuk menentukan lokasi tumbukan dengan dinding
        //ini belum bener ding. belum mempertimbangkan waktu tumbukan, harusnya dipertimbangkan juga. (jadi partikel bergerak lumayan setelah tumbukan jauh gitu.)
        if(!boundaryCheck(posawal, posakhir, dt_awal, part.vel, part.mpw, dinding_awal))
        {
            //null collision
            //nullCollision(part.vel, part.mpw, part.pos, nelastic);
            continue; //jadi di bawah ini, dilakukan lagi boris push menggunakan sisa waktu yang belum habis
            //kalau udah jelas di dalam domain, nggak perlu dikerjain, jadi tinggal kasih continue.
        }
        part.pos[0] = posakhir[0];
        part.pos[1] = posakhir[1];
        
        //jadi misalkan udah numbuk dinding, set mpw=0. kalau melewati batas simetri, lanjutkan gerakan dengan dt_awal yang sudah diupdate
        if (part.mpw == 0)
        {
            continue;
        }
        else
        {
            borisPush(part.pos, part.vel, dt_awal, t, s);
            negativeRotation(part.pos,part.vel); //diputar lagi. mestinya sekali aja cukup sih ini.
        }

    }
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw > 0)
        {
            continue;
        }
        particles[p] = particles[np-1]; //yang paling belakang dikopi
        np--; //kecilkan panjang list, artinya yang paling belakang dihapus
        p--; //supaya urutan yang sama dicek sekali lagi
    }
    //selanjutnya tinggal mengecilkan ukuran std::vector<Particle>
    particles.erase(particles.begin()+np, particles.end()); //yang kosong (selisih yang terisi dengan awal) dihapus
    //std::cout << "\t Jumlah tumbukan elastik = " << nelastic << "\t";
}

void Species::advanceM(double dt)
{
    //ambil time step yang digunakan pada world
    //double dt = world.getDt(); //jangan lupa didefinisikan
    
    //sementara medan magnet ditulis di sini, mestinya sih nanti di potensialsolver. tapi berhubung homogen kan, asal masukin angka aja.
    double Bz = 0.0;

    //loop terhadap seluruh partikel
    for(Particle &part:particles)
    {
        //tentukan posisi partikel dengan grid ternormalisasi
        double2 lc = world.XtoL({part.pos[0],part.pos[1]});
        //tentukan medan listrik pada posisi partikel
        double3 ef_part = world.ef.gather(lc,world.getX0()[1],world.getDh()[1]);
        //update kecepatan partikel
        part.vel[0] = part.vel[0] + ef_part[0]*(dt*charge/mass) ;
        part.vel[1] = part.vel[1] + ef_part[1]*(dt*charge/mass) + part.vel[2]*Bz*charge*dt/mass;
        part.vel[2] = part.vel[2] + ef_part[2]*(dt*charge/mass) - part.vel[1]*Bz*charge*dt/mass;
        //update posisi partikel
        part.pos += part.vel*dt;
        
        double r = sqrt(part.pos[1]*part.pos[1]+part.pos[2]*part.pos[2]);
        part.pos = {part.pos[0],r,0};
        //kalo nggak terpenuhi ya biarkan aja, lanjut ke bawah. jadi nggak perlu else ya?

        //untuk kasus simulasi ini, kalau partikel keluar domain, bobot partikel diset nol
        if(!world.inBounds({part.pos[0],part.pos[1]}))
        {
            //addTrail(part.pos-part.vel*dt,part.mpw);
            part.mpw = 0; //set bobot sama dengan nol
            continue; //harus pake ini ya?
        }
    }
    //yang bobotnya nol, kemudian dihapus. Partikel tetap disusun menggunakan vektor
    //sementara nggak pakai linked list, kebanyakan yang diubah
    size_t np = particles.size();
    for (size_t p=0; p<np; p++)
    {
        if(particles[p].mpw>0)
        {
            continue;
        }
        //di bawah ini cuma berlaku untuk yang mpw = 0
        particles[p] = particles[np-1]; //kayaknya ini nggak perlu overload, yang disamadengankan member vektor, bukan vec3<double> nya
        np--; //batasnya dikurangi (soalnya entri matriksnya berkurang)
        p--; //supaya titik yang sama dicek lagi.
    }
    particles.erase(particles.begin()+np, particles.end() ); //ini buat ngecilin matriksnya
}