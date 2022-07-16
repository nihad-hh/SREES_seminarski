#ifndef CIRCUITSOLVER_H
#define CIRCUITSOLVER_H

#endif // CIRCUITSOLVER_H


#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <cmath>


using namespace std::complex_literals;
const double PI = 4 * std::atan(1);
// operator a = exp(i*2*pi/3)
const std::complex<double> a = std::exp(2i * PI/3.);
extern bool inverzija;

// matrica T
Eigen::Matrix3cd T  {{1/3., 1/3., 1/3.},
                     {1/3., a/3., a*a/3.},
                     {1/3., a*a/3., a/3.},};

// Unos EMS, impedansi generatora, impedansi linija, impedansi potrosaca

// Unos jednog EMS-a, pa se ostali pomjere za 120 stepeni (Eg)
Eigen::Vector3cd unesi_ems(std::complex<double> Eg) {
    Eigen::Vector3cd Eg_V;
    // direktno
    Eg_V << Eg,
            a * a * Eg,
            a *  Eg;
    //inverzno
    /*Eg_V << Eg,
              a * a Eg,
              a * Eg */
    return Eg_V;
}

// Unos svakog EMS-a pojedinačno ( Eg1, Eg2, Eg3 )
Eigen::Vector3cd unesi_ems(std::complex<double> Eg1, std::complex<double> Eg2, std::complex<double> Eg3) {
    Eigen::Vector3cd Eg_V;
    Eg_V << Eg1,
            Eg2,
            Eg3;
    return Eg_V;
}

// Pojedinačni unos impedansi generatora ( Zg1, Zg2, Zg3 )
Eigen::Matrix3cd unesi_Z_gen(std::complex<double> Zg1, std::complex<double> Zg2, std::complex<double> Zg3) {
    Eigen::Matrix3cd m;
    m << Zg1, 0, 0,
         0, Zg2, 0,
         0, 0, Zg3;
    return m;
}

// Unos impedansi generatora ( Zg = Zg1 = Zg2 = Zg3 )
Eigen::Matrix3cd unesi_Z_gen(std::complex<double> Zg) {
    return unesi_Z_gen(Zg, Zg, Zg);
}

// Pojedinačni unos impedansi linija ( Zl1, Zl2, Zl3)
Eigen::Matrix3cd unesi_Z_lin(std::complex<double> Zl1, std::complex<double> Zl2, std::complex<double> Zl3) {
    Eigen::Matrix3cd m;
    m << Zl1, 0, 0,
         0, Zl2, 0,
         0, 0, Zl3;
    return m;
}

// Unos impedansi linija ( Zl= Zl1 = Zl2 = Zl3)
Eigen::Matrix3cd unesi_Z_lin(std::complex<double> Zl) {
    return unesi_Z_lin(Zl, Zl, Zl);
}

// Pojedinačni unos impedansi potrosaca za zvijezdu sa uzemljenjem preko impedanse i bez impedanse ( Zp1, Zp2, Zp3, Zn )
Eigen::Matrix3cd unesi_Z_zvijezde_sa_uzemljenjem(std::complex<double> Zp1, std::complex<double> Zp2, std::complex<double> Zp3, std::complex<double> Zn = 0) {
    Eigen::Matrix3cd m;
    m << Zp1 + Zn, Zn, Zn,
        Zn, Zp2 + Zn, Zn,
        Zn, Zn, Zp3 + Zn;
    return m;
}

// Unos impedansi potrosaca za zvijezdu sa uzemljenjem preko impedanse i bez impedanse ( Zp, Zn )
Eigen::Matrix3cd unesi_Z_zvijezde_sa_uzemljenjem(std::complex<double> Zp, std::complex<double> Zn  = 0) {
    return unesi_Z_zvijezde_sa_uzemljenjem(Zp, Zp, Zp, Zn);
}

// Pojedinacni unos impedansi potrosaca za zvijezdu sa izolovanim cvorom (Zp1, Zp2, Zp3)
Eigen::Matrix3cd unesi_Y_zvijezde_izolovane(std::complex<double> Zp1, std::complex<double> Zp2, std::complex<double> Zp3) {
    std::complex<double> Yp1 = 1. / Zp1;
    std::complex<double> Yp2 = 1. / Zp2;
    std::complex<double> Yp3 = 1. / Zp3;
    Eigen::Matrix3cd m;
    m << Yp1 * Yp2 + Yp1 * Yp3, -Yp1 * Yp2, -Yp1 * Yp3,
        -Yp2 * Yp1, Yp2* Yp1 + Yp2 * Yp3, -Yp2 * Yp3,
        -Yp3 * Yp1, -Yp3 * Yp2, Yp3* Yp1 + Yp3 * Yp2;
    m = 1. / (Yp1 + Yp2 + Yp3) * m;

    return m;
}

// Unos impedansi potrosaca za zvijezdu sa izolovanim cvorom (Zp)
Eigen::Matrix3cd unesi_Y_zvijezde_izolovane(std::complex<double> Zp) {
    return unesi_Y_zvijezde_izolovane(Zp, Zp, Zp);
}

// Pojedinacni unos impedansi potrosaca za trokut
Eigen::Matrix3cd unesi_Y_trokut(std::complex<double> Zp1, std::complex<double> Zp2, std::complex<double> Zp3) {
    Eigen::Matrix3cd m = 3 * unesi_Y_zvijezde_izolovane(Zp1, Zp2, Zp3);
    return m;
}

// Unos impedansi potrosaca za trokut
Eigen::Matrix3cd unesi_Y_trokut(std::complex<double> Zp) {
    Eigen::Matrix3cd m = unesi_Y_trokut(Zp, Zp, Zp);
    return m;
}

//Funkcija koja daje vektor struja ekvivalentnog Nortonovog strujnog izvora
Eigen::Vector3cd GeneratorNortonCurrentSource(Eigen::Vector3cd Eg, Eigen::Matrix3cd Zg) {
    return Zg.inverse() * Eg;
}

//Funkcija koja vraca matricu admitansi(6x6) pi-cetveropola
Eigen::MatrixXcd LineModel(double c0, double c1, double r0, double r1, double L0, double L1, double l, double f) {
    std::complex<double> Ysh0(0, 2*PI*f*c0*l), Ysh1(0, 2*PI*f*c1*l);
    std::complex<double> Zs0(r0*l, 2*PI*f*L0*l), Zs1(r1*l, 2*PI*f*L1*l);

    Eigen::Matrix3cd Yl11 = T.inverse()*(Eigen::Matrix3cd {{1./Zs0 + Ysh0/2., 0, 0},
                                                           {0, 1./Zs1 + Ysh1/2., 0},
                                                           {0, 0, 1./Zs1 + Ysh1/2.},})*T;

    Eigen::Matrix3cd Yl12 = T.inverse()*(Eigen::Matrix3cd {{-1./Zs0, 0, 0},
                                                            {0, -1./Zs1, 0},
                                                            {0, 0, -1./Zs1}})*T;
    Eigen::MatrixXcd Yl(6,6);

    Yl << Yl11, Yl12,
          Yl12, Yl11;

    return Yl;
}

// Funkcija koja vraca matricu admitansi cetveropola Yii
Eigen::MatrixXcd LineModel_ii(double c0, double c1, double r0, double r1, double L0, double L1, double l, double f) {
    std::complex<double> Ysh0(0, 2*PI*f*c0*l), Ysh1(0, 2*PI*f*c1*l);
    std::complex<double> Zs0(r0*l, 2*PI*f*L0*l), Zs1(r1*l, 2*PI*f*L1*l);

    Eigen::Matrix3cd Yl11 = T.inverse()*(Eigen::Matrix3cd {{1./Zs0 + Ysh0/2., 0, 0},
                                                           {0, 1./Zs1 + Ysh1/2., 0},
                                                           {0, 0, 1./Zs1 + Ysh1/2.},})*T;
    return Yl11;
}

// Funkcija koja vraca matricu admitansi cetveropola Yij
Eigen::MatrixXcd LineModel_ij(double c0, double c1, double r0, double r1, double L0, double L1, double l, double f) {
    std::complex<double> Ysh0(0, 2*PI*f*c0*l), Ysh1(0, 2*PI*f*c1*l);
    std::complex<double> Zs0(r0*l, 2*PI*f*L0*l), Zs1(r1*l, 2*PI*f*L1*l);

    Eigen::Matrix3cd Yl12 = T.inverse()*(Eigen::Matrix3cd {{-1./Zs0, 0, 0},
                                                            {0, -1./Zs1, 0},
                                                            {0, 0, -1./Zs1}})*T;
    return Yl12;
}

//Funkcija za modeliranje kondenzatorske baterije kada su kapacitivnosti iste
Eigen::Matrix3cd CondenserBatteryModel(double c, double f){
    std::complex<double> Z(0, -1./2*PI*f*c);
    return unesi_Y_trokut(Z);
}

//Funkcija za modeliranje kondenzatorske baterije u slucaju razlicitih kapacitivnosti
Eigen::Matrix3cd CondenserBatteryModel(double c1, double c2, double c3, double f) {
    std::complex<double> Y1(0, (2.*PI*f*c1)), Y2(0, (2.*PI*f*c2)), Y3(0, (2.*PI*f*c3));

    Eigen::Matrix3cd Yc;

    Yc << Y1 + Y2, -Y1, -Y3,
          -Y1, Y1 + Y3, -Y2,
          -Y3, -Y2, Y2 + Y3;

    return Yc;
}

Eigen::VectorXcd NetworkSolver(std::complex<double> E, std::complex<double> Zg1, std::complex<double> Zg2, std::complex<double> Zg3, // parametri generatora
                               double c0, double c1, double r0, double r1, double L0, double L1, double l, double f,// parametri linija + frekvencija
                               double C1, double C2, double C3, // parametri kondenzatorske baterije
                               Eigen::Matrix3cd Zp) // parametri potrosaca u vidu matrice potrosaca (prakticnije je obaviti unos potrosaca prije poziva solvera)
{
    // modeliranje generatora
    Eigen::Vector3cd Eg = unesi_ems(E);
    Eigen::Matrix3cd Zg = unesi_Z_gen(Zg1, Zg2, Zg3);
    Eigen::Vector3cd In = GeneratorNortonCurrentSource(Eg, Zg); //nortonov izvor
    Eigen::Matrix3cd Yg = Zg.inverse();

    // modeliranje linija - pi ekvivalent
    Eigen::MatrixXcd Yl = LineModel(c0, c1, r0, r1, L0, L1, l, f);

    // modeliranje kondenzatorske baterije
    Eigen::Matrix3cd Yc = CondenserBatteryModel(C1, C2, C3, f);

    // modeliranje potrosaca
    Eigen::Matrix3cd Yp;

    if(inverzija == true)
        Yp = Zp.inverse();
    else
        Yp = Zp;

    // RJESAVANJE MREZE
    Eigen::Matrix3cd zeroM = Eigen::Matrix3cd::Zero(3,3); // matrica nula
    Eigen::Matrix3cd eyeM = Eigen::Matrix3cd::Identity(3,3); // jedinicna matrica
    Eigen::Vector3cd zeroV(0, 0, 0); // nula-vektor

    // matrica sistema - A
    Eigen::MatrixXcd A(21,21);

    A << Yl.topLeftCorner(3,3), Yl.topRightCorner(3,3), -1.*eyeM, zeroM, zeroM, zeroM, zeroM,
         Yl.bottomLeftCorner(3,3), Yl.bottomRightCorner(3,3), zeroM, -1.*eyeM, zeroM, zeroM, zeroM,
         zeroM, zeroM, eyeM, zeroM, eyeM, zeroM, zeroM,
         zeroM, zeroM, zeroM, eyeM, zeroM, eyeM, eyeM,
         Yg, zeroM, zeroM, zeroM, -1.*eyeM, zeroM, zeroM,
         zeroM, Yc, zeroM, zeroM, zeroM, -1.*eyeM, zeroM,
         zeroM, Yp, zeroM, zeroM, zeroM, zeroM, -1.*eyeM;

    // vektor kolona "poznatih"
    Eigen::VectorXcd b(21);
    b << zeroV, zeroV, In, zeroV, zeroV, zeroV, zeroV;

    // rjesenje: [V1, V2, I12, I21, I10, I20, Ip].transpose()
   return A.inverse()*b;

}

Eigen::Vector3cd DajNaponV1(Eigen::VectorXcd x){
    return Eigen::Vector3cd {x(0),x(1),x(2)};
}

Eigen::Vector3cd DajNaponV2(Eigen::VectorXcd x){
    return Eigen::Vector3cd {x(3),x(4),x(5)};
}

Eigen::Vector3cd DajStrujuI12(Eigen::VectorXcd x){
    return Eigen::Vector3cd {x(6),x(7),x(8)};
}

Eigen::Vector3cd DajStrujuI21(Eigen::VectorXcd x){
    return Eigen::Vector3cd {x(9),x(10),x(11)};
}

Eigen::Vector3cd DajStrujuI10(Eigen::VectorXcd x){
    return Eigen::Vector3cd {x(12),x(13),x(14)};
}

Eigen::Vector3cd DajStrujuI20(Eigen::VectorXcd x){
    return Eigen::Vector3cd {x(15),x(16),x(17)};
}

Eigen::Vector3cd DajStrujuIp(Eigen::VectorXcd x){
    return Eigen::Vector3cd {x(18),x(19),x(20)};
}

Eigen::Vector3cd DajStruje(Eigen::Vector3cd Eg, Eigen::Matrix3cd Zg, Eigen::Matrix3cd Zl, Eigen::Matrix3cd Zp) {
    return (Zg + Zl + Zp).inverse() * Eg;
}

// Ispis u Polarnim Koordinatama
void PolarPrint(Eigen::VectorXcd x){
    for (int i=0;i<x.rows();i++){
        std::cout<<"("<<abs(x(i))<<" <"<<arg(x(i)) *(180./PI)<<")"<<std::endl;
    }
}

//Dijagonalizacija vektora 3x1
Eigen::Matrix3cd Diagonalize(Eigen::Vector3cd v){
    return Eigen::Matrix3cd{{v(0), 0, 0},
                            {0, v(1), 0},
                            {0, 0,  v(2)}, };
}

// Funkcija za racunanje snage
Eigen::Vector3cd Snaga(Eigen::Vector3cd V, Eigen::Vector3cd I){
    return Diagonalize(V) * I.conjugate();
}

//Funkcija za racunanje ukupne snage iz Vektora 3x1
std::complex<double> UkupnaSnaga(Eigen::Vector3cd V){
    return V.sum();
}

int main(){
    // Primjer zadatka s c2 za provjeru tačnosti solvera
    double l = 80, f = 50;
    std::complex<double>  E = 225./std::sqrt(3)*std::exp(1i*PI/6.), Zg(0.1,3);
    double r0 = 0.089 , r1 = 0.017, L0 = 0.347/(2.*PI*f), L1 = 0.122/(2.*PI*f), c0 = 1.5e-9, c1 = 2.3e-9, c = 1e-6;
    std::complex<double> zp(200, 80), zn(1,0);
    Eigen::Matrix3cd Zp = unesi_Z_zvijezde_sa_uzemljenjem(zp, zn);

    std::cout << NetworkSolver(E, Zg, Zg, Zg, c0, c1, r0, r1, L0, L1, l, f, c, c, c, Zp);

    return 0;
}
