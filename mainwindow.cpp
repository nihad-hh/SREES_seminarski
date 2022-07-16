#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "circuitSolver.h"

#include <QMainWindow>
#include <QCoreApplication>
#include <QFile>
#include <QTextStream>
#include <QRegularExpression>
#include <iostream>
#include <string>
#include <fstream>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    //generator parameters
    ui->Eg->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Eg->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Phi_g->setValidator( new QDoubleValidator(-360, 360, 100, this) );
    ui->Rg->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Xg->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->f->setValidator( new QDoubleValidator(0, 100, 100, this) );
    //line parameters
    ui->x0->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->r0->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->x1->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->r1->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->c0->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->c1->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->L->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->C->setValidator( new QDoubleValidator(0, 100, 100, this) );
    //consumer parameters
    ui->Ra->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Xa->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Rb->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Xb->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Rc->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Xc->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Rn->setValidator( new QDoubleValidator(0, 100, 100, this) );
    ui->Xn->setValidator( new QDoubleValidator(0, 100, 100, this) );
}

MainWindow::~MainWindow()
{

    delete ui;
}

//globalne varijable
QString naziv_potrosaca;
bool inverzija = true;

QString read_File(QString file_name)
{
    QFile file(file_name); //application current directory


    if(!file.exists())
    {
        qCritical() << "File not found";
        return "";
    }


    if(!file.open(QIODevice::OpenModeFlag::ReadOnly))
    {
        qCritical() << "Could not open file";
        qCritical() << file.errorString();
        return "";
    }

    QTextStream stream(&file); //handles encoding

    QString file_content;

    while(!stream.atEnd())
    {
        QString line = stream.readLine();
        file_content.append(line);
    }

    file.close();
    return file_content;
}

void tex_creator_ulazni_podaci(Eigen::VectorXd ulazni_podaci){

    QString file_content = read_File("src.txt");

    for(int i = 10; i<31; i++)
    {
         QString str = QString("ulaz") + QString("%1").arg(i);

         file_content.replace(str, QString::number(ulazni_podaci(i-10), 'f', 3));
    }

    file_content.replace(QString("NAZIVPOTROSACA"), naziv_potrosaca);

    QFile file("tex.txt");

    if(!file.open(QIODevice::OpenModeFlag::WriteOnly))
    {
        qCritical() << "Could not open file";
        qCritical() << file.errorString();
        return;
    }

    QTextStream stream(&file);

    stream << file_content;

    file.close();
}

void tex_creator_vector_x(Eigen::VectorXcd x){



    QString file_content = read_File("tex.txt");
    for(int i = 10; i<31; i++)
    {
         QString str1 = QString("vX%1").arg(i);
         QString str2 = QString("xPh%1").arg(i);

         file_content.replace(str1, QString::number(abs(x(i-10)), 'f', 6));
         file_content.replace(str2, QString::number(arg(x(i-10))*180/PI, 'f', 2));
    }

    QFile file("tex.txt");

    if(!file.open(QIODevice::OpenModeFlag::WriteOnly))
    {
        qCritical() << "Could not open file";
        qCritical() << file.errorString();
        return;
    }

    QTextStream stream(&file);

    stream << file_content;
    qInfo() << file_content;
    file.close();

}

void tex_creator_vector_I(Eigen::VectorXcd I, QString str){


    QString file_content = read_File("tex.txt");
    for(int i = 10; i<13; i++)
    {
         QString str1 = str + QString("%1").arg(i);
         QString str2 = str + QString("Ph%1").arg(i);

         file_content.replace(str1, QString::number(abs(I(i-10)), 'f', 3));
         file_content.replace(str2, QString::number(arg(I(i-10))*180/PI, 'f', 2));
    }

    QFile file("tex.txt");

    if(!file.open(QIODevice::OpenModeFlag::WriteOnly))
    {
        qCritical() << "Could not open file";
        qCritical() << file.errorString();
        return;
    }

    QTextStream stream(&file);

    stream << file_content;

    file.close();

}

void tex_creator_matrix_Y(Eigen::MatrixXcd Y, QString str){

    QString file_content = read_File("tex.txt");

    int max = 19;

    Y.transposeInPlace();

    for(int i = 10; i<max; i++)
    {
         QString str_oznaka = str + QString("%1").arg(i);
         QString str_zamjena = QString::number(Y.data()[i-10].real(), 'f', 6);
         if(Y.data()[i-10].imag() < 0)
             str_zamjena = str_zamjena + "-j" + QString::number(abs(Y.data()[i-10].imag()), 'f', 6);
         else
             str_zamjena = str_zamjena + "+j" + QString::number(abs(Y.data()[i-10].imag()), 'f', 6);

         file_content.replace(str_oznaka, str_zamjena);

         str_zamjena = "";
         str_oznaka = "";
    }


    QFile file("tex.txt");

    if(!file.open(QIODevice::OpenModeFlag::WriteOnly))
    {
        qCritical() << "Could not open file";
        qCritical() << file.errorString();
        return;
    }

    QTextStream stream(&file);

    stream << file_content;

    file.close();

}

void tex_creator_matrix_Yc(Eigen::MatrixXcd Y, QString str){

    QString file_content = read_File("tex.txt");

    QString replacment;

    Y.transposeInPlace();
    for(int i = 10; i<19; i++)
    {
         QString str1 = str + QString("%1").arg(i);

         if(Y.data()[i-10].imag()<0)
             replacment = "-j";
         else
            replacment = "j";

         replacment = replacment + QString::number(abs(Y.data()[i-10].imag()), 'f', 6);

         file_content.replace(str1, replacment);

    }

    QFile file("tex.tex");

    if(!file.open(QIODevice::OpenModeFlag::WriteOnly))
    {
        qCritical() << "Could not open file";
        qCritical() << file.errorString();
        return;
    }

    QTextStream stream(&file);

    stream << file_content;

    file.close();
}

void tex_creator_matrix_Zg(Eigen::MatrixXcd Y, QString str){

    QString file_content = read_File("tex.txt");

    Y.transposeInPlace();
    auto data = Y.data();

    QString str_oznaka = str + QString("%1").arg(10);
    QString str_zamjena = QString::number(data[0].real(), 'f', 3);

    if(Y.data()[1].imag() < 0)
        str_zamjena = str_zamjena + "-j" + QString::number(abs(data[0].imag()), 'f', 3);
    else
        str_zamjena = str_zamjena + "+j" + QString::number(abs(data[0].imag()), 'f', 3);

    for(int i = 10; i<13; i++)
    {
         file_content.replace(str_oznaka, str_zamjena);
    }


    QFile file("tex.txt");

    if(!file.open(QIODevice::OpenModeFlag::WriteOnly))
    {
        qCritical() << "Could not open file";
        qCritical() << file.errorString();
        return;
    }

    QTextStream stream(&file);

    stream << file_content;

    file.close();

}

void tex_creator_vector_S(Eigen::VectorXcd S, QString str){

    QString file_content = read_File("tex.txt");

    for(int i = 10; i<13; i++)
    {
         QString str_oznaka = str + QString("%1").arg(i);
         QString str_zamjena = QString::number(S(i-10).real(), 'f', 2);

         if(S(i-10).imag() < 0)
             str_zamjena = str_zamjena + "-j" + QString::number(abs(S(i-10).imag()), 'f', 3);
         else
             str_zamjena = str_zamjena + "+j" + QString::number(abs(S(i-10).imag()), 'f', 3);

         file_content.replace(str_oznaka, str_zamjena);

         str_zamjena = "";
         str_oznaka = "";
    }

    QFile file("tex.txt");

    if(!file.open(QIODevice::OpenModeFlag::WriteOnly))
    {
        qCritical() << "Could not open file";
        qCritical() << file.errorString();
        return;
    }

    QTextStream stream(&file);

    stream << file_content;

    file.close();

}

Eigen::MatrixXcd determine_consumer(QString consumer,double ra,double xa,double rb,double xb,double rc,double xc,
                                    double rn,double xn){
    Eigen::MatrixXcd Zp;

    if(consumer == "Wye grounded through impedance")
    {
        Zp = unesi_Z_zvijezde_sa_uzemljenjem(std::complex<double> (ra,xa), std::complex<double> (rb, xb),
                                             std::complex<double> (rc, xc), std::complex<double> (rn, xn));
        naziv_potrosaca = "zvijezdu uzemljenu preko impedanse";
    }
    else if(consumer == "Wye directly grounded")
    {
        Zp = unesi_Z_zvijezde_sa_uzemljenjem(std::complex<double> (ra,xa), std::complex<double> (rb, xb),
                                             std::complex<double> (rc, xc));
        naziv_potrosaca = "direktno uzemljenu zvijezdu";
    }
    else if(consumer == "Wye without the ground")
    {
        Zp = unesi_Y_zvijezde_izolovane(std::complex<double> (ra,xa), std::complex<double> (rb, xb),
                                             std::complex<double> (rc, xc));
        naziv_potrosaca = "izolovanu zvijezdu";
        inverzija = false;
    }
    else if(consumer == "Delta")
    {
        auto Za = std::complex<double> (ra,xa);
        auto Zb = std::complex<double> (rb, xb);
        auto Zc = std::complex<double> (rc, xc);
        Zp = unesi_Y_trokut(Za, Zb, Zc);

        naziv_potrosaca = "delta spoj";
        inverzija = false;
    }
    return Zp;
}

void MainWindow::on_pushButton_clicked()
{

    //generator parameters
    double eg = ui->Eg->text().toDouble();
    double phig = ui->Phi_g->text().toDouble();
    double rg = ui->Rg->text().toDouble();
    double xg = ui->Xg->text().toDouble();
    double fg = ui->f->text().toDouble();
    double w = 2*PI*fg;

    //line parameters
    double x0 = ui->x0->text().toDouble();
    double r0 = ui->r0->text().toDouble();
    double x1 = ui->x1->text().toDouble();
    double r1 = ui->r1->text().toDouble();
    double c0 = ui->c0->text().toDouble();
    double c1 = ui->c1->text().toDouble();
    double l = ui->L->text().toDouble();
    double c = ui->C->text().toDouble();

    //consumer parameters
    double ra= ui->Ra->text().toDouble();
    double xa = ui->Xa->text().toDouble();
    double rb = ui->Rb->text().toDouble();
    double xb = ui->Xb->text().toDouble();
    double rc = ui->Rc->text().toDouble();
    double xc = ui->Xc->text().toDouble();
    double rn = ui->Rn->text().toDouble();
    double xn = ui->Xn->text().toDouble();

    eg = 225;
      fg = 50;
      phig = 30;
      rg = 0.1;
      xg = 3;
      l = 80 ;
      r0 = 0.0089;
      x0 = 0.347;
      r1 = 0.017;
      x1 = 0.122;
      c0 = 1.5; //nano
      c1 = 2.3; // nano
      ra = 100;
      rb = 100;
      rc = 100;
      xa = 30;
      xb = 30;
      xc = 30;
      rn = 1;
      xn = 0;
      c = 1; //mikro
      w = 2 * PI * fg;

    //consumer type
    QString consumer = ui->consumerType->currentText();

    //Zp je Yp ako je potrosac delta ili izolovana zvijezda
    Eigen::Matrix3cd Zp = determine_consumer(consumer, ra, xa, rb, xb, rc, xc, rn, xn);
    Eigen::Matrix3cd Yp;

    if(inverzija == 1)
        Yp = Zp.inverse();
    else
        Yp = Zp;

    Yp = Yp * pow(10, 3); // 10^-3 u izvjestaju pa mnozimo sa 10^3

    //rjesavanje kola
    std::complex<double> Zgcd = std::complex<double>(rg, xg);
    Eigen::VectorXcd solved_vector = NetworkSolver(polar(eg/sqrt(3), phig*PI/180), Zgcd, Zgcd, Zgcd,
                                                       c0*pow(10, -9), c1*pow(10, -9), r0, r1, x0/w, x1/w, l, fg, c*pow(10, -6), c*pow(10, -6), c*pow(10, -6), Zp);

    //Nortronov ekvivalent
    Eigen::Matrix3cd Zn = unesi_Z_gen(Zgcd, Zgcd, Zgcd);
    Eigen::Matrix3cd Yn = Zn.inverse();

    //racunanje napona i struja
    Eigen::Vector3cd V1 = DajNaponV1(solved_vector);
    Eigen::Vector3cd V2 = DajNaponV2(solved_vector);
    Eigen::Vector3cd I12 = DajStrujuI12(solved_vector);
    Eigen::Vector3cd I21 = DajStrujuI21(solved_vector);
    //Eigen::Vector3cd I10 = DajStrujuI10(solved_vector);
    Eigen::Vector3cd I20 = DajStrujuI20(solved_vector);
    Eigen::Vector3cd Ip = DajStrujuIp(solved_vector);

    //proracun kondenzator
    Eigen::Matrix3cd Yc = CondenserBatteryModel(c*pow(10, -3), c*pow(10, -3), c*pow(10, -3), fg);

    //proracun generator
    Eigen::Matrix3cd Zg = unesi_Z_gen(Zgcd);
    Eigen::Vector3cd Eg = unesi_ems(polar(eg/sqrt(3), phig*PI/180));
    Eigen::Vector3cd In = GeneratorNortonCurrentSource(Eg, Zg);

    //proracun linija
    Eigen::MatrixXcd Y11 = LineModel_ii(c0*pow(10, -9), c1*pow(10, -9), r0, r1, x0/w, x1/w, l, fg);
    Eigen::MatrixXcd Y12 = LineModel_ij(c0, c1, r0, r1, x0/w, x1/w, l, fg);

    //racunanje snaga
    Eigen::Vector3cd Sp = Snaga(V2, Ip);
    Eigen::Vector3cd Sc = Snaga(V2, I20);
    Eigen::Vector3cd Sg = Snaga(V1, I12);
    Eigen::Vector3cd S12 = Snaga(V1, I12);
    Eigen::Vector3cd S21 = Snaga(V2, I21);
    Eigen::Vector3cd Sloss = S12 + S21;


    // ucitavanja ulaznih podataka
    Eigen::VectorXd ulazni_podaci(21);
    ulazni_podaci << eg, phig, rg, xg, fg, x0, x1, r0, r1, c0, c1, l, c, ra, xa, rb, xb, rc, xc, rn, xn;
    tex_creator_ulazni_podaci(ulazni_podaci);
    tex_creator_vector_x(solved_vector);

    //popunjavanje napona i struja
    tex_creator_vector_I(Eg, "Eg");
    tex_creator_vector_I(V1, "V1");
    tex_creator_vector_I(V2, "V2");
    Ip = Ip * pow(10, 3);
    tex_creator_vector_I(Ip, "Ip");
    In = In * pow(10, 3);
    tex_creator_vector_I(In, "In");
    I12 = I12 * pow(10, 3);
    tex_creator_vector_I(I12, "I12");
    I21 = I21 * pow(10, 3);
    tex_creator_vector_I(I21, "I21");
    I20 = I20 * pow(10, 3);
    tex_creator_vector_I(I20, "I20");

    //popunjavanje snaga
    tex_creator_vector_S(Sp, "Sp");
    tex_creator_vector_S(Sc, "Sc");
    tex_creator_vector_S(Sloss, "Sloss");
    tex_creator_vector_S(Sg, "Sg");

    //popunjavanje impedansi i admitansi
    tex_creator_matrix_Zg(Zg, "Zg");
    //tex_creator_matrix_Y(Zp, "Zp");
    tex_creator_matrix_Y(Y11, "Y11");
    tex_creator_matrix_Y(Y12, "Y12");
    tex_creator_matrix_Y(Yp, "Yp");
    tex_creator_matrix_Y(Yn, "Yn");
    tex_creator_matrix_Yc(Yc, "Yc");

    system("del /f tex.pdf");
    system("pdflatex -interaction=nonstopmode tex.tex");
    system("tex.pdf");
}

