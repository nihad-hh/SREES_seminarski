#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "circuitSolver.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


//void MainWindow::on_btnOK_clicked()
//{
//    QString firstName = ui->leFirstName->text();
//    QString lastName = ui->leLastName->text();
//    ui->lbWelcome->setText("Welcome" + firstName + " " + lastName);
//}


//void MainWindow::on_btnClear_clicked()
//{
//    ui->leFirstName->setText("");
//    ui->leLastName->setText("");
//    ui->lbWelcome->setText("");
//}


void MainWindow::on_pushButton_clicked()
{
    //generator parameters
    double eg = ui->Eg->text().toDouble();
    double rg = ui->Rg->text().toDouble();
    double xg = ui->Xg->text().toDouble();
    double f = ui->f->text().toDouble();
    double w = 2*PI*f;

    std::complex<double> Eg(eg, 0), Zg(rg, xg);

    //line parameters
    double x0 = ui->x0->text().toDouble();
    double r0 = ui->r0->text().toDouble();
    double x1 = ui->x1->text().toDouble();
    double r1 = ui->r1->text().toDouble();
    double c0 = ui->c0->text().toDouble();
    double c1 = ui->c1->text().toDouble();
    double l = ui->L->text().toDouble();
    double c = ui->C->text().toDouble();

    std::complex<double> Z0(r0*l, x0*l), Z1(r1*l, x1*l), Ysh0(0, 2*PI*f*c0*l), Ysh1(0, 2*PI*f*c1*l), Yc(2*PI*f*c);

    //consumer parameters
    double ra= ui->Ra->text().toDouble();
    double xa = ui->Xa->text().toDouble();
    double rb = ui->Rb->text().toDouble();
    double xb = ui->Xb->text().toDouble();
    double rc = ui->Rc->text().toDouble();
    double xc = ui->Xc->text().toDouble();
    double rn = ui->Rn->text().toDouble();
    double xn = ui->Xn->text().toDouble();

    std::complex<double> Za(ra, xa), Zb(rb, xb), Zc(rc, xc), Zn(rn, xn);

    //consumer type

    QString consumer = ui->consumerType->currentText();

    // napisati logiku za odabir potrosaca
    //if();

    Eigen::Matrix3cd Zp;
    qInfo() << consumer;
    qInfo() << Zg.imag();
    qInfo() << Zg.real();


    //if(consumer.toString() == "Wye")

    //za probu - potrosac u zvijezdu
    Zp << Za, 0,  0,
          0,  Zb, 0,
          0,  0,  Zc;

//    Eigen::VectorXcd NetworkSolver(std::complex<double> E, std::complex<double> Zg1, std::complex<double> Zg2, std::complex<double> Zg3, // parametri generatora
//                                   double c0, double c1, double r0, double r1, double L0, double L1, double l, double f,// parametri linija + frekvencija
//                                   double C1, double C2, double C3, // parametri kondenzatorske baterije
//                                   Eigen::Matrix3cd Zp) // parametri potrosaca u vidu matrice potrosaca (prakticnije je obaviti unos potrosaca prije poziva solvera)
//    {

    Eigen::VectorXcd solved_vector = NetworkSolver(Eg, Zg, Zg, Zg, c0, c1, r0, r1, x0/w, x1/w, l, f, c, c, c, Zp);

    qInfo() << solved_vector(0,0).real();
    qInfo() << solved_vector(0,0).imag();
    auto arr = solved_vector.array();
    //qInfo() << arr[0].;
}

