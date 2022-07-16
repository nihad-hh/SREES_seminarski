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

//    Eigen::VectorXcd NetworkSolver(std::complex<double> E, std::complex<double> Zg1, std::complex<double> Zg2, std::complex<double> Zg3, // parametri generatora
//                                   double c0, double c1, double r0, double r1, double L0, double L1, double l, double f,// parametri linija + frekvencija
//                                   double C1, double C2, double C3, // parametri kondenzatorske baterije
//                                   Eigen::Matrix3cd Zp) // parametri potrosaca u vidu matrice potrosaca (prakticnije je obaviti unos potrosaca prije poziva solvera)
//    {

//  Eigen::MatrixXcd A(3,2);
//    A << 1,2,3,4,5,6;
//    A.transposeInPlace();
//    Eigen::VectorXcd B;
//    Eigen::VectorXcd Y11_V(9);

//    for(int i = 0; i++; i<9)
//        Y11_V << Zp.data();

//    qInfo() << Y11_V.size();
//    for(int i=0; i<6;i++){
//        B(i) = A.data()[i];
//        qInfo() << B(i).real();
//    }

//VectorXd B(Map<VectorXd>(A.data(), A.cols()*A.rows()));
//tex_creator_vector_I20(I20);

//QString str1 = QString("vX%1").arg(12);

//void create_tex_(){
//    QString file_content = read_File("tex.txt");
//    //QStringList file_content_list = file_content.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);

//    file_content.replace(QString("BROJ1"), QString("12"));
//    file_content.replace("BROJ2", "13");

//    QFile file("tex.tex"); //application current directory

//    if(!file.open(QIODevice::OpenModeFlag::Write
//    {
//        qCritical() << "Could not open file";
//        qCritical() << file.errorString();
//        return;
//    }

//    QTextStream stream(&file);


//    stream << file_content;
//    file.close();
//    return;

//}
