/*Author: Marco Montecchi
Department of Energy Technologies
ENEA C.R. Casaccia
Roma - Italy

SoilMirMod is a Software for modeling
multi-angular near-specular reflectance-spectra of soiled solar mirrors


   Copyright (C) 2022  Marco Montecchi

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include "soilmirmod.h"
#include "ui_soilmirmod.h"
#include <iostream>
#include<QtWidgets>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_series_data.h>
#include <qwt_plot_grid.h>
#include <qwt_legend.h>
#include <qwt_plot_marker.h>
#include <qwt_symbol.h>
#include <qwt_interval.h>
#include <QwtPickerMachine>
#include <QwtPlotZoomer>
#include <QPen>
#include <cminpack.h>

using namespace std;

//global variables*****************
int iLock=0;
int Ndat=0,iGph=0,nTheta;
double thetas[8];
double Rmat[500][13];
double ParFit[8][2];
QwtPlot *Gph;
QString pathRoot="/home/marco/Workspace/SoiledMirror";
QString fileExp=pathRoot+"/Results";

//invoked functions ********************************
double fmodel(double wl,double theta,double A,double S,double sigma,double C,double D,double k,double L,double B,double E);
int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag);

struct pointToFit2{
 int NdatFit;
 int Nmin;
 int Nmax;
}pTF2[1];

SoilMirMod::SoilMirMod(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::SoilMirMod)
{
    ui->setupUi(this);

    connect(ui->pushButton_LoadMeas,SIGNAL(clicked()),this, SLOT(LoadMeas()));
    connect(ui->pushButton_BF,SIGNAL(clicked()),this, SLOT(bestFit()));
    connect(ui->doubleSpinBox_S, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_sigma, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_C, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_D, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_k, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_L, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_B, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_E, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_Xmin, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_Xmax, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_Ymin, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
    connect(ui->doubleSpinBox_Ymax, SIGNAL(valueChanged(double)),this, SLOT( calcPlot()));
}

SoilMirMod::~SoilMirMod()
{
    delete ui;
}


void SoilMirMod::LoadMeas(){
    QString fnam= QFileDialog::getOpenFileName(
                this,
                "Choose a R matrix ", //titolo della finestra
                pathRoot, //directory iniziale
                "matrix (*.txt)"); //tipi di file da cercare
    QFile file(fnam);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        cout<<"Attention: error opening "+fnam.toStdString()<< "\n";
        return;
    }
    QTextStream stream(&file);
    QString line;
    line=stream.readLine();
    line =line.simplified();
    ui->lineEdit_measFile->setText(line);
    line=stream.readLine();
    line=line.simplified();
    cout<<line.toStdString()<<"\n";
    QStringList List;
    List =line.split(" ");
    int nV=List.count();
    nTheta=(nV-3)/2;
    printf("nV=%d nTheta=%d\n",nV,nTheta);
    QString pezzo;
    for(int iv=1;iv<nV;iv++){
        pezzo=List.at(iv).toLocal8Bit().constData();
        thetas[iv-1]=pezzo.toDouble();
        printf("thetas[%d]=%f\t",iv-1,thetas[iv-1]);
    }
    line=stream.readLine();
    Ndat=-1;
    do{
        Ndat++;
        //printf("\n");
        for(int i=0;i<nV;i++){
            stream>>Rmat[Ndat][i];
            //printf("Rmat[%d][%d]=%f\t",Ndat,i,Rmat[Ndat][i]);
        }
    }while(!stream.atEnd());
    printf("\nNdat= %d\n",Ndat);
    file.close();
    ui->doubleSpinBox_Xmin->setMinimum(Rmat[0][0]);
    ui->doubleSpinBox_Xmin->setMaximum(Rmat[Ndat-1][0]);
    ui->doubleSpinBox_Xmax->setMinimum(Rmat[0][0]);
    ui->doubleSpinBox_Xmax->setMaximum(Rmat[Ndat-1][0]);
    calcPlot();
}


void SoilMirMod::calcPlot(){
    if(iLock!=0)
        return;
//graph-window initialization
    if(iGph==0){
        Gph=new QwtPlot();
        iGph=1;
        QwtPlotPicker *m_picker1 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                                 QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                                 Gph->canvas() );
    }
    else
        Gph->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
//experimental curves
    QwtPlotCurve *curve1=new QwtPlotCurve("Curve 1");
    QwtPlotCurve *curve2=new QwtPlotCurve("Curve 2");
    QwtPlotCurve *curve3=new QwtPlotCurve("Curve 3");
    QwtPlotCurve *curve4=new QwtPlotCurve("Curve 4");
    double wl[Ndat],dRh[Ndat],dRns1[Ndat],dRns2[Ndat],dRns3[Ndat],dRns1cal[Ndat],dRns2cal[Ndat],dRns3cal[Ndat];
    for(int i=0;i<Ndat;i++){
        wl[i]=Rmat[i][0];
        dRh[i]=(Rmat[i][2]-Rmat[i][1])/Rmat[i][1];
        dRns1[i]=(Rmat[i][4]-Rmat[i][3])/Rmat[i][3];
        Rmat[i][9]=dRh[i];
        Rmat[i][10]=dRns1[i];
        if(nTheta>=2){
            dRns2[i]=(Rmat[i][6]-Rmat[i][5])/Rmat[i][5];
            Rmat[i][11]=dRns2[i];
        }
        if(nTheta==3){
            dRns3[i]=(Rmat[i][8]-Rmat[i][7])/Rmat[i][7];
            Rmat[i][12]=dRns3[i];
        }
    }
    Gph -> setAxisTitle(0,"normalized decrement");
    Gph -> setAxisTitle(2,"wavelength (nm)");
    Gph -> setAutoReplot();
    curve1->setSamples(wl, dRh, Ndat);
    curve1->setPen(QPen(Qt::black,2,Qt::SolidLine));
    curve1->attach(Gph);
    curve2->setSamples(wl, dRns1, Ndat);
    curve2->setPen(QPen(Qt::blue,2,Qt::SolidLine));
    curve2->attach(Gph);
    if(nTheta>=2){
        curve3->setSamples(wl, dRns2, Ndat);
        curve3->setPen(QPen(Qt::green,2,Qt::SolidLine));
        curve3->attach(Gph);
    }
    if(nTheta==3){
        curve4->setSamples(wl, dRns3, Ndat);
        curve4->setPen(QPen(Qt::red,2,Qt::SolidLine));
        curve4->attach(Gph);
    }

//model
    double S=ui->doubleSpinBox_S->value();
    double sigma=ui->doubleSpinBox_sigma->value();
    double C=ui->doubleSpinBox_C->value();
    double D=ui->doubleSpinBox_D->value();
    double k=ui->doubleSpinBox_k->value();
    double L=ui->doubleSpinBox_L->value();
    double B=ui->doubleSpinBox_B->value();
    double E=ui->doubleSpinBox_E->value();
    double Ymin=ui->doubleSpinBox_Ymin->value();
    double Ymax=ui->doubleSpinBox_Ymax->value();
    double WLmin=ui->doubleSpinBox_Xmin->value();
    double WLmax=ui->doubleSpinBox_Xmax->value();
    QwtPlotCurve *curve5=new QwtPlotCurve("Curve 5");
    QwtPlotCurve *curve6=new QwtPlotCurve("Curve 6");
    QwtPlotCurve *curve7=new QwtPlotCurve("Curve 7");
    double WL,A;
    double chi2=0.;
    int Ndafi=0;
    QFile fexp(fileExp);
    if (!fexp.open(QIODevice::WriteOnly | QIODevice::Text)){
        cout<<"Attention: error opening "+fileExp.toStdString()<< endl;
        return;
    }
    QTextStream stream(&fexp);
    for(int i=0;i<Ndat;i++){
        WL=Rmat[i][0];
        A=dRh[i];
        dRns1cal[i]=fmodel(WL,thetas[2],A,S,sigma,C,D,k,L,B,E);
        dRns2cal[i]=fmodel(WL,thetas[4],A,S,sigma,C,D,k,L,B,E);
        dRns3cal[i]=fmodel(WL,thetas[6],A,S,sigma,C,D,k,L,B,E);
        if(WL>=WLmin && WL<=WLmax){
            chi2=chi2+pow(dRns1cal[i]-dRns1[i],2.)
                     +pow(dRns2cal[i]-dRns2[i],2.)
                     +pow(dRns3cal[i]-dRns3[i],2.);
            Ndafi++;
        }
        stream<<WL<<"\t"<<dRns1[i]<<"\t"<<dRns1cal[i];
        if(nTheta>=2)
            stream<<"\t"<<dRns2[i]<<"\t"<<dRns2cal[i];
        if(nTheta==3)
            stream<<"\t"<<dRns3[i]<<"\t"<<dRns3cal[i];
        stream<<"\n";
    }
    chi2=chi2/Ndafi;
    ui->lineEdit_chi2->setText(QString::number(chi2));
    curve5->setSamples(wl, dRns1cal, Ndat);
    curve5->setPen(QPen(Qt::blue,2,Qt::DotLine));
    curve5->attach(Gph);
    if(nTheta>=2){
        curve6->setSamples(wl, dRns2cal, Ndat);
        curve6->setPen(QPen(Qt::green,2,Qt::DotLine));
        curve6->attach(Gph);
    }
    if(nTheta==3){
        curve7->setSamples(wl, dRns3cal, Ndat);
        curve7->setPen(QPen(Qt::red,2,Qt::DotLine));
        curve7->attach(Gph);
    }
    QString title=ui->lineEdit_measFile->text();
    Gph->setTitle(title);
    Gph->setAxisScale(0,Ymin,Ymax,0);
    Gph->setAxisScale(2,Rmat[0][0],Rmat[Ndat-1][0],0);
    Gph->replot();
    Gph->show();

}


void SoilMirMod::bestFit(){
    iLock=1;
    double WLmin=ui->doubleSpinBox_Xmin->value();
    double WLmax=ui->doubleSpinBox_Xmax->value();
    int Ndafit=0;
    int Nmin=Ndat;
    int Nmax=0;
    for(int i=0;i<Ndat;i++){
        if(Rmat[i][0]>=WLmin && Rmat[i][0]<=WLmax){
            Ndafit++;
            Nmin=min(Nmin,i);
            Nmax=max(Nmax,i);
        }
    }
    printf("Ndafit= %d Nmin= %d Nmax= %d\n",Ndafit,Nmin,Nmax);
    Qt::CheckState state;
    for(int i=0; i<7; i++)
        ParFit[i][1]=0.;
    //fit parameters
    ParFit[0][0]=ui->doubleSpinBox_S->value();
    ParFit[1][0]=ui->doubleSpinBox_sigma->value();
    ParFit[2][0]=ui->doubleSpinBox_D->value();
    ParFit[3][0]=ui->doubleSpinBox_k->value();
    ParFit[4][0]=ui->doubleSpinBox_L->value();
    ParFit[5][0]=ui->doubleSpinBox_B->value();
    ParFit[6][0]=ui->doubleSpinBox_E->value();
    ParFit[7][0]=ui->doubleSpinBox_C->value();
    int n=0;
    state = ui->checkBox_S -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[0][1]=1.;
    }
    state = ui->checkBox_sigma -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[1][1]=1.;
    }
    state = ui->checkBox_D -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[2][1]=1.;
    }
    state = ui->checkBox_k -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[3][1]=1.;
    }
    state = ui->checkBox_L -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[4][1]=1.;
    }
    state = ui->checkBox_B -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[5][1]=1.;
    }
    state = ui->checkBox_E -> checkState();
    if(state==Qt::Checked){
        n++;
        ParFit[6][1]=1.;
    }
    int m=Ndafit*nTheta;
    int lwa=m*n+5*n+m;
    int iwa[n];
    double x[n],fvec[m],wa[lwa];
    double tol=sqrt(dpmpar(1));
    int j=0;
    for(int i=0; i<7; i++){
        if(ParFit[i][1]>0.5){
            x[j]=ParFit[i][0];
            printf("Initial x[%d]= %f Nparameter= %d\n",j,x[j],i);
            j++;
        }
    }
    pTF2[0].NdatFit=Ndafit;
    pTF2[0].Nmin=Nmin;
    pTF2[0].Nmax=Nmax;
    int info=0;
    if(m>0 && n>0 && Ndafit>0){
        printf("Launch bestFit with %d parameters and %d data .... ",n,m);
        info=lmdif1(fcn, &pTF2, m, n, x,fvec, tol, iwa, wa, lwa);
    }
    else{
        printf("Nmis=0 || m=0 || n=0 then RETURN!");
        iLock=0;
        return;
    }
    printf(" done with info=%d\n",info);
    j=0;
    for(int i=0; i<7; i++){
        if(ParFit[i][1]>0.5){
            ParFit[i][0]=x[j];
            printf("Final x[%d]= %f Nparameter= %d\n",j,x[j],i);
            j++;
        }
    }
    ui->doubleSpinBox_S    ->setValue(ParFit[0][0]);
    ui->doubleSpinBox_sigma->setValue(ParFit[1][0]);
    ui->doubleSpinBox_D    ->setValue(ParFit[2][0]);
    ui->doubleSpinBox_k    ->setValue(ParFit[3][0]);
    ui->doubleSpinBox_L    ->setValue(ParFit[4][0]);
    ui->doubleSpinBox_B    ->setValue(ParFit[5][0]);
    ui->doubleSpinBox_E    ->setValue(ParFit[6][0]);
    iLock=0;
    calcPlot();
}


double fmodel(double wl,double theta,double A,double S,double sigma,double C,double D,double k,double L,double B,double E){
    double pi=3.14159265359;
    double theRad=theta/180.*pi;
    double Le=L*sqrt(pow(cos(theRad),2.)+pow(E*sin(theRad),2.));
    double y=(A-S*(1-exp(-pow((4*pi*sigma*pow(cos(theRad),C)/wl),2.)))-D/(1+exp(-2*(wl-Le)/k)))/pow(cos(theRad),B);
    return(y);
}


int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag){
    /* calculate the functions at x and return the values in fvec[0] through fvec[m-1] */
    struct pointToFit2 *pTF2 = (struct pointToFit2 *)p;
    int Nmin=pTF2[0].Nmin;
    int Nmax=pTF2[0].Nmax;
    double S,sigma,C,D,k,L,B,E;
    int j=0;
    if(ParFit[0][1]>0.5){
        S=x[j];
        j++;
    }
    else
        S=ParFit[0][0];
    if(ParFit[1][1]>0.5){
        sigma=x[j];
        j++;
    }
    else
        sigma=ParFit[1][0];
    if(ParFit[2][1]>0.5){
        D=x[j];
        j++;
    }
    else
        D=ParFit[2][0];
    if(ParFit[3][1]>0.5){
        k=x[j];
        j++;
    }
    else
        k=ParFit[3][0];
    if(ParFit[4][1]>0.5){
        L=x[j];
        j++;
    }
    else
        L=ParFit[4][0];
    if(ParFit[5][1]>0.5){
        B=x[j];
        j++;
    }
    else
        B=ParFit[5][0];
    if(ParFit[6][1]>0.5){
        E=x[j];
        j++;
    }
    else
        E=ParFit[6][0];
    C=ParFit[7][0];
    int iv=0;
    double WL;
    for(j=Nmin;j<=Nmax;j++){
        WL=Rmat[j][0];
        fvec[iv]=Rmat[j][10]-fmodel(WL,thetas[2],Rmat[j][9],S,sigma,C,D,k,L,B,E);
        iv++;
        if(nTheta>=2){
            fvec[iv]=Rmat[j][11]-fmodel(WL,thetas[4],Rmat[j][9],S,sigma,C,D,k,L,B,E);
            iv++;
        }
        if(nTheta==3){
            fvec[iv]=Rmat[j][12]-fmodel(WL,thetas[6],Rmat[j][9],S,sigma,C,D,k,L,B,E);
            iv++;
        }
    }
    return(0);
}
