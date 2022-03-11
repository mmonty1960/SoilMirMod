#ifndef SOILMIRMOD_H
#define SOILMIRMOD_H

#include <QWidget>

QT_BEGIN_NAMESPACE
namespace Ui { class SoilMirMod; }
QT_END_NAMESPACE

class SoilMirMod : public QWidget
{
    Q_OBJECT

public:
    SoilMirMod(QWidget *parent = nullptr);
    ~SoilMirMod();

private:
    Ui::SoilMirMod *ui;

public slots:
    void LoadMeas();
    void calcPlot();
    void bestFit();
};
#endif // SOILMIRMOD_H
