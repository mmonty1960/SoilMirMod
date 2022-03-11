#include "soilmirmod.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    SoilMirMod w;
    w.show();
    return a.exec();
}
