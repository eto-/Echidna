#include <QApplication>
#include "dbMon.h"        // constants and classes common for the whole project are here
#include "mainWindow.h"

int main(int argc, char *argv[])
{
    //streamsRedirection();

    QApplication app(argc, argv);

    mainWindow mainWin;
    mainWin.show();

    return app.exec();
}
