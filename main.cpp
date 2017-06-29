#include "img_viewer.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    ImageViewer *imageViewer = new ImageViewer();
    imageViewer ->show();
    return a.exec();
}
