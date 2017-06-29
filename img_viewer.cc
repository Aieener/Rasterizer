#include <QApplication>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QSlider>
#include <QSpinBox>
#include <QDebug>
#include <QMessageBox>

#include "img_viewer.h"

ImageViewer::ImageViewer(QWidget *parent) :
    QMainWindow(parent), channels(0) {

    //------checkboxes ------------------------------
        rCheck = new QCheckBox("&Red", this);
        gCheck = new QCheckBox("&Green", this);
        bCheck = new QCheckBox("&Blue", this);

        grays_scaleCheck = new QCheckBox("&gray",this);
        flipCheck  = new QCheckBox("&flip",this);
        flopCheck  = new QCheckBox("&flop",this);
        transposeCheck  = new QCheckBox("&transpose",this);
        sobelCheck  = new QCheckBox("&sobel",this);

        boxblurCheck  = new QCheckBox("&box_blur");
        medianCheck  = new QCheckBox("&median_blur");
        gaussianCheck  = new QCheckBox("&gaussian_blur");
        scaleCheck  = new QCheckBox("&rescale");
        sizeCheck  = new QCheckBox("&resize");

    //------radiobuttons-----------------------------
        dCheck = new QRadioButton("&diffuse",this);
        dCheck->setChecked(true);
        wCheck = new QRadioButton("&white",this);
        nfCheck = new QRadioButton("&norm_flat",this);
        ngCheck = new QRadioButton("&norm_gouraud",this);
        nbCheck = new QRadioButton("&norm_bary",this);
        ngzCheck = new QRadioButton("&norm_gouraud_z",this);
        nbzCheck = new QRadioButton("&norm_bary_z",this);
    //-----------a special trivial checkbox just to stuck the key------------
        stuckkeyrotate = new QRadioButton("&Rotate",this);

        imgLabel = new QLabel(this);

        QVBoxLayout *layout = new QVBoxLayout();
        setCentralWidget(new QWidget());
        centralWidget()->setLayout(layout);

        layout->addWidget(imgLabel);
       //-------------------------------------------

        createActions();
        createMenus();
        createtoolbar();
        createSpinBoxes();
        createCheckBoxes();
        createDockWindows();
        createButtons();
        //--
        //--

        setWindowTitle(tr("GUI"));
        //show the status
        statusBar();
//        setUnifiedTitleAndToolBarOnMac(true);
//        adjustSize();
        setFixedSize(1600,1000);

}

ImageViewer::~ImageViewer() {
}

// -------- Keypress ----------------------------------------
void ImageViewer::keyPressEvent(QKeyEvent *event)
{
    event = static_cast<QKeyEvent *>(event);
    if(event->key() == Qt::Key_Left)
    {
        qDebug()<<"Rotate left"<<"\n";
        leftrot(5, center, up,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_Right)
    {
        qDebug()<<"Rotate right"<<"\n";
        leftrot(-5, center, up,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_Up)
    {
        qDebug()<<"Rotate up"<<"\n";
        uprot(5, center, up,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_Down)
    {
        qDebug()<<"Rotate down"<<"\n";
        uprot(-5, center, up,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_A)
    {
        qDebug()<<"Trans left"<<"\n";
        translateleft(0.5, center,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_D)
    {
        qDebug()<<"Trans right"<<"\n";
        translateleft(-0.5, center,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_Q)
    {
        qDebug()<<"Trans up"<<"\n";
        translateup(0.5, center,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_E)
    {
        qDebug()<<"Trans down"<<"\n";
        translateup(-0.5, center,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_W)
    {
        qDebug()<<"zoom in"<<"\n";
        zoomin(1.1, center,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_S)
    {
        qDebug()<<"zoom out"<<"\n";
        zoomout(1.1, center,eye);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_1)
    {
        qDebug()<<"increase FOV"<<"\n";
        FOVadjust(0.0001,frustumpara);
        updatefrust();
        raster();
    }
    else if(event->key() == Qt::Key_2)
    {
        qDebug()<<"decrease FOV"<<"\n";
        FOVadjust(-0.0001,frustumpara);
        updatefrust();
        raster();
    }
}

// --------base code slots ----------------------------------
void ImageViewer::redStateChanged(int state) {
    if (state == Qt::Unchecked) {
        channels &= ~RED;
    } else {
        channels |=  RED;
    }

    qDebug() << channels;

    updateLabel();
}

void ImageViewer::greenStateChanged(int state) {
    if (state == Qt::Unchecked) {
        channels &= ~GREEN;
    } else {
        channels |=  GREEN;
    }

    qDebug() << channels;

    updateLabel();
}

void ImageViewer::blueStateChanged(int state) {
    if (state == Qt::Unchecked) {
        channels &= ~BLUE;
    } else {
        channels |=  BLUE;
    }

    qDebug() << channels;

    updateLabel();
}

void ImageViewer::updateLabel() {
    imgLabel->setText(QString::number(channels));
}

// ------------------slots ----------------------------------
void ImageViewer::diffuseStateChanged(){
    roption = DIFFUSE;
}

void ImageViewer::whiteStateChanged(){
    roption = WHITE;
}

void ImageViewer::norm_flatStateChanged(){
    roption = NORM_FLAT;
}

void ImageViewer::norm_gourandStateChanged(){
    roption = NORM_GOURAND;
}

void ImageViewer::norm_baryStateChanged(){
    roption = NORM_BARY;
}

void ImageViewer::norm_gourand_zStateChanged(){
    roption = NORM_GOURAND_Z;
}
void ImageViewer::norm_bary_zStateChanged(){
    roption = NORM_BARY_Z;
}

void ImageViewer::open(){
    QString filename = QFileDialog::getOpenFileName(
                this,
                tr("Open File"),
                "~/",
                "All files (*.*);;ppm files (*.ppm);;text files (*.txt);;obj files (*.obj)"
                );

    QFileInfo fi(filename);
    QString ext = fi.completeSuffix();
    QString bn = fi.baseName();
    if(ext == "ppm"){
        img.load(filename);
        imgLabel->setPixmap(img);
        temp = filename.toStdString();
    }
    else if(ext =="txt"){

        QString mes =bn + " opened";
        QMessageBox::information(this,tr("File Name"), mes);
        //update the cam_opened indicator;
        // for this GUI, if it is a txt file, it has to be the camera file
        // convert qstring to char*
        const char* cfile = filename.toLatin1().data();
        rasterize::Readcamera(cfile,frustumpara,eye,center,up);
        // now the camera info is loaded
        // reset the value on the SpinBoxes:

        cam_change = true;
        eyexSpinBox->setValue(eye[0]);
        eyeySpinBox->setValue(eye[1]);
        eyezSpinBox->setValue(eye[2]);

        centerxSpinBox->setValue(center[0]);
        centerySpinBox->setValue(center[1]);
        centerzSpinBox->setValue(center[2]);

        upxSpinBox->setValue(up[0]);
        upySpinBox->setValue(up[1]);
        upzSpinBox->setValue(up[2]);

        LSpinBox->setValue(frustumpara[0]);
        RSpinBox ->setValue(frustumpara[1]);
        BSpinBox ->setValue(frustumpara[2]);
        TSpinBox ->setValue(frustumpara[3]);
        NSpinBox ->setValue(frustumpara[4]);
        FSpinBox ->setValue(frustumpara[5]);

        //use camerafiles to calc the frustrum matrix
        frus = rasterize::frustrum(frustumpara,eye,center,up);

        cam_opened=1;

        cam_change = false;
    }
    else if(ext =="obj"){
        QString mes =bn + " opened";
        QMessageBox::information(this,tr("File Name"), mes);
        obj_opened = 1;
        //then use tinyobj loader to load the objfile
        const char* cfile = filename.toLatin1().data();
        tinyobj::LoadObj(shapes,materials,cfile,NULL);
    }
    else{
        qDebug()<<"Nothing is opened!";
    }
}

void ImageViewer::save(){
    QString filename = QFileDialog::getSaveFileName(this, tr("Save File"),
                                "~/untitled.ppm",
                                tr("Images (*.ppm)"));
    img.save(filename);
}

void ImageViewer::raster(){
    if( obj_opened == 1){
        img_t *cimg = rasterize::render(imgwitdh,imgheight,shapes,frus,materials,roption);
        write_ppm(cimg,"untitled.ppm");
        img.load("untitled.ppm");
        imgLabel->setPixmap(img);
//        filtfilename = "untitled.ppm";
        temp = "untitled.ppm";
    }
    else{
        qDebug()<<"No obj file has been loaded yet!";
    }
}

void ImageViewer::updatefrust(){
        cam_change = true;
        eyexSpinBox->setValue(eye[0]);
        eyeySpinBox->setValue(eye[1]);
        eyezSpinBox->setValue(eye[2]);

        centerxSpinBox->setValue(center[0]);
        centerySpinBox->setValue(center[1]);
        centerzSpinBox->setValue(center[2]);

        upxSpinBox->setValue(up[0]);
        upySpinBox->setValue(up[1]);
        upzSpinBox->setValue(up[2]);

        LSpinBox->setValue(frustumpara[0]);
        RSpinBox ->setValue(frustumpara[1]);
        BSpinBox ->setValue(frustumpara[2]);
        TSpinBox ->setValue(frustumpara[3]);
        NSpinBox ->setValue(frustumpara[4]);
        FSpinBox ->setValue(frustumpara[5]);

        //use camerafiles to calc the frustrum matrix
        frus = rasterize::frustrum(frustumpara,eye,center,up);

        cam_change = false;
}

void ImageViewer::refrust(){
    if(!cam_change){
        //only refrust it the change is directly from the spin box
        // do not refrust it if the change is due to camera change
        eye[0] = eyexSpinBox->value();
        eye[1] = eyeySpinBox->value();
        eye[2] = eyezSpinBox->value();

        center[0] = centerxSpinBox->value();
        center[1] = centerySpinBox->value();
        center[2] = centerzSpinBox->value();

        up[0] = upxSpinBox->value();
        up[1] = upySpinBox->value();
        up[2] = upzSpinBox->value();

        frustumpara[0] = LSpinBox->value();
        frustumpara[1] = RSpinBox->value();
        frustumpara[2] = BSpinBox->value();
        frustumpara[3] = TSpinBox->value();
        frustumpara[4] = NSpinBox->value();
        frustumpara[5] = FSpinBox->value();
        frus = rasterize::frustrum(frustumpara,eye,center,up);
    }
}

void ImageViewer::rasresize(){
    imgwitdh = WidthBox->value();
    imgheight = HeightBox->value();
}

void ImageViewer::filter_para(){
    nbblur = bblurBox->value();
    nmblur = medianBox->value();
    filtwidth = filtwBox->value();
    filtheight = filthBox->value();
    ngblur = ngaussianBox->value();
    sgblur = sgaussianBox->value();
    filtscale = scaleBox->value();
}

void ImageViewer::filter(){

    // read in the input image
    filtfilename = temp.c_str();
    img_t *cimg = read_ppm(filtfilename);

    //check if there is a tag of "-grayscale"
    if(grays_scaleCheck->isChecked()){
        filters::grays_scale(cimg);
    }
    //check if there is a tag of "-flipimage"
    if(flipCheck->isChecked()){
        filters::flipimage(cimg);
    }
    //check if there is a tag of "-flopimage"
    if(flopCheck->isChecked()){
        filters::flopimage(cimg);
    }
    //check if there is a tag of "-transpose"
    if(transposeCheck->isChecked()){
        filters::transpose(cimg);
    }
    //check if there is a tag of "-boxblur"
    if(boxblurCheck->isChecked()){
        filters::boxblur(cimg,nbblur);
    }
    //check if there is a tag of "-median"
    if(medianCheck->isChecked()){
        filters::median(cimg,nmblur);
    }
    //check if there is a tag of "-gaussian"
    if(gaussianCheck->isChecked()){
        filters::gaussian(cimg,ngblur,sgblur);
    }
    //check if there is a tag of "-sobel"
    if(sobelCheck->isChecked()){
        filters::sobel(cimg);
    }

    //check if there is a tag of "-scale"
    if(scaleCheck->isChecked()&&!sizeCheck->isChecked()){
        int scalewidth, scaleheight;
        scalewidth = (int)(filtscale * cimg->w);
        scaleheight = (int)(filtscale * cimg->h);
        filters::resize(cimg,scalewidth,scaleheight);
    }
    //check if there is a tag of "-size"
    else if(sizeCheck->isChecked()&&!scaleCheck->isChecked()){
        filters::resize(cimg,filtwidth,filtheight);
    }

    //the case if both size and scale is checked
    //then rescale the img based on the resized img
    else if(sizeCheck->isChecked()&&scaleCheck->isChecked()){
        int scalewidth, scaleheight;
        filters::resize(cimg,filtwidth,filtheight);
        scalewidth = (int)(filtscale * filtwidth);
        scaleheight = (int)(filtscale * filtheight);
        filters::resize(cimg,scalewidth,scaleheight);
    }

    //format_RGB888
    // now write the image
    write_ppm(cimg, "filuntitled.ppm");
    img.load("filuntitled.ppm");
    imgLabel->setPixmap(img);
    filtfilename = "filuntitled.ppm";

    // free up memory
    destroy_img(&cimg); // &img is the address in memory where the img variable is stored
}


void ImageViewer::leftrot(float degrees, vec4&center, vec4&up,const vec4&eye){
    vec4 forwar = center - eye;
    // the lenth of forward should not change
    double flenth = forwar.length();

    mat4 rotmatrix;
    rotmatrix = rotmatrix.rot(degrees,up[0],up[1],up[2]);

    forwar = rotmatrix*forwar;
    double nflenth = forwar.length();
    forwar[0] *= flenth/nflenth;
    forwar[1] *= flenth/nflenth;
    forwar[2] *= flenth/nflenth;
    forwar[3] *= flenth/nflenth;
    center = forwar + eye;
}


void ImageViewer::uprot(float degrees, vec4&center, vec4&up, const vec4&eye){

    vec4 forward = center - eye;
    double flenth = forward.length();

    vec4 right = cross(forward,up);
    right.norm();

    mat4 rotmatrix;
    rotmatrix = rotmatrix.rot(degrees,right[0],right[1],right[2]);

    up = rotmatrix*up;
    forward = rotmatrix*forward;
    double nflenth = forward.length();
    forward[0] *= flenth/nflenth;
    forward[1] *= flenth/nflenth;
    forward[2] *= flenth/nflenth;
    forward[3] *= flenth/nflenth;

    center = forward + eye;
}

void ImageViewer::translateleft(float dis, vec4&center,vec4&eye){
    mat4 translate;
    translate = translate.trans(dis,0,0);
    vec4 forward = center - eye;
    eye = translate*eye;
    center = forward + eye;

}

void ImageViewer::translateup(float dis, vec4&center,vec4 &eye ){
    mat4 translate;
    translate = translate.trans(0,dis,0);
    vec4 forward = center - eye;
    eye = translate*eye;
    center = forward + eye;

}

void ImageViewer::zoomin(float dis, vec4 &center, vec4 &eye){
    // eye should move towards to center
    mat4 translate;
    vec4 forward = center - eye;
    double flenth = forward.length();
    forward.norm();

    translate = translate.trans(dis*forward[0],dis*forward[1],dis*forward[2]);

    if(flenth !=0){
        eye = translate*eye;
    }
    else{
        qDebug()<<"Can't zoom in anyfurther!\n";
    }
}
void ImageViewer::zoomout(float dis, vec4 &center, vec4 &eye){
    // eye should move away from center
    mat4 translate;
    vec4 forward =  eye-center;
    forward.norm();
    translate = translate.trans(dis*forward[0],dis*forward[1],dis*forward[2]);
    eye = translate*eye;
}
void ImageViewer::FOVadjust(float degrees,std::array<float,6> &frustumpara ){
    frustumpara[0] -=degrees;
    frustumpara[1] +=degrees;
    frustumpara[2] -=degrees;
    frustumpara[3] +=degrees;
}

// ------------------create things--------------------------
void ImageViewer::createCheckBoxes(){

    //------------Rasterize sections----------------------------
    checkBoxesGroup = new QGroupBox();
    checkBoxLayout = new QGridLayout;

    //--- add the raster options---
    checkBoxLayout->addWidget(dCheck,0,0);
    connect(dCheck, SIGNAL(toggled(bool)),this, SLOT(diffuseStateChanged()));
    checkBoxLayout->addWidget(wCheck,0,1);
    connect(wCheck, SIGNAL(toggled(bool)),this, SLOT(whiteStateChanged()));
    checkBoxLayout->addWidget(nfCheck,0,2);
    connect(nfCheck, SIGNAL(toggled(bool)),this, SLOT(norm_flatStateChanged()));
    checkBoxLayout->addWidget(ngCheck,1,0);
    connect(ngCheck, SIGNAL(toggled(bool)),this, SLOT(norm_gourandStateChanged()));
    checkBoxLayout->addWidget(nbCheck,1,1);
    connect(nbCheck, SIGNAL(toggled(bool)),this, SLOT(norm_baryStateChanged()));
    checkBoxLayout->addWidget(ngzCheck,1,2);
    connect(ngzCheck, SIGNAL(toggled(bool)),this, SLOT(norm_gourand_zStateChanged()));
    checkBoxLayout->addWidget(nbzCheck,1,3);
    connect(nbzCheck, SIGNAL(toggled(bool)),this, SLOT(norm_bary_zStateChanged()));


    //--- add the size options ----
    QLabel *widthlabel = new QLabel(tr("Width of img"));
    QLabel *heightlabel = new QLabel(tr("Height of img"));
    checkBoxLayout->addWidget(widthlabel,0,4);
    checkBoxLayout->addWidget(WidthBox,0,5);
    checkBoxLayout->addWidget(heightlabel,1,4);
    checkBoxLayout->addWidget(HeightBox,1,5);

    checkBoxesGroup->setLayout(checkBoxLayout);

    //------------Filters   sections----------------------------
    dcheckBoxesGroup = new QGroupBox();
    dcheckBoxLayout = new QGridLayout();

    dcheckBoxLayout->addWidget(grays_scaleCheck ,0,0);
    dcheckBoxLayout->addWidget(flipCheck ,0,1);
    dcheckBoxLayout->addWidget(flopCheck ,1,0);
    dcheckBoxLayout->addWidget(transposeCheck ,1,1);
    dcheckBoxLayout->addWidget(sobelCheck ,2,0);


    QLabel *bblurlabel = new QLabel(tr("n of box_blur"));
    dcheckBoxLayout->addWidget(boxblurCheck ,0,3);
    dcheckBoxLayout->addWidget(bblurlabel,1,3);
    dcheckBoxLayout->addWidget(bblurBox,2,3);

    QLabel *medianlabel = new QLabel(tr("n of median_blur"));
    dcheckBoxLayout->addWidget(medianCheck ,0,4);
    dcheckBoxLayout->addWidget(medianlabel,1,4);
    dcheckBoxLayout->addWidget(medianBox,2,4);


    QLabel *ngaulabel = new QLabel(tr("n of gaussian_blur"));
    QLabel *sgaulabel = new QLabel(tr("s of gaussian_blur"));
    dcheckBoxLayout->addWidget(gaussianCheck ,0,5);
    dcheckBoxLayout->addWidget(ngaulabel,1,5);
    dcheckBoxLayout->addWidget(sgaulabel,3,5);
    dcheckBoxLayout->addWidget(ngaussianBox,2,5);
    dcheckBoxLayout->addWidget(sgaussianBox,4,5);

    QLabel *fwlabel = new QLabel(tr("width for resize"));
    QLabel *fhlabel = new QLabel(tr("height for resize"));
    dcheckBoxLayout->addWidget(sizeCheck ,0,6);
    dcheckBoxLayout->addWidget(fwlabel,1,6);
    dcheckBoxLayout->addWidget(filtwBox,2,6);
    dcheckBoxLayout->addWidget(fhlabel,3,6);
    dcheckBoxLayout->addWidget(filthBox,4,6);

    QLabel *scalelabel = new QLabel(tr("scale for rescale"));
    dcheckBoxLayout->addWidget(scaleCheck ,0,7);
    dcheckBoxLayout->addWidget(scalelabel,1,7);
    dcheckBoxLayout->addWidget(scaleBox,2,7);

    dcheckBoxesGroup->setLayout(dcheckBoxLayout);
}

void ImageViewer::createMenus(){
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openAct);
    fileMenu->addAction(saveAct);
}

void ImageViewer::createtoolbar(){
    fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(openAct);
    fileToolBar->addAction(saveAct);

}

void ImageViewer::createActions(){
    //open Act
    const QIcon openIcon = QIcon::fromTheme("document-open",QIcon("open.png"));
    openAct = new QAction(tr("&Open..."), this);
    openAct = new QAction(openIcon, tr("&Open..."), this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip(tr("Open an existing file"));
    connect(openAct, &QAction::triggered, this, &ImageViewer::open);

    //save Act
    const QIcon saveIcon = QIcon::fromTheme("document-save", QIcon("save.png"));
    saveAct = new QAction(saveIcon, tr("&Save"), this);
    saveAct->setShortcuts(QKeySequence::Save);
    saveAct->setStatusTip(tr("Save the img to disk"));
    connect(saveAct, &QAction::triggered, this, &ImageViewer::save);

    //raster Act
    rasterAct = new QAction(this);
    connect(rasterAct, &QAction::triggered, this, &ImageViewer::raster);

    //filter Act
    filterAct = new QAction(this);
    connect(filterAct, &QAction::triggered, this, &ImageViewer::filter);
}

void ImageViewer::createSpinBoxes(){
//    spinBoxesGroup = new QGroupBox(tr("Spinboxes"));
    spinBoxesGroup = new QGroupBox();

    //------------eyex----------------------------
    QLabel *eyexlabel = new QLabel(tr("Enter camera eye_x"));
    eyexSpinBox = new QDoubleSpinBox;
    eyexSpinBox ->setRange(-1000,1000);
    eyexSpinBox ->setSingleStep(1.0);
    eyexSpinBox ->setValue(eye[0]);
    connect(eyexSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    //------------eyey----------------------------
    QLabel *eyeylabel = new QLabel(tr("Enter camera eye_y"));
    eyeySpinBox = new QDoubleSpinBox;
    eyeySpinBox ->setRange(-1000,1000);
    eyeySpinBox ->setSingleStep(1.0);
    eyeySpinBox ->setValue(eye[1]);
    connect(eyeySpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    //------------eyez----------------------------
    QLabel *eyezlabel = new QLabel(tr("Enter camera eye_z"));
    eyezSpinBox = new QDoubleSpinBox;
    eyezSpinBox ->setRange(-1000,1000);
    eyezSpinBox ->setSingleStep(1.0);
    eyezSpinBox ->setValue(eye[2]);
    connect(eyezSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

     //------------centerx-------------------------
    QLabel *centerxlabel = new QLabel(tr("Enter camera center_x"));
    centerxSpinBox = new QDoubleSpinBox;
    centerxSpinBox ->setRange(-1000,1000);
    centerxSpinBox ->setSingleStep(1.0);
    centerxSpinBox ->setValue(center[0]);
    connect(centerxSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    //------------centery-------------------------
    QLabel *centerylabel = new QLabel(tr("Enter camera center_y"));
    centerySpinBox = new QDoubleSpinBox;
    centerySpinBox  ->setRange(-1000,1000);
    centerySpinBox  ->setSingleStep(1.0);
    centerySpinBox ->setValue(center[1]);
    connect(centerySpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    //------------centerz-------------------------
    QLabel *centerzlabel = new QLabel(tr("Enter camera center_z"));
    centerzSpinBox = new QDoubleSpinBox;
    centerzSpinBox ->setRange(-1000,1000);
    centerzSpinBox ->setSingleStep(1.0);
    centerzSpinBox ->setValue(center[2]);
    connect(centerzSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    //-------------upx----------------------------
    QLabel *upxlabel = new QLabel(tr("Enter camera up_x"));
    upxSpinBox = new QDoubleSpinBox;
    upxSpinBox ->setRange(-1000,1000);
    upxSpinBox ->setSingleStep(1.0);
    upxSpinBox ->setValue(up[0]);
    connect(upxSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    //------------up_y----------------------------
    QLabel *upylabel = new QLabel(tr("Enter camera up_y"));
    upySpinBox = new QDoubleSpinBox;
    upySpinBox ->setRange(-1000,1000);
    upySpinBox  ->setSingleStep(1.0);
    upySpinBox  ->setValue(up[1]);
    connect(upySpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    //------------up_z----------------------------
    QLabel *upzlabel = new QLabel(tr("Enter camera up_z"));
    upzSpinBox = new QDoubleSpinBox;
    upzSpinBox ->setRange(-1000,1000);
    upzSpinBox ->setSingleStep(1.0);
    upzSpinBox ->setValue(up[2]);
    connect(upzSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    //------------L R B T N F--------------------
    QLabel *Llabel = new QLabel(tr("Enter camera left"));
    LSpinBox = new QDoubleSpinBox;
    LSpinBox ->setRange(-1000,1000);
    LSpinBox ->setSingleStep(1.0);
    LSpinBox ->setValue(frustumpara[0]);
    LSpinBox ->setDecimals(15);
    connect(LSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    QLabel *Rlabel = new QLabel(tr("Enter camera right"));
    RSpinBox = new QDoubleSpinBox;
    RSpinBox ->setRange(-1000,1000);
    RSpinBox ->setSingleStep(1.0);
    RSpinBox ->setValue(frustumpara[1]);
    RSpinBox ->setDecimals(15);
    connect(RSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    QLabel *Blabel = new QLabel(tr("Enter camera bottom"));
    BSpinBox = new QDoubleSpinBox;
    BSpinBox ->setRange(-1000,1000);
    BSpinBox ->setSingleStep(1.0);
    BSpinBox ->setValue(frustumpara[2]);
    BSpinBox ->setDecimals(15);
    connect(BSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    QLabel *Tlabel = new QLabel(tr("Enter camera top"));
    TSpinBox = new QDoubleSpinBox;
    TSpinBox ->setRange(-1000,1000);
    TSpinBox ->setSingleStep(1.0);
    TSpinBox ->setValue(frustumpara[3]);
    TSpinBox ->setDecimals(15);
    connect(TSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    QLabel *Nlabel = new QLabel(tr("Enter camera near"));
    NSpinBox = new QDoubleSpinBox;
    NSpinBox ->setRange(-1000,1000);
    NSpinBox ->setSingleStep(1.0);
    NSpinBox ->setValue(frustumpara[4]);
    NSpinBox ->setDecimals(15);
    connect(NSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));

    QLabel *Flabel = new QLabel(tr("Enter camera far"));
    FSpinBox = new QDoubleSpinBox;
    FSpinBox ->setRange(-1000,1000);
    FSpinBox ->setSingleStep(1.0);
    FSpinBox ->setValue(frustumpara[5]);
    FSpinBox ->setDecimals(15);
    connect(FSpinBox,SIGNAL(valueChanged(double)),this, SLOT(refrust()));



    //--------------------------------------------
    spinBoxLayout = new QGridLayout;
    spinBoxLayout->addWidget(Llabel,0,0);
    spinBoxLayout->addWidget(LSpinBox,1,0);

    spinBoxLayout->addWidget(Rlabel,0,1);
    spinBoxLayout->addWidget(RSpinBox,1,1);

    spinBoxLayout->addWidget(Blabel,0,2);
    spinBoxLayout->addWidget(BSpinBox,1,2);

    spinBoxLayout->addWidget(Tlabel,0,3);
    spinBoxLayout->addWidget(TSpinBox,1,3);

    spinBoxLayout->addWidget(Nlabel,0,4);
    spinBoxLayout->addWidget(NSpinBox,1,4);

    spinBoxLayout->addWidget(Flabel,0,5);
    spinBoxLayout->addWidget(FSpinBox,1,5);

    spinBoxLayout->addWidget(eyexlabel,2,0);
    spinBoxLayout->addWidget(eyexSpinBox,3,0);

    spinBoxLayout->addWidget(eyeylabel,2,1);
    spinBoxLayout->addWidget(eyeySpinBox,3,1);

    spinBoxLayout->addWidget(eyezlabel,2,2);
    spinBoxLayout->addWidget(eyezSpinBox,3,2);

    //---
    spinBoxLayout->addWidget(centerxlabel,2,3);
    spinBoxLayout->addWidget(centerxSpinBox,3,3);

    spinBoxLayout->addWidget(centerylabel,2,4);
    spinBoxLayout->addWidget(centerySpinBox,3,4);

    spinBoxLayout->addWidget(centerzlabel,2,5);
    spinBoxLayout->addWidget(centerzSpinBox,3,5);

    //---
    spinBoxLayout->addWidget(upxlabel,4,0);
    spinBoxLayout->addWidget(upxSpinBox,5,0);

    spinBoxLayout->addWidget(upylabel,4,1);
    spinBoxLayout->addWidget(upySpinBox,5,1);

    spinBoxLayout->addWidget(upzlabel,4,2);
    spinBoxLayout->addWidget(upzSpinBox,5,2);

    QLabel *rotatelabel = new QLabel(tr("Click me then press keyboard to rotate"));
    spinBoxLayout->addWidget(rotatelabel,4,4);
    spinBoxLayout->addWidget(stuckkeyrotate,5,4);
    //------------------------------------
    spinBoxesGroup->setLayout(spinBoxLayout);
    //--- set rastersize ----------------------------------------

    WidthBox = new QDoubleSpinBox;
    WidthBox->setRange(0,1024);
    WidthBox->setSingleStep(8);
    WidthBox->setValue(imgwitdh);
    connect(WidthBox,SIGNAL(valueChanged(double)),this, SLOT(rasresize()));

    HeightBox = new QDoubleSpinBox;
    HeightBox ->setRange(0,1024);
    HeightBox ->setSingleStep(8);
    HeightBox ->setValue(imgheight);
    connect(HeightBox,SIGNAL(valueChanged(double)),this, SLOT(rasresize()));

    //--- set blur parameters -----------------------------------
    bblurBox = new QSpinBox;
    bblurBox->setRange(0,100);
    bblurBox->setSingleStep(1);
    bblurBox->setValue(nbblur);
    connect(bblurBox,SIGNAL(valueChanged(int)),this, SLOT(filter_para()));

    medianBox = new QSpinBox;
    medianBox->setRange(0,100);
    medianBox->setSingleStep(1);
    medianBox->setValue(nmblur);
    connect(medianBox,SIGNAL(valueChanged(int)),this, SLOT(filter_para()));

    ngaussianBox = new QDoubleSpinBox;
    ngaussianBox->setRange(0,100);
    ngaussianBox->setSingleStep(1);
    ngaussianBox->setValue(ngblur);
    connect(ngaussianBox,SIGNAL(valueChanged(double)),this, SLOT(filter_para()));

    sgaussianBox = new QDoubleSpinBox;
    sgaussianBox->setRange(0,100);
    sgaussianBox->setSingleStep(1);
    sgaussianBox->setValue(sgblur);
    connect(sgaussianBox,SIGNAL(valueChanged(double)),this, SLOT(filter_para()));

    //--- set resize/rescale parameters -----------------------------------
    filtwBox = new QSpinBox;
    filtwBox->setRange(0,1024);
    filtwBox->setSingleStep(8);
    filtwBox->setValue(filtwidth);
    connect(filtwBox,SIGNAL(valueChanged(int)),this, SLOT(filter_para()));

    filthBox = new QSpinBox;
    filthBox->setRange(0,1024);
    filthBox->setSingleStep(8);
    filthBox->setValue(filtheight);
    connect(filthBox,SIGNAL(valueChanged(int)),this, SLOT(filter_para()));

    scaleBox = new QDoubleSpinBox;
    scaleBox->setRange(-50,50);
    scaleBox->setSingleStep(0.5);
    scaleBox->setValue(filtscale);
    connect(scaleBox,SIGNAL(valueChanged(double)),this, SLOT(filter_para()));
}

void ImageViewer::createDockWindows(){
    // the camera dock
    QDockWidget *dock_cam = new QDockWidget(tr("Camera"),this);
//    dock_cam->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    dock_cam ->setWidget(spinBoxesGroup);
    addDockWidget(Qt::TopDockWidgetArea, dock_cam);
    //add memu action
    fileMenu->addAction(dock_cam->toggleViewAction());

    // the filter option dock
    QDockWidget *dock_filter = new QDockWidget(tr("Filters"),this);
    dock_filter->setWidget(dcheckBoxesGroup);
    addDockWidget(Qt::BottomDockWidgetArea, dock_filter);
    fileMenu->addAction(dock_filter->toggleViewAction());

    // the raster options dock
    QDockWidget *dock_rast = new QDockWidget(tr("Raster Options"),this);
    dock_rast->setWidget(checkBoxesGroup);
    addDockWidget(Qt::BottomDockWidgetArea, dock_rast);
    fileMenu->addAction(dock_rast->toggleViewAction());

    tabifyDockWidget(dock_rast, dock_filter);
    dock_rast->raise();

}

void ImageViewer::createButtons(){
    QPushButton *button = new QPushButton("Rasterize it!");
    connect(button,SIGNAL(clicked()), rasterAct, SLOT(trigger()));
    checkBoxLayout->addWidget(button,1,9);

    QPushButton *filtbutton = new QPushButton("Filter it!");
    connect(filtbutton,SIGNAL(clicked()), filterAct, SLOT(trigger()));
    dcheckBoxLayout->addWidget(filtbutton,4,7);
}
