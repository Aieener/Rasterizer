#ifndef IMG_VIEWER_H
#define IMG_VIEWER_H

#include "ppm_filters.h"
#include "rasterize.h"
#include <QMainWindow>
#include <QCheckBox>
#include <QLabel>
#include <QPixmap>
#include <QtWidgets>
#include <QFileDialog>
#include <QImage>
#include <QHBoxLayout>
#include <QSlider>
#include <QSpinBox>
#include <QPushButton>

class ImageViewer : public QMainWindow
{
    Q_OBJECT

    public:
        enum {RED = 1, GREEN = 2, BLUE = 3};
        enum {DIFFUSE,WHITE,NORM_FLAT,NORM_GOURAND,NORM_BARY,NORM_GOURAND_Z,NORM_BARY_Z};
        enum {GREY,FLIP,FLOP,TRANS,SOBEL,BBLUR,MBLUR,GBLUR,SIZE,SCALE};
        explicit ImageViewer(QWidget *parent = 0);
        virtual ~ImageViewer();
        //----------------------------
        bool loadFile(const QString &);
        //----------------------------

//    protected:

//        virtual void keyPressEvent(QKeyEvent *event);

    private slots:
        //--------- open and save slots ----------------
        void open();
        void save();
        //---------state changed slots -----------------
        void diffuseStateChanged();
        void whiteStateChanged();
        void norm_flatStateChanged();
        void norm_gourandStateChanged();
        void norm_baryStateChanged();
        void norm_gourand_zStateChanged();
        void norm_bary_zStateChanged();

        //--------function slots------------------------
        void raster();
        void refrust();
        void rasresize();
        void filter();
        void filter_para();
        void updatefrust();
        // ----------------------------------------------

        //------------from basecode----------------------
        void redStateChanged(int state);
        void greenStateChanged(int state);
        void blueStateChanged(int state);
        // ----------------------------------------------

        //------extra credit-----------------------------
        void leftrot(float degrees, vec4&center, vec4&up ,const vec4&eye);
        void uprot(float degrees, vec4&center, vec4&up ,const vec4&eye);
        void translateleft(float dis, vec4&center, vec4&eye);
        void translateup(float dis, vec4&center,vec4 &eye);
        void zoomin(float dis, vec4&center,vec4 &eye);
        void zoomout(float dis, vec4&center,vec4 &eye);
        void FOVadjust(float degrees,std::array<float,6> &frustumpara );

    private:
        //------------from basecode----------------------
        void updateLabel(void);

        virtual void keyPressEvent(QKeyEvent *event);
        int channels;
        QCheckBox *rCheck, *gCheck, *bCheck;
        QCheckBox *grays_scaleCheck , *flipCheck , *flopCheck , *transposeCheck , *sobelCheck ;
        QCheckBox *boxblurCheck , *medianCheck , *gaussianCheck , *scaleCheck , *sizeCheck;
        QRadioButton *dCheck,*wCheck, *nfCheck, *ngCheck,*nbCheck,*ngzCheck,*nbzCheck, *stuckkeyrotate ;
        QLabel *imgLabel;
        QPixmap img;
        // ------------------------------------
        // following the documentation
        // from Qt's Image Viewer/manu Class Example
        // ------------------------------------
        void createActions();
        void createMenus();
        // create dockwindow
        void createDockWindows();
        // create spinbox for camera
        void createSpinBoxes();
        // create butto;
        void createButtons();
        // create checkboxes
        void createCheckBoxes();
        // create toolbar
        void createtoolbar();


        //--spinbox elements
        QGridLayout *spinBoxLayout;
        QGridLayout *checkBoxLayout;
        QGridLayout  *dcheckBoxLayout;

        QDoubleSpinBox *eyexSpinBox;
        QDoubleSpinBox *eyeySpinBox;
        QDoubleSpinBox *eyezSpinBox;
        QDoubleSpinBox *centerxSpinBox;
        QDoubleSpinBox *centerySpinBox;
        QDoubleSpinBox *centerzSpinBox;
        QDoubleSpinBox *upxSpinBox;
        QDoubleSpinBox *upySpinBox;
        QDoubleSpinBox *upzSpinBox;
        QDoubleSpinBox *LSpinBox;
        QDoubleSpinBox *RSpinBox;
        QDoubleSpinBox *BSpinBox;
        QDoubleSpinBox *TSpinBox;
        QDoubleSpinBox *NSpinBox;
        QDoubleSpinBox *FSpinBox;
        QDoubleSpinBox *WidthBox;
        QDoubleSpinBox *HeightBox;
        QDoubleSpinBox *ngaussianBox;
        QDoubleSpinBox *sgaussianBox;
        QDoubleSpinBox *scaleBox;
        QSpinBox *bblurBox;
        QSpinBox *medianBox;
        QSpinBox *filtwBox;
        QSpinBox *filthBox;


        QGroupBox *spinBoxesGroup;   //cam
        QGroupBox *checkBoxesGroup;  //raster
        QGroupBox *dcheckBoxesGroup; //filter

        //---
        //----------------------------------
        QMenu *fileMenu;
        QToolBar *fileToolBar;
        QAction *openAct;
        QAction *saveAct;
        QAction *rasterAct;
        QAction *filterAct;

        QLabel *infoLabel;
        //----------------------------------
        // --- opened indicator-------------
        bool cam_change = false;
        int cam_opened = 0;
        int obj_opened = 0;
        // --camera parameters-----------
        vec4 eye = vec4(0,0,0,1);
        vec4 center=vec4(0,0,0,1);
        vec4 up=vec4(0,0,0,1);

        std::array<float,6> frustumpara;
        std::array<mat4,3> frus;
        //---img size parameter----------
        int imgwitdh = 512; //default value to be 512
        int imgheight = 512; //default value to be 512
        // ---tinyobj files--------------
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;
        int roption = 0;

        // ---- filter parameters -------
        const char *filtfilename;
        std::string temp;
        int nbblur = 0;
        int nmblur = 0;
        int filtwidth = 512;
        int filtheight = 512;
        double ngblur = 0;
        double sgblur = 1;
        double filtscale = 1;
};

#endif /* IMG_VIEWER_H */
