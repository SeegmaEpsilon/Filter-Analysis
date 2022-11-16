/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Wed 16. Nov 22:35:37 2022
**      by: Qt User Interface Compiler version 4.8.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QSlider>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QGridLayout *gridLayout;
    QPushButton *pushButtonFilter;
    QCustomPlot *canvas_impulse;
    QSlider *horizontalSlider;
    QLineEdit *lineEdit_sineFreq;
    QCustomPlot *canvas_raw;
    QCustomPlot *canvas_filtred;
    QLineEdit *lineEdit_sampleRate;
    QLabel *label_sampleRate;
    QCheckBox *checkBox_magnutideResponse;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1202, 638);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        gridLayout = new QGridLayout(centralWidget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        pushButtonFilter = new QPushButton(centralWidget);
        pushButtonFilter->setObjectName(QString::fromUtf8("pushButtonFilter"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(pushButtonFilter->sizePolicy().hasHeightForWidth());
        pushButtonFilter->setSizePolicy(sizePolicy);

        gridLayout->addWidget(pushButtonFilter, 4, 0, 1, 1);

        canvas_impulse = new QCustomPlot(centralWidget);
        canvas_impulse->setObjectName(QString::fromUtf8("canvas_impulse"));
        canvas_impulse->setStyleSheet(QString::fromUtf8(""));

        gridLayout->addWidget(canvas_impulse, 0, 2, 2, 1);

        horizontalSlider = new QSlider(centralWidget);
        horizontalSlider->setObjectName(QString::fromUtf8("horizontalSlider"));
        sizePolicy.setHeightForWidth(horizontalSlider->sizePolicy().hasHeightForWidth());
        horizontalSlider->setSizePolicy(sizePolicy);
        horizontalSlider->setMinimum(1);
        horizontalSlider->setMaximum(2000);
        horizontalSlider->setValue(50);
        horizontalSlider->setOrientation(Qt::Horizontal);

        gridLayout->addWidget(horizontalSlider, 2, 0, 1, 1);

        lineEdit_sineFreq = new QLineEdit(centralWidget);
        lineEdit_sineFreq->setObjectName(QString::fromUtf8("lineEdit_sineFreq"));
        sizePolicy.setHeightForWidth(lineEdit_sineFreq->sizePolicy().hasHeightForWidth());
        lineEdit_sineFreq->setSizePolicy(sizePolicy);
        lineEdit_sineFreq->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(lineEdit_sineFreq, 3, 0, 1, 1);

        canvas_raw = new QCustomPlot(centralWidget);
        canvas_raw->setObjectName(QString::fromUtf8("canvas_raw"));

        gridLayout->addWidget(canvas_raw, 0, 0, 1, 1);

        canvas_filtred = new QCustomPlot(centralWidget);
        canvas_filtred->setObjectName(QString::fromUtf8("canvas_filtred"));

        gridLayout->addWidget(canvas_filtred, 1, 0, 1, 1);

        lineEdit_sampleRate = new QLineEdit(centralWidget);
        lineEdit_sampleRate->setObjectName(QString::fromUtf8("lineEdit_sampleRate"));
        sizePolicy.setHeightForWidth(lineEdit_sampleRate->sizePolicy().hasHeightForWidth());
        lineEdit_sampleRate->setSizePolicy(sizePolicy);
        lineEdit_sampleRate->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(lineEdit_sampleRate, 3, 2, 1, 1);

        label_sampleRate = new QLabel(centralWidget);
        label_sampleRate->setObjectName(QString::fromUtf8("label_sampleRate"));
        label_sampleRate->setAlignment(Qt::AlignCenter);

        gridLayout->addWidget(label_sampleRate, 2, 2, 1, 1);

        checkBox_magnutideResponse = new QCheckBox(centralWidget);
        checkBox_magnutideResponse->setObjectName(QString::fromUtf8("checkBox_magnutideResponse"));

        gridLayout->addWidget(checkBox_magnutideResponse, 4, 2, 1, 1);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1202, 20));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Filter Analysis", 0, QApplication::UnicodeUTF8));
        pushButtonFilter->setText(QApplication::translate("MainWindow", "\320\237\320\276\321\201\321\202\321\200\320\276\320\270\321\202\321\214", 0, QApplication::UnicodeUTF8));
        lineEdit_sineFreq->setText(QApplication::translate("MainWindow", "50", 0, QApplication::UnicodeUTF8));
        lineEdit_sampleRate->setText(QApplication::translate("MainWindow", "6400", 0, QApplication::UnicodeUTF8));
        label_sampleRate->setText(QApplication::translate("MainWindow", "Sample rate [Hz]", 0, QApplication::UnicodeUTF8));
        checkBox_magnutideResponse->setText(QApplication::translate("MainWindow", "\320\237\320\276\321\201\321\202\321\200\320\276\320\270\321\202\321\214 \320\220\320\247\320\245 (\320\276\320\263\320\270\320\261\320\260\321\216\321\211\320\260\321\217 \321\201\320\277\320\265\320\272\321\202\321\200\320\260)", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
