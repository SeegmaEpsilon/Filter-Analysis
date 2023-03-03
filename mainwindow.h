#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qcustomplot.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void sweep(double f_start, double f_end, double interval, int n_steps);

private slots:
    void on_pushButtonFilter_clicked();

    void on_horizontalSlider_valueChanged(int value);

    void slotMouseMove(QMouseEvent * event);
    void slotMousePress(QMouseEvent* event);
    void slotMouseRelease(QMouseEvent* event);

private:
    Ui::MainWindow *ui;
    void plotGrapics();
    void plotMagnitudeResponse();
    QVector<double> magnitudeResponseMaximumes;
    QVector<double> sweepValues;
    QVector<double> x_all, x_3db, y_3db;

    bool mouse_hold;
    QCPRange DragStartHorzRange;
    QCPRange DragStartVertRange;
    QPoint mDragStart;
};

#endif // MAINWINDOW_H
