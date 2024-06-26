#ifndef AMCSET_MAIN_WINDOW_H
#define AMCSET_MAIN_WINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui {
class amcset_main_window;
}
QT_END_NAMESPACE

class amcset_main_window : public QMainWindow
{
    Q_OBJECT

public:
    amcset_main_window(QWidget *parent = nullptr);
    ~amcset_main_window();

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

private:
    Ui::amcset_main_window *ui;
    void promptServer();
};
#endif // AMCSET_MAIN_WINDOW_H
