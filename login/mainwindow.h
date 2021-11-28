#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void newFile();
    bool maybeSave();
    bool save();
    bool saveAs();
    bool saveFile(const QString &fileName);
    bool loadFile(const QString &fileName);

private slots:
    void on_action_N_triggered();

    void on_action_S_triggered();

    void on_action_A_triggered();

    void on_action_Open_triggered();

    void on_action_Close_triggered();

    void on_action_Exit_triggered();

    void on_action_Undo_triggered();

    void on_action_Cut_triggered();

    void on_action_Copy_triggered();

    void on_action_Paste_triggered();

private:
    Ui::MainWindow *ui;
    bool isUntitled;
    QString curFile;

protected:
    void closeEvent(QCloseEvent *event);

};

#endif // MAINWINDOW_H
