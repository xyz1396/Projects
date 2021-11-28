#include "logindialog.h"
#include "ui_logindialog.h"

LoginDialog::LoginDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::LoginDialog)
{
    ui->setupUi(this);
    ui->pushButton->setText("登录到主窗口");
}

LoginDialog::~LoginDialog()
{
    delete ui;
}
