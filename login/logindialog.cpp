#include "logindialog.h"
#include "ui_logindialog.h"
#include <QMessageBox>
#include <QLineEdit>

loginDialog::loginDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::loginDialog)
{
    ui->setupUi(this);
    this->setWindowTitle("Login");
    ui->pwdLineEdit->setEchoMode(QLineEdit::Password);
    ui->userLineEdit->setPlaceholderText("请输入账号");
    ui->pwdLineEdit->setPlaceholderText("请输入密码");
}



loginDialog::~loginDialog()
{
    delete ui;
}

void loginDialog::on_loginButton_clicked()
{
    if(ui->userLineEdit->text().trimmed() == tr("xyz") &&
            ui->pwdLineEdit->text() == tr("1234"))
    {
        accept();
    }
    else
    {
        QMessageBox::warning(this,tr("Warning"),
                             tr("User name or password error!"),
                             QMessageBox::Yes);
        ui->userLineEdit->clear();
        ui->pwdLineEdit->clear();
        ui->userLineEdit->setFocus();
    }
}
