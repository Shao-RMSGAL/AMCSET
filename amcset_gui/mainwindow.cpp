#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <iostream>
#include <exception>
#include <QCoreApplication>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    std::cout << "Hello!" << std::endl;
}


void MainWindow::on_exitButton_clicked()
{
    std::terminate();
}

