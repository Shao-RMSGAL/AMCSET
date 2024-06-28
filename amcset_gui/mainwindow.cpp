#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include <iostream>
#include <exception>
#include <QCoreApplication>
#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/process.hpp>
#include <filesystem>
#include <string>
#include <cstdlib>
#include <iostream>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    try {
        std::filesystem::path currentPath = std::filesystem::current_path();
        std::cout << "Directory: " << currentPath.string() << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(1);
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_promptServer_clicked()
{
    std::cout << "Hello" << std::endl;
    boost::interprocess::message_queue client_to_server_queue(boost::interprocess::open_only, "client_to_server_message_queue");
    boost::interprocess::message_queue server_to_client_queue(boost::interprocess::open_only, "server_to_client_message_queue");

    std::string sent_message("Prompting server...");
    client_to_server_queue.send(sent_message.c_str(), sent_message.size(), 0);

    char recieved_message[256];
    size_t recieved_size;
    unsigned int priority;
    server_to_client_queue.receive(&recieved_message, sizeof(recieved_message), recieved_size, priority);
    recieved_message[recieved_size] = '\0';

    if(std::string(recieved_message) == "exit") {
        QCoreApplication::quit();
    }

    QString displayMessage(recieved_message);
    ui->textEdit->setPlainText(displayMessage);
    QTextCursor cursor = ui->textEdit->textCursor();
    cursor.movePosition(QTextCursor::Start, QTextCursor::MoveAnchor, 1);
}


void MainWindow::on_exitButton_clicked()
{
    boost::interprocess::message_queue client_to_server_queue(boost::interprocess::open_only, "client_to_server_message_queue");

    std::string exit_message("exit");
    client_to_server_queue.send(exit_message.c_str(), exit_message.size(), 0);
    QCoreApplication::quit();
}
