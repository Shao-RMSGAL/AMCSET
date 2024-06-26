#include "amcset_main_window.h"
#include "./ui_amcset_main_window.h"
#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/process.hpp>
#include <filesystem>
#include <string>
#include <cstdlib>
#include <iostream>

amcset_main_window::amcset_main_window(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::amcset_main_window)
{
    ui->setupUi(this);
    try {
        std::filesystem::path currentPath = std::filesystem::current_path();
        std::cout << "Directory: " << currentPath.string() << std::endl;
        std::string executable = "../server/amcset_server";
        std::vector<std::string> args = {};
        boost::process::spawn(executable, args);
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(1);
    }
}

amcset_main_window::~amcset_main_window()
{
    delete ui;
}

void amcset_main_window::on_pushButton_clicked()
{
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
        std::exit(0);
    }

    QString displayMessage(recieved_message);
    ui->textEdit->setPlainText(displayMessage);
    QTextCursor cursor = ui->textEdit->textCursor();
    cursor.movePosition(QTextCursor::Start, QTextCursor::MoveAnchor, 1);
}

void amcset_main_window::on_pushButton_2_clicked()
{
    boost::interprocess::message_queue client_to_server_queue(boost::interprocess::open_only, "client_to_server_message_queue");

    std::string exit_message("exit");
    client_to_server_queue.send(exit_message.c_str(), exit_message.size(), 0);
    std::exit(0);
}

