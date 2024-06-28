// Copyright 2024, Texas A&M University
//
// This file is part of AMCSET.
//
// AMCSET is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with AMCSET.
// If not, see <https://www.gnu.org/licenses/>.

#include "amcset_server.h"
#include "amcset_common.h"

#include <iostream>
#include <string>
#include <boost/math/constants/constants.hpp>
#include <boost/interprocess/ipc/message_queue.hpp>

int main(const int argc, const char* argv[]) {
    std::cout << amcset::server::server_greeting() << std::endl;
    std::cout << amcset::common::common_greeting() << std::endl;
    std::cout << "Value of pi: " << boost::math::constants::pi<double>()
              << std::endl;

    boost::interprocess::message_queue::remove("server_to_client_message_queue");
    boost::interprocess::message_queue::remove("client_to_server_message_queue");

    boost::interprocess::message_queue server_to_client_message_queue(boost::interprocess::create_only, "server_to_client_message_queue", 100, 256);
    boost::interprocess::message_queue client_to_server_message_queue(boost::interprocess::create_only, "client_to_server_message_queue", 100, 256);

    size_t count = 0;

    while(true) {
        char recieved_message[256];
        size_t recieved_size;
        unsigned int priority;
        client_to_server_message_queue.receive(&recieved_message, sizeof(recieved_message), recieved_size, priority);
        recieved_message[recieved_size] = '\0';

        if(std::string(recieved_message) == "exit") {
            std::cout << "Exiting..." << std::endl;
            return 0;
        }

        std::cout << "Recieved: " << recieved_message << std::endl;

        count++;

        std::string message = "Server call count: " + std::to_string(count);

        server_to_client_message_queue.send(message.c_str(), message.size(), 0);
    }

    boost::interprocess::message_queue::remove("server_to_client_message_queue");
    boost::interprocess::message_queue::remove("client_to_server_message_queue");

    return 0;
}
