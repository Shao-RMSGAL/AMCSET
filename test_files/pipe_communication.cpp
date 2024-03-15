#include <iostream>
#include <unistd.h>

int main() {
    int pipefd[2];
    if (pipe(pipefd) == -1) {
        std::cerr << "Pipe creation failed." << std::endl;
        return 1;
    }

    pid_t pid = fork();
    if (pid == -1) {
        std::cerr << "Fork failed." << std::endl;
        return 1;
    } else if (pid == 0) {  // Child process (Python script)
        close(pipefd[1]);    // Close the write end of the pipe

        dup2(pipefd[0], STDIN_FILENO);  // Redirect stdin to the read end of the pipe
        close(pipefd[0]);               // Close the read end of the pipe

        execlp("python3", "python3", "pipe_communication.py", nullptr);  // Execute the Python script
    } else {  // Parent process (C++ program)
            close(pipefd[0]);  // Close the read end of the pipe

            std::string message = "Hello from C++!";
            write(pipefd[1], message.c_str(), message.size());  // Write message to the pipe

            close(pipefd[1]);  // Close the write end of the pipe
    return 0;
    }
}
