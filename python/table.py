# Example tkinter app with a Hello, World! button
import tkinter as tk


class Application:
    def __init__(self, master=None):
        self.master = master
        self.create_widgets()

    def create_widgets(self):
        self.hi_there = tk.Button(self.master)  # Change this line
        self.hi_there["text"] = "Hello, World!\n(click me)"
        self.hi_there["command"] = self.say_hi
        self.hi_there.pack(side="top")

        self.quit = tk.Button(self.master, text="QUIT", fg="red",  # And this line
                              command=self.master.destroy)
        self.quit.pack(side="bottom")

    def say_hi(self):
        print("hi there, everyone!")


def main():
    root = tk.Tk()
    app = Application(master=root)
    root.mainloop()

if __name__ == "__main__":
    main()