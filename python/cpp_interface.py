import pandas as pd
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import os
import subprocess
import time

# Create initial empty plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sc = ax.plot([], [], [], c='r', marker='o')[0]

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Scatter Plot')
ax.set_xlim(-600, 600)
ax.set_ylim(-600, 600)
ax.set_zlim(0, 800)


plt.ion()  # Turn on interactive mode

while True:
    # Generate data
    os.chdir("./cpp")
    bash_command = "./bin/eSRIM"
    subprocess.run(bash_command, shell=True)
    os.chdir("..")

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv("./cpp/output/coordinateOutput.csv")

    # Extract x, y, z data
    x = df["x"]
    y = df["y"]
    z = df["z"]

    # Update scatter plot with new data
    sc.set_data(x, y)
    sc.set_3d_properties(z)
    plt.draw()
    plt.pause(0.01)

    # Wait for one second before the next iteration
    time.sleep(1)