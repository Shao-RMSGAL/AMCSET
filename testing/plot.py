import pandas as pd
import matplotlib.pyplot as plt
#  from mpl_toolkits.mplot3d import Axes3D

if __name__ == "__main__":

    # Read the CSV file
    df = pd.read_csv('output.csv')

    # Create a 3D plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the x, y, z coordinates
    ax.scatter(df['x'], df['y'], df['z'], c=df['scalar'],
               cmap='viridis', s=0.1)
    #  ax.plot(df['x'], df['y'], df['z'], color='red', alpha=0.5)
    ax.set_xlim([-1000, 1000])
    ax.set_ylim([-1000, 1000])
    ax.set_zlim([0, 1000])

    # Add labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Visualization of Quaternion Data')

    # Add a color bar to show scalar values
    plt.colorbar(ax.scatter(df['x'], df['y'], df['z'],
                 c=df['scalar'], cmap='viridis'))

    # Show the plot
    plt.tight_layout()
    plt.show()
