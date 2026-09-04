import matplotlib.pyplot as plt
import pandas as pd


def plot_collision_data(filename):
    # Read the CSV file
    df = pd.read_csv(filename)

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(df['Depth'], df['Collisions'], 'b.-', linewidth=2, markersize=8)

    # Customize the plot
    plt.xlabel('Depth (m)')
    plt.ylabel('Collisions')
    plt.title('Collisions vs Depth')
    plt.grid(True, linestyle='--', alpha=0.7)

    # Format x-axis to use scientific notation
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Show the plot
    plt.show()


# Usage example
if __name__ == "__main__":
    plot_collision_data("collisions.csv")
