import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the CSV file
df = pd.read_csv('./output/coordinateOutput.csv')

step = 1
startidx = 0  # Python uses 0-based indexing
endidx = len(df)
data = df.iloc[startidx:endidx:step]

# Extract relevant columns
bombardmentID = data['bombardmentID']
x = data['x']
y = data['y']
z = data['z']
# Remote any NaN values
x = x.dropna()
y = y.dropna()
z = z.dropna()

# Group data by bombardmentID
uniqueIDs = bombardmentID.unique()
groupedData = {uid: data[bombardmentID == uid][['x', 'y', 'z']].values for uid in uniqueIDs}

# Plot data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for uid, group in groupedData.items():
    ax.scatter(group[:, 0], group[:, 1], group[:, 2], marker='o', s=0.1, label=f'Bombardment ID {uid}')

# ax.set_xlim([-6E8, 6E8])
# ax.set_ylim([-6E8, 6E8])
# ax.set_zlim([0, 6E8])
ax.set_xlabel('X (angstrom)')
ax.set_ylabel('Y (angstrom)')
ax.set_zlabel('Z (angstrom)')
ax.set_title('Bombardment Data for 1000 particles (by eSRIM)')
ax.view_init(6, 45)  # Set the viewing angle

plt.tight_layout() 

# ax.legend()
plt.show()