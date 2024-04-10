import pandas as pd
from mayavi import mlab
import numpy as np

# User settings
step = 1  # How many data points to skip when rendering. 1 to show all data
electron = False  # Is this electron data? False if ion data
numRender = 1  # How many simulations to render?
scatter = True  # True for scatter plot, False for line plot
filename = "./output/coordinateOutput.csv"  # Filename and directory 

# Read the CSV file
df = pd.read_csv(filename, dtype={0: str})

startidx = 0  # Python uses 0-based indexing
endidx = len(df)
data = df.iloc[startidx:endidx:step]

# Extract relevant columns and remove NaN values
bombardmentID = data['bombardmentID']
x = data['x'].dropna()
y = data['y'].dropna()
z = data['z'].dropna()


# mlab.points3d(x,y,z, scale_factor=5)

# Group data by bombardmentID
uniqueIDs = bombardmentID.unique()
groupedData = {uid: data[bombardmentID == uid][['x', 'y', 'z']].values for uid
               in uniqueIDs}

# Extract data that has a depth=1 from the "data" object
# groupedData = {uid: data[bombardmentID == uid][['x', 'y', 'z']].values for uid
#                in uniqueIDs if data[bombardmentID == uid]['depth'].values[0] == 1}

if electron:
    size = 5000000
else:
    size = 5

numberOfParticles = bombardmentID.nunique()

print(f'Number of particles: {numberOfParticles}')

# Add axes with labels
mlab.axes(xlabel='X (Angstrom)', ylabel='Y (Angstrom)', zlabel='Z (Angstrom)')

if electron:
    # Add a title with numberOfParticles
    mlab.title(f'Bombardment Data for {numberOfParticles}'
               'Electrons (by eSRIM)')
else:
    # Add a title with numberOfParticles
    mlab.title(f'Bombardment Data for {numberOfParticles} Ions (by eSRIM)')

# Seed the random number generator with a constant
np.random.seed(0)

# Plot data
for uid, group in groupedData.items():
    # Generate a random color
    color = (np.random.random(), np.random.random(), np.random.random())
    if scatter:
        mlab.points3d(group[:, 0],
                      group[:, 1],
                      group[:, 2],
                      scale_factor=size,
                      color=color)
    else:
        mlab.plot3d(group[:, 0],
                    group[:, 1],
                    group[:, 2],
                    tube_radius=size,
                    color=color)
    numRender -= 1
    if numRender == 0:
        break
mlab.axes(xlabel='X (Angstrom)', ylabel='Y (Angstrom)', zlabel='Z (Angstrom)')

mlab.show()