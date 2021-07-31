import os

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib._color_data as mcd
import matplotlib.colors

from scipy import interpolate

from math import sqrt
from math import pi

rootdir = 'C:/Users/t.kuipers/OneDrive - Ultimaker B.V/20210702-TiKui Test-Main Folder/'

sample_files = {'A': '20210701-S5-TplaPP-XY FLat-A1-5/20210701-S5-TplaPP-XY Flat-A1-5.is_tens_Exports/'
    , 'J': '20210701-S5-TplaPP-XY FLat-J1-5/20210701-S5-TplaPP-XY Flat-J1-5.is_tens_Exports/'
    , 'V': '20210716-S5-TplaPP-XY Flat-V1-5/20210716-S5-TplaPP-XY Flat-V1-5.is_tens_Exports/'
    , 'jigsaw': '20210707-S5-TplaPP-XY Flat-Jigsaw/20210701-S5-TplaPP-XY Flat-JigSaw.is_tens_Exports/'
    , 'suture': '20210707-S5-TplaPP-XY Flat-Suture WIDE/20210701-S5-TplaPP-XY Flat-Suture Wide.is_tens_Exports/'}
#    , 'suture': '20210707-S5-TplaPP-XY Flat-Suture/20210701-S5-TplaPP-XY Flat-Suture.is_tens_Exports/'}

sample_name = {'A': 'whole', 'J': 'broken wb+', 'V': 'diagonal 1.2', 'jigsaw': 'jigsaw 1.8', 'suture': 'suture 2.7'}
sample_color = {'A': 'red', 'J': 'green', 'V': 'blue', 'jigsaw': 'magenta', 'suture': 'orange'}

ideal_area = {'A': 75, 'J': 40.25, 'V': 48.644, 'jigsaw': 92.3995, 'suture': 90}
fig = plt.figure()
ax = plt.axes()

legend_elements = []

for i, (sample, sample_dirname) in enumerate(sample_files.items()):
    print('handling sample ' + sample + '...')
    color = sample_color[sample]
    color = plt.get_cmap('tab10')(i / 10)
    legend_elements.append(Patch(facecolor=color, edgecolor='none', label=sample_name[sample]))
    for subdir, dirs, files in os.walk(rootdir + sample_dirname):
        for file in files:
            filename = os.path.join(subdir, file)

            data_file = np.genfromtxt(filename, delimiter=';', skip_header=2, dtype=str)

            data_file = np.char.replace(data_file, ',', '.')
            data = data_file.astype(np.float64)


            displacement = data[:, 1]
            force = data[:, 2]
            displacement = displacement - displacement[np.argmax(force > 0)]
            stress = force * 1000 / ideal_area[sample]
            ax.plot(displacement, stress, color=color, alpha = .75)


ax.legend(handles=legend_elements)

ax.set_xlabel("displacement (mm)")
ax.set_ylabel("cell stress (MPa)")

ax.set_xlim(0, 10)
plt.show()
