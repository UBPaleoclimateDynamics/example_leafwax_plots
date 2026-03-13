"""
Example for creating stacked bar plots of plant leaf wax carbon chain-length
relative/fractional abundances.

Author: Kurt R. Lindberg
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import leafwaxtools as leafwax


# Use these parameters for easier editing in Inkscape/Illustrator
plt.rcParams['svg.fonttype'] = 'none'    # make text recognized as editable text
# plt.rcParams['font.size'] = 12    # changes all figure text font sizes
plt.rcParams['font.family'] = "Liberation Sans"    # changes all figure text fonts

# Import plant leaf wax data .csv file
wax_df = pd.read_csv('input_data/qpt_16_3a_gorbey2022.csv')

# Extract columns of only even-chain n-alkanoic acids
wax_even = wax_df[
    [
        'c22_fconc',
        'c24_fconc',
        'c26_fconc',
        'c28_fconc',
    ]
]

# Convert filtered Pandas DataFrame to Numpy array to work with leafwaxtools
wax_even_arr = np.array(wax_even)

# Calculate chain-length relative abundances from original concentration data
# calculate_percent=True (chain-length percent abundances out of 100)
# calculate_percent=False (chain-length fractional abundances out of 1)
wax_even_relabd = leafwax.Chain(wax_even_arr).relative_abd(calculate_percent=True)

# Create a new a array where the each chain-length column is the sum of the
# current chain-length plus all previous ones (e.g., new C24 = C24 + C22)
wax_sum = np.zeros(shape=np.shape(wax_even_relabd))
wax_sum[:,0] = wax_even_relabd[:,0]
for row in range(0, len(wax_even_relabd[:,0])):
    for col in range(1, len(wax_even_relabd[1,:])):
        wax_sum[row,col] = np.sum(wax_even_relabd[row,0:col+1])

# List of hexadecimal color values, make sure they're colorblind-friendly :)
colors = [
    '#1b7837',
    '#7fbf7b',
    '#af8dc3',
    '#762a83'
]


# Figure script

# I create two panel rows then delete the bottom one to force the figure 
# to be more wider than tall, kind of a hacky fix
fig, axs = plt.subplots(2,1,layout='constrained')

ax = axs[0]
ax.set_title("Chain-length stacked bar plot example")

for c in range(0, len(colors)):
    sns.barplot(
        ax=ax,
        x=wax_df.age,    # can replace with core depths, if desired
        y=wax_sum[:,c],
        color=colors[c],    # sets barplot segment color from 'colors' list
        edgecolor='black',    # creates black border around each bar
        width=1,    # sets width of each bar
        native_scale=True,    # makes x-axis distances scale with age/depth values
        zorder=(len(colors)-c)    # prevents bars from blocking each other while overlapping
    )
ax.set_xlim([6000,-100])
ax.set_xticks(
    ticks=[6000,5500,5000,4500,4000,3500,3000,2500,2000,1500,1000,500,0],    # set tick locations on the axis
    labels=[6,"",5,"",4,"",3,"",2,"",1,"",0]   # set individual tick labels
)
ax.set_xlabel("Age (cal kyr BP)")
ax.set_ylim([0,100])
ax.set_yticks(
    ticks=[0,25,50,75,100],
    labels=[0,25,50,75,100]
)
ax.set_ylabel("% Abundance")
ax.grid(visible=False)

# Create rectangular color patches for the legend
c28 = mpatches.Patch(color=colors[3], label="C28")
c26 = mpatches.Patch(color=colors[2], label="C26")
c24 = mpatches.Patch(color=colors[1], label="C24")
c22 = mpatches.Patch(color=colors[0], label="C22")

ax.legend(
    handles=[c28,c26,c24,c22],    # add the color patches we just created
    loc='center left',    # set position of the legend
    bbox_to_anchor=(1,0.5)    # more legend positioning
)

fig.delaxes(axs[1])

figure_chain_stack = plt.gcf()
figure_chain_stack.savefig('figures/chain_stack_example.png', dpi=300)
# Can change saved figure to .png, .jpeg or other image file formats
# Recommended to save figure with at least 300 dpi
