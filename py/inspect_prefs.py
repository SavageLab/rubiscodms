# Native Modules
import io
# Installed Modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def get_label_rotation(angle, offset):
    # Rotation must be specified in degrees :(
    rotation = np.rad2deg(angle + offset)
    if angle <= np.pi:
        alignment = "right"
        rotation = rotation + 180
    else:
        alignment = "left"
    return rotation, alignment


def add_labels(angles, values, labels, offset, ax):
	# This is the space between the end of the bar and the label
	padding = 4

	# Iterate over angles, values, and labels, to add all of them.
	for angle, value, label, in zip(angles, values, labels):
		angle = angle

		# Obtain text rotation and alignment
		rotation, alignment = get_label_rotation(angle, offset)

		# And finally add the text
		ax.text(
			x=angle,
			y=value + padding,
			s=label,
			ha=alignment,
			va="center",
			rotation=rotation,
			rotation_mode="anchor",
			fontsize=20
		)


def plot_grouped_radar(df, title, values_series, label_series, group_series):
	fontsize = 70
	# Grab the group values
	df['interp'] = np.interp(values_series, (values_series.min(), values_series.max()), (0, 100))
	GROUP = group_series.values
	VALUES = df['interp'] .values
	LABELS = label_series.values
	# Determines where to place the first bar.
	# By default, matplotlib starts at 0 (the first bar is horizontal)
	# but here we say we want to start at pi/2 (90 deg)
	OFFSET = np.pi / 2
	# Add three empty bars to the end of each group
	PAD = 3
	ANGLES_N = len(VALUES) + PAD * len(np.unique(GROUP))
	ANGLES = np.linspace(0, 2 * np.pi, num=ANGLES_N, endpoint=False)
	WIDTH = (2 * np.pi) / len(ANGLES)

	# Obtain size of each group
	GROUPS_SIZE = [len(i[1]) for i in df.groupby("group")]

	# Obtaining the right indexes is now a little more complicated
	offset = 0
	IDXS = []
	for size in GROUPS_SIZE:
		IDXS += list(range(offset + PAD, offset + size + PAD))
		offset += size + PAD

	# Same layout as above
	fig, ax = plt.subplots(figsize=(15, 15), subplot_kw={"projection": "polar"})
	ax.set_title(title, fontsize=fontsize)
	ax.tick_params(axis='both', labelsize=fontsize)  # Set the tick label font size
	ax.set_theta_offset(OFFSET)
	ax.set_ylim(-40, 100)
	ax.set_frame_on(False)
	ax.xaxis.grid(False)
	ax.yaxis.grid(False)
	ax.set_xticks([])
	ax.set_yticks([])

	# Use different colors for each group!
	GROUPS_SIZE = [len(i[1]) for i in df.groupby("group")]
	COLORS = [f"C{i}" for i, size in enumerate(GROUPS_SIZE) for _ in range(size)]

	# And finally add the bars.
	# Note again the `ANGLES[IDXS]` to drop some angles that leave the space between bars.
	ax.bar(
		ANGLES[IDXS], VALUES, width=WIDTH, color=COLORS,
		edgecolor="white", linewidth=2
	)
	add_labels(ANGLES[IDXS], VALUES, LABELS, OFFSET, ax)
	return fig, ax
	# plt.show()


def main():
	# Snakemake I/O
	# Inputs
	dms_prefs = str(snakemake.input.dms_out)
	# Outputs
	radial_prefs = str(snakemake.output.radial_prefs)
	# Params
	site_offset = int(snakemake.params.site_offset)
	# DEBUG INPUT
	# dms_prefs = "/Users/bellieny/projects/team_resources/toolbox/phylo_dms_workflow/scratch/edited_aa_preferences.csv"

	aa_groups = {'G': 'np',
	             'A': 'np',
	             'C': 'np',
	             'D': 'mch',
	             'E': 'mch',
	             'F': 'np',
	             'H': 'pch',
	             'I': 'np',
	             'K': 'pch',
	             'L': 'np',
	             'M': 'np',
	             'N': 'p',
	             'P': 'np',
	             'Q': 'p',
	             'R': 'pch',
	             'S': 'p',
	             'T': 'p',
	             'V': 'np',
	             'W': 'np',
	             'Y': 'p'
	             }

	df_aapref = pd.read_csv(dms_prefs)
	df_aapref = df_aapref.head(18)
	trns_df = df_aapref.T
	# Create a list to store each subplot
	subplots = []
	for col in trns_df:
		graph_title = col + site_offset
		trns_df['group'] = trns_df.index.map(aa_groups)
		b = trns_df[[col, 'group']].dropna(axis=0)
		b = b.rename(columns={col: 'value'})

		# Create a subplot
		# fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(8, 8))

		# Plot on the subplot
		(fig, ax) = plot_grouped_radar(b, graph_title, b['value'], b.index, b['group'])

		# Store the subplot in the list
		subplots.append((fig, ax))
		plt.close('all')
	# Create a grid of subplots
	cols = 5  # Change this based on how many columns you want
	rows = int(np.round(len(df_aapref) / cols))  # Change this based on how many rows you want

	fig, axes = plt.subplots(rows, cols, figsize=(15, 15))

	# Iterate through the subplots and place them in the grid
	for i in range(rows):
		for j in range(cols):
			idx = i * cols + j
			if idx < len(subplots):
				print(idx)
				fig, ax = subplots[idx]

				# Save the figure to a BytesIO object
				img_io = io.BytesIO()
				fig.savefig(img_io, format='png', dpi=300)
				img_io.seek(0)

				# Read the image from BytesIO and display it on the subplot
				img = plt.imread(img_io, format='png')
				axes[i, j].imshow(img)

				axes[i, j].axis("off")
			else:
				# If no subplot, remove the empty axes
				axes[i, j].axis("off")

	plt.tight_layout()
	plt.show()
	plt.savefig('delete.png', format='png')


if __name__ == "__main__":
	main()
