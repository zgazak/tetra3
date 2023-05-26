import os
import pdb
import json
from glob import glob
from astropy.io import fits
from shs.eval.image import zscale
import matplotlib.pyplot as plt

all_files = glob("/data/zgazak/astrometry/star_annots/*/ImageFiles/*.json")


def plates_stars(annot):
    stars = []
    for obj in annot["objects"]:
        if obj["type"] == "line":
            stars.append([obj["y_center"], obj["x_center"], obj["isoflux"]])
    stars.sort(key=lambda x: x[-1])
    return stars


def prep_axes(width, height, dpi=150):
    fig = plt.figure(
        frameon=False, figsize=(1 * width / dpi, 1 * height / dpi), dpi=dpi
    )
    ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    fig.add_axes(ax)


def plot_annot(stars, ax):
    ax.plot([s[1] for s in stars], [s[0] for s in stars], marker="x", c="tab:red", lw=0)


for gt_annot in all_files:
    raw = os.path.join(
        "/data/avdberg/SatNet.v.1.4.0.0/raw_data",
        os.path.sep.join(gt_annot.split(os.path.sep)[-3:]),
    ).replace(".json", ".fits")
    if os.path.exists(raw):
        annot = json.load(open(gt_annot, "r"))
        gt_plates = plates_stars(annot)

        raw_data = fits.open(raw)[0]
        arr = zscale(raw_data.data)

        ax = prep_axes(arr.shape[0], arr.shape[1])
        ax.imshow(arr, cmap="Greys_r")
        plot_annot(gt_plates, ax)

        plt.show()
