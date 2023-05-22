import os
import pdb
import json
import tetra3
import numpy as np
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt
import sstr7 as cat
import astropy.units as u
import seaborn as sns
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
import arrow
from shs.eval.wcs import astrometry_server

base_path = "../spectranet_hyperspectral/.scratch/shs_oft/"

lines = json.load(
    open(
        "../spectranet_hyperspectral/.scratch/shs_oft/annotations/lines_train.json",
        "r",
    )
)
result = {}


def plot_annot(stars, ax):
    ax.plot([s[1] for s in stars], [s[0] for s in stars], marker="x", c="tab:red", lw=0)


def sol2wcs(solution):
    cd = solution["FOV"] / solution["width"]

    wc = {
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "CRVAL1": solution["RA"],
        "CRVAL2": solution["Dec"],
        "CUNIT1": "deg",
        "CUNIT2": "deg",
        "CRPIX1": solution["width"] / 2,
        "CRPIX2": solution["height"] / 2,
        # "CDELT1": np.deg2rad(solution["FOV"] / solution["width"]),
        # "CDELT2": np.deg2rad(solution["FOV"] / solution["height"]),
        "CD1_1": cd * np.cos(np.deg2rad(solution["Roll"] - 180)),
        "CD1_2": cd * -np.sin(np.deg2rad(solution["Roll"] - 180)),
        "CD2_1": cd * np.sin(np.deg2rad(solution["Roll"] - 180)),
        "CD2_2": cd * np.cos(np.deg2rad(solution["Roll"] - 180)),
        "NAXIS": 2,
        "NAXIS1": solution["width"],
        "NAXIS2": solution["height"],
    }
    return wc


def plot_cat(solution, ax):
    # print(json.dumps(solution, indent=2))
    wc = sol2wcs(solution)
    wf = wcs.WCS(wc)

    starsf = cat.query_by_los_radec(
        solution["FOV"] * 2,
        solution["FOV"] * 2,
        solution["RA"],
        solution["Dec"],
        rootPath="/home/zach.gazak/wd/sstrc7",
    )
    starxw = []
    staryw = []
    for star in starsf:
        stc = SkyCoord(
            ra=np.rad2deg(star["ra"] + (star["ra_pm"] * 3.154e7 * (2018.0 - 1991.25)))
            * u.deg,
            dec=np.rad2deg(
                star["dec"] + (star["dec_pm"] * 3.154e7 * (2018.0 - 1991.25))
            )
            * u.deg,
            frame="icrs",
            obstime="2018-01-02T12:34:56",
        )
        xf, yf = wf.world_to_pixel(stc)
        xf, yf = wf.wcs_world2pix(
            [[np.rad2deg(star["ra"]), np.rad2deg(star["dec"])]], 0
        )[0]
        if (
            star["mv"] < 19
            and xf > 0
            and xf < solution["width"]
            and yf > 0
            and yf < solution["height"]
        ):
            staryw.append(yf)
            starxw.append(xf)

    # plt.scatter(stary, starx, s=7, color="red")
    ax.plot(starxw, staryw, marker="+", c="black", lw=0)
    catx = []
    caty = []

    for star in solution["catmatch"]:
        x, y = wf.wcs_world2pix([[np.rad2deg(star[0]), np.rad2deg(star[1])]], 0)[0]
        catx.append(x)
        caty.append(y)

    ax.plot(catx, caty, marker="o", c="blue", lw=0, markerfacecolor="None")
    colors = sns.color_palette("Spectral", n_colors=len(solution["quads"]))

    for idxx, quad in enumerate(solution["quads"]):
        qua = [wf.wcs_world2pix(np.rad2deg(s[0]), np.rad2deg(s[1]), 0) for s in quad]

        for idx in range(len(qua) - 1):
            ax.plot(
                [qua[0][0], qua[idx + 1][0]],
                [qua[0][1], qua[idx + 1][1]],
                c=colors[idxx],
                lw=1,
            )


width = 512
height = 512
dpi = 150

# convert lines to centroids
for idx in range(len(lines["images"])):
    idx += 0
    lim_mag = 19
    fname = "%.9i_ast.png" % idx
    result[fname] = {}
    # star list, ordered by decreasing flux
    stars = [
        [annot["line"][1], annot["line"][0], annot["mag"]]
        for annot in lines["annotations"]
        if annot["image_id"] == idx
        and annot["mag"] <= lim_mag
        and annot["line"][0] > 0
        and annot["line"][1] > 0
        and annot["line"][0] < width
        and annot["line"][1] < height
    ]

    image = [img for img in lines["images"] if img["id"] == idx][0]

    arr = np.array(
        Image.open(os.path.join(base_path, "train_sum_zs", image["file_name"]))
    )
    ra = image["coords"]["CRVAL1"]
    dec = image["coords"]["CRVAL2"]

    ra_boresight = image["coords"]["CRVAL1"]  # + np.random.uniform(-0.5, 0.5)
    dec_boresight = image["coords"]["CRVAL2"]  # + np.random.uniform(-0.5, 0.5)
    print(json.dumps(image["coords"], indent=2))
    grid_rad = 1
    fov = 0.2
    ra_grid = int(round(ra_boresight))
    dec_grid = int(round(dec_boresight))

    stars.sort(key=lambda x: x[-1])
    print(stars[:4])
    # flip y axis??
    # stars = [[s[0], 256 + (-1 * (s[1] - 256))] for s in stars]

    # are x and y switched?
    # stars = [[s[1], s[0]] for s in stars]
    anet_list = [
        {
            "xpix": annot[1],
            "ypix": annot[0],
            "flux": 1e7 * 10 ** (annot[2] / -2.5),
        }
        for annot in stars
    ]
    stars = np.array([s[:2] for s in stars])

    print(len(stars))

    print(stars[:4])

    t3 = tetra3.Tetra3()
    pattern_size = 4
    num_per_fov = 6
    db_name = "tinyFOVallstar_%.1f_%i_%i_%i_%i_%i" % (
        fov,
        int(ra_grid),
        int(dec_grid),
        int(grid_rad),
        pattern_size,
        num_per_fov,
    )

    db_path = Path(__file__).parent / (db_name + ".npz")

    if not os.path.exists(db_path):
        print("generating database " + db_name)
        try:
            t3.generate_database(
                fov,
                save_as=db_path,
                star_catalog="sstrc7",
                pattern_stars_per_fov=num_per_fov,
                verification_stars_per_fov=400,
                star_max_magnitude=lim_mag,
                star_min_separation=0.001,
                pattern_max_error=0.01,
                pattern_size=pattern_size,
                temporal_corr=True,
                center_radec=[ra_grid, dec_grid],
                radec_radius_degrees=grid_rad,
            )
        except:
            print("database build failed")
            pdb.set_trace()
    else:
        print(db_name + " exists")
    t3.load_database(db_path)

    # solution = t3.solve_from_centroids(stars, 512, 512)
    # solution = t3.solve_from_centroids(stars, 512, 512)
    print(ra, dec)
    solution = {"RA": None}
    count = 20
    # while solution["RA"] is None and count < 21:
    pcs = count
    t0 = arrow.now()
    solution2 = t3.solve_from_centroids2(
        stars,
        512,
        512,
        pattern_checking_stars=10,
        match_threshold=1e-4,
        match_radius=0.001,
        fov_estimate=0.18,
        fov_max_error=0.01,
    )
    # pdb.set_trace()
    result[fname]["zg"] = {"tfit": (arrow.now() - t0).total_seconds()}
    solution = t3.solve_from_centroids(
        stars,
        512,
        512,
        pattern_checking_stars=10,
        match_threshold=1e-4,
        match_radius=0.001,
        fov_estimate=0.18,
        fov_max_error=0.01,
    )
    # pdb.set_trace()
    result[fname]["t"] = {"tfit": (arrow.now() - t0).total_seconds()}

    t0 = arrow.now()
    astronetsol = astrometry_server(
        "http://localhost:8080/",
        anet_list,
        width=width,
        height=height,
        fov_upper_bound=0.2,
    )
    result[fname]["a"] = {
        "tfit": (arrow.now() - t0).total_seconds(),
        "resp": astronetsol[1],
    }

    if astronetsol[1] == 200:
        offset_error_deg = (180 / np.pi) * np.arctan2(
            np.sin(np.deg2rad(ra) - np.deg2rad(astronetsol[0][0].header["CRVAL1"])),
            np.cos(np.deg2rad(dec) - np.deg2rad(astronetsol[0][0].header["CRVAL2"])),
        )
        result[fname]["a"]["err"] = np.abs(offset_error_deg)

    # print(count, json.dumps(solution, indent=1))
    print(ra, dec)

    fig = plt.figure(
        frameon=False, figsize=(1 * width / dpi, 1 * height / dpi), dpi=dpi
    )
    ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(arr, cmap="Greys_r")  # , aspect="auto")

    if solution2["RA"] is not None:
        result[fname]["zg"]["resp"] = 200
        print(solution2["RA"], solution2["Dec"], solution2["Roll"])
        offset_error_deg = (180 / np.pi) * np.arctan2(
            np.sin(np.deg2rad(ra) - np.deg2rad(solution2["RA"])),
            np.cos(np.deg2rad(dec) - np.deg2rad(solution2["Dec"])),
        )
        result[fname]["zg"]["err"] = np.abs(offset_error_deg)

        solution2["width"] = width
        solution2["height"] = height

        # ax = fig.add_subplot(1, 1, 1)
        plot_annot(stars, ax)
        plot_cat(solution2, ax)
    else:
        result[fname]["zg"]["resp"] = 204

    ax.set_axis_off()
    # ax.tight_layout()
    fig.savefig(fname, dpi=dpi)
    if solution["RA"] is not None:
        result[fname]["t"]["resp"] = 200
        print(solution["RA"], solution["Dec"], solution["Roll"])
        offset_error_deg = (180 / np.pi) * np.arctan2(
            np.sin(np.deg2rad(ra) - np.deg2rad(solution["RA"])),
            np.cos(np.deg2rad(dec) - np.deg2rad(solution["Dec"])),
        )
        result[fname]["t"]["err"] = np.abs(offset_error_deg)

        solution["width"] = width
        solution["height"] = height

    else:
        result[fname]["t"]["resp"] = 204

    print(json.dumps(result, indent=1))

    with open("comp_result.json", "w") as fp:
        fp.write(json.dumps(result, indent=4))
