import os
import pdb
import json
from glob import glob
from astropy.io import fits
from shs.eval.image import zscale
import matplotlib.pyplot as plt
import tetra3
from astropy import wcs
from io import BytesIO
import sstr7 as cat
import astropy.units as u
import base64
import seaborn as sns
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.coordinates import angular_separation

all_files = glob("/data/zgazak/astrometry/star_annots/*/ImageFiles/*.json")


def plates_catalog(
    fov=0.2, ra=0, dec=0, grid_coarse_deg=1, grid_rad=1, pattern_size=4, num_per_fov=6
):
    ra = round(ra, -1)
    dec = round(dec, -1)

    best_dist = None
    for r in sorted(set([round(r, -1) + grid_coarse_deg for r in range(355)])):
        for d in sorted(
            set([round(r, -1) + grid_coarse_deg for r in np.arange(-90, 85)])
        ):
            dist = angular_separation(
                np.deg2rad(ra), np.deg2rad(dec), np.deg2rad(r), np.deg2rad(d)
            )
            if best_dist is None or dist < best_dist:
                best_dist = dist
                cra = r
                cdec = d
                print(np.rad2deg(best_dist), cra, cdec)

    db_name = "tFOV_fullgrid_%.1f_%i_%i_%i_%i_%i.npz" % (
        fov,
        int(cra),
        int(cdec),
        int(grid_rad),
        pattern_size,
        num_per_fov,
    )

    db_gen_params = {
        "max_fov": fov * 1.1,
        "save_as": str(db_name),
        "star_catalog": "sstrc7",
        "catalog_location": "/data/shared/sstrc7",
        "pattern_stars_per_fov": int(num_per_fov),
        "verification_stars_per_fov": int(400),
        "star_max_magnitude": 17,
        "star_min_separation": 0.001,
        "pattern_max_error": 0.01,
        "pattern_size": int(pattern_size),
        "temporal_corr": True,
        "center_radec": [int(cra), int(cdec)],
        "radec_radius_degrees": grid_rad,
    }
    return db_name, db_gen_params


def plates_stars(annot, width, height):
    stars = []
    for obj in annot["objects"]:
        if obj["type"] == "line":
            stars.append(
                [obj["y_center"] * height, obj["x_center"] * width, obj["iso_flux"]]
            )
    stars.sort(key=lambda x: x[-1])
    return stars


def dotnet_stars(annot, width, height):
    anet_list = []
    for obj in annot["objects"]:
        if obj["type"] == "line":
            anet_list.append(
                {
                    "xpix": obj["x_center"] * width,
                    "ypix": obj["y_center"] * height,
                    "flux": obj["iso_flux"],
                }
            )
    anet_list.sort(key=lambda x: x["flux"])
    pdb.set_trace()
    return anet_list


def prep_axes(width, height, dpi=150):
    fig = plt.figure(
        frameon=False, figsize=(1 * width / dpi, 1 * height / dpi), dpi=dpi
    )
    ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    fig.add_axes(ax)
    return ax


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


def plot_annot(stars, ax):
    ax.plot([s[1] for s in stars], [s[0] for s in stars], marker="x", c="tab:red", lw=0)


def plot_cat(solution, ax):
    # print(json.dumps(solution, indent=2))
    wc = sol2wcs(solution)
    wf = wcs.WCS(wc)

    starsf = cat.query_by_los_radec(
        solution["FOV"] * 2,
        solution["FOV"] * 2,
        solution["RA"],
        solution["Dec"],
        rootPath="/data/shared/sstrc7",
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
            star["mv"] < 17
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


for gt_annot in all_files:
    raw = os.path.join(
        "/data/avdberg/SatNet.v.1.4.0.0/raw_data",
        os.path.sep.join(gt_annot.split(os.path.sep)[-3:]),
    ).replace(".json", ".fits")
    if os.path.exists(raw):
        t3 = tetra3.Tetra3()

        annot = json.load(open(gt_annot, "r"))
        raw_data = fits.open(raw)[0]
        arr = zscale(raw_data.data)
        width = arr.shape[0]
        height = arr.shape[1]

        gt_plates = plates_stars(annot, width, height)
        stars = np.array([s[:2] for s in gt_plates])

        expected_ra_hms = np.array(
            [float(x) for x in raw_data.header["OBJCTRA"].split(" ")]
        )
        expected_dec_dms = [float(x) for x in raw_data.header["OBJCTDEC"].split(" ")]

        expected_ra_deg = (expected_ra_hms / np.array([1, 60, 3600])).sum()
        expected_dec_deg = (expected_dec_dms / np.array([1, 60, 3600])).sum()

        plate_db, plate_cfg = plates_catalog(
            fov=0.5, ra=expected_ra_deg, dec=expected_dec_deg, grid_rad=2
        )
        if not os.path.exists(plate_db):
            t3.generate_database(**plate_cfg)

        ax = prep_axes(width, height)
        ax.imshow(arr, cmap="Greys_r")

        plot_annot(gt_plates, ax)

        t3.load_database(plate_db)
        solution = t3.solve_from_centroids2(
            stars,
            width,
            height,
            pattern_checking_stars=10,
            match_threshold=1e-4,
            match_radius=0.001,
        )

        if solution["RA"] is not None:
            # result[fname]["zg"]["resp"] = 200
            print(solution["FOV"], solution["RA"], solution["Dec"], solution["Roll"])
            offset_error_deg = (180 / np.pi) * np.arctan2(
                np.sin(np.deg2rad(expected_ra_deg) - np.deg2rad(solution["RA"])),
                np.cos(np.deg2rad(expected_dec_deg) - np.deg2rad(solution["Dec"])),
            )
            # result[fname]["zg"]["err"] = np.abs(offset_error_deg)

            solution["width"] = width
            solution["height"] = height

            # ax = fig.add_subplot(1, 1, 1)
            plot_cat(solution, ax)

        plt.show()
        if solution["wcs"]:
            astro = fits.open(BytesIO(base64.b64decode(solution["wcs"].encode())))
            pdb.set_trace()
