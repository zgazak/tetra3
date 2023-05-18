import os
import pdb
import json
import tetra3
import numpy as np
from shs.eval.image import zscale
from pathlib import Path

base_path = "../spectranet_hyperspectral/.scratch/StarNet-67k/"
lines = json.load(
    open(
        "../spectranet_hyperspectral/.scratch/StarNet-67k/annotations/lines_test.json",
        "r",
    )
)
result = {}
# convert lines to centroids
for idx in range(len(lines["images"])):
    lim_mag = 18
    image = [img for img in lines["images"] if img["id"] == idx][0]
    height = image["height"]
    width = image["width"]
    print(json.dumps(image, indent=2))
    import pdb

    pdb.set_trace()
    # star list, ordered by decreasing flux
    stars = [
        [annot["line"][1], annot["line"][0], annot["mag"]]
        for annot in lines["annotations"]
        if annot["image_id"] == idx
        and annot["mag"] <= lim_mag
        and annot["line"][0] > 0
        and annot["line"][1] > 0
        and annot["line"][0] < image["height"]
        and annot["line"][1] < image["width"]
    ]

    ra = image["coords"]["CRVAL1"]
    dec = image["coords"]["CRVAL2"]

    ra_boresight = image["coords"]["CRVAL1"]  # + np.random.uniform(-0.5, 0.5)
    dec_boresight = image["coords"]["CRVAL2"]  # + np.random.uniform(-0.5, 0.5)
    print(json.dumps(image["coords"], indent=2))

    grid_rad = 5
    fov = 0.2
    ra_grid = int(round(2 * ra_boresight) / 2)
    dec_grid = int(round(2 * dec_boresight) / 2)

    stars.sort(key=lambda x: x[-1])
    print(stars[:4])
    # flip y axis??
    # stars = [[s[0], 256 + (-1 * (s[1] - 256))] for s in stars]

    # are x and y switched?
    # stars = [[s[1], s[0]] for s in stars]

    stars = np.array([s[:2] for s in stars])

    print(len(stars))

    print(stars[:4])

    t3 = tetra3.Tetra3()
    pattern_size = 4
    db_name = "tinyFOV_%.1f_%i_%i_%i_%i" % (
        fov,
        int(ra_grid),
        int(dec_grid),
        int(grid_rad),
        pattern_size,
    )

    db_path = Path(__file__).parent / (db_name + ".npz")

    if not os.path.exists(db_path):
        print("generating database " + db_name)
        try:
            t3.generate_database(
                fov,
                save_as=db_path,
                star_catalog="sstrc7",
                pattern_stars_per_fov=10,
                verification_stars_per_fov=20,
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
    solution = t3.solve_from_centroids(
        stars,
        512,
        512,
        pattern_checking_stars=pcs,
        match_threshold=1e-2,
        # fov_estimate=0.168,
        # fov_max_error=0.005,
    )
    count += 1
    print(count, json.dumps(solution, indent=1))
    print(ra, dec)
    if solution["RA"] is not None:
        offset_error_deg = (180 / np.pi) * np.arctan2(
            np.sin(np.deg2rad(ra) - np.deg2rad(solution["RA"])),
            np.cos(np.deg2rad(dec) - np.deg2rad(solution["Dec"])),
        )
        result[idx] = offset_error_deg
    else:
        result[idx] = None

    print(json.dumps(result, indent=1))
