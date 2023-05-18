import os
import pdb
import json
import tetra3
import numpy as np
from pathlib import Path

lines = json.load(
    open(
        "../spectranet_hyperspectral/.scratch/shs_oft/annotations/lines_train.json",
        "r",
    )
)

# convert lines to centroids
for idx in range(len(lines["images"])):
    lim_mag = 25

    # star list, ordered by decreasing flux
    stars = [
        [annot["line"][1], annot["line"][0], annot["mag"]]
        for annot in lines["annotations"]
        if annot["image_id"] == idx
        and annot["mag"] <= lim_mag
        and annot["line"][0] > 0
        and annot["line"][1] > 0
        and annot["line"][0] < 512
        and annot["line"][1] < 512
    ]

    image = [img for img in lines["images"] if img["id"] == idx][0]
    ra = image["coords"]["CRVAL1"]
    dec = image["coords"]["CRVAL2"]

    ra_boresight = image["coords"]["CRVAL1"]  # + np.random.uniform(-0.5, 0.5)
    dec_boresight = image["coords"]["CRVAL2"]  # + np.random.uniform(-0.5, 0.5)

    grid_rad = 3
    ra_grid = round(ra_boresight, 0)
    dec_grid = round(dec_boresight, 0)

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

    db_name = "tinyFOV_14_%i_%i_%i" % (int(ra_grid), int(dec_grid), int(grid_rad))

    db_path = (Path(__file__).parent / db_name).with_suffix(".npz")
    if not os.path.exists(db_path):
        print("generating database " + db_name)
        t3.generate_database(
            0.2,
            save_as=db_name,
            star_catalog="b1",
            pattern_stars_per_fov=50,
            verification_stars_per_fov=10,
            star_max_magnitude=lim_mag,
            # star_min_separation=0.005,
            # pattern_max_error=0.01,
            temporal_corr=False,
            center_radec=[ra_grid, dec_grid],
            radec_radius_degrees=grid_rad,
        )
    else:
        print(db_name + " exists")
    t3.load_database(db_name)

    # solution = t3.solve_from_centroids(stars, 512, 512)
    # solution = t3.solve_from_centroids(stars, 512, 512)
    print(ra, dec)
    solution = t3.solve_from_centroids(
        stars, 512, 512, pattern_checking_stars=10, match_threshold=1e-1
    )

    print(json.dumps(solution, indent=1))
    print(ra, dec)
    import pdb

    pdb.set_trace()
