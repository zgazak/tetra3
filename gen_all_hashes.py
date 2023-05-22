import os
import pdb
import json
import tetra3
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import sstr7 as cat
import argparse
import os
import subprocess
from multiprocessing import Process, Pool, Manager
import signal
from glob import glob


try:
    import daemon

except Exception:
    daemon = None


def init_worker():
    """
    Pool worker initializer for keyboard interrupt on Windows
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def writer(config_file_path, q):
    q.put(config_file_path)


def reader(q):
    config = q.get()
    print("processing " + json.dumps(config, indent=1))
    if not os.path.exists(config["save_as"]):
        t3 = tetra3.Tetra3()
        try:
            t3.generate_database(**config)
            print("database created, releasing thread")
        except:
            print("database failed")
    else:
        print("Database exists, passing")


def process_configs(configs, n_processes):
    # Create manager
    m = Manager()

    # Create multiprocessing queue
    q = m.Queue()

    # Create a group of parallel writers and start them
    for config in configs:
        Process(
            target=writer,
            args=(
                config,
                q,
            ),
        ).start()

    # Create multiprocessing pool
    p = Pool(n_processes, init_worker)

    # Create a group of parallel readers and start them
    # Number of readers is matching the number of writers
    # However, the number of simultaneously running
    # readers is constrained to the pool size

    readers = []
    for i in range(len(configs)):
        readers.append(p.apply_async(reader, (q,)))
    # Wait for the asynchrounous reader threads to finish

    [r.get() for r in readers]


if __name__ == "__main__":
    # 10 deg grid
    # for ra in (np.range(0,360))
    ras = sorted(set([round(r, -1) + 5 for r in range(355)]))
    decs = sorted(set([round(r, -1) + 5 for r in np.arange(-90, 85)]))
    grid_rad = 6
    pattern_size = 4
    num_per_fov = 6

    lim_mag = 19
    fov = 0.2

    configs = []
    for ra in ras:
        for dec in decs:
            db_name = "tFOV_fullgrid_%.1f_%i_%i_%i_%i_%i" % (
                fov,
                int(ra),
                int(dec),
                int(grid_rad),
                pattern_size,
                num_per_fov,
            )

            db_path = Path(__file__).parent / (db_name + ".npz")

            configs.append(
                {
                    "max_fov": fov,
                    "save_as": str(db_path),
                    "star_catalog": "sstrc7",
                    "catalog_location": "/data/shared/sstrc7",
                    "pattern_stars_per_fov": int(num_per_fov),
                    "verification_stars_per_fov": int(400),
                    "star_max_magnitude": lim_mag,
                    "star_min_separation": 0.001,
                    "pattern_max_error": 0.01,
                    "pattern_size": int(pattern_size),
                    "temporal_corr": True,
                    "center_radec": [int(ra), int(dec)],
                    "radec_radius_degrees": grid_rad,
                }
            )
    process_configs(configs, 50)
