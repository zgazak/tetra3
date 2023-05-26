import pdb
import json
from glob import glob
from astropy.io import fits

all_files = glob("/data/zgazak/astrometry/star_annots/*/ImageFiles/*.json")

for gt_annot in all_files:
    annot = json.load(open(gt_annot, "r"))
    print(gt_annot)
    pdb.set_trace()
    file = fits.open()
