import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import pdb

res = json.load(open("comp_result.json", "r"))

result = {
    "frame_count": 0,
    "failure": {},
    "success": {},
    "astminuszg": [],
    "mean_err": {},
}

for item in res:
    # print(json.dumps(res[item], indent=1))

    result["frame_count"] += 1
    for meth in res[item]:
        if meth not in result["success"]:
            result["success"][meth] = 0
        if meth not in result["failure"]:
            result["failure"][meth] = 0
        if meth not in result["mean_err"]:
            result["mean_err"][meth] = []

        if res[item][meth]["resp"] == 200:
            if res[item][meth]["err"] < 0.01:
                result["success"][meth] += 1
            else:
                result["failure"][meth] += 1
            result["mean_err"][meth].append(res[item][meth]["err"])

        elif meth == "zg":
            print(item)

    if 204 not in [res[item]["a"]["resp"], res[item]["zg"]["resp"]]:
        diff = res[item]["a"]["err"] - res[item]["zg"]["err"]
        if np.abs(res[item]["a"]["err"]) < 0.01 and np.abs(
            res[item]["zg"]["err"] < 0.01
        ):
            result["astminuszg"].append(diff)
        else:
            print(item)

plt.hist(
    result["astminuszg"],
    bins=30,
)
plt.xlabel("Astrometry.net - MultiQuad center error [degrees]")


plt.xlim([-0.0008, 0.0008])
plt.axvline(0, color="black")
plt.savefig("comp_err.png", dpi=75)
# print(json.dumps(result, indent=1))
pdb.set_trace()
