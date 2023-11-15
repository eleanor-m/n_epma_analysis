"""Helper functions for use in notebooks to compile summaries"""

import pandas as pd
import numpy as np
from IPython.display import display


def get_n_and_detlim(results, results_detlim, suffix):
    """For a given suffix (i.e. correction level) compile N wt and detlims for each spot"""

    results = results.copy()
    results_detlim = results_detlim.copy()

    # Get the N wt%
    results_n = {
        k: v.loc[["N"], :].T for k, v in results["wtdata"].items() if suffix in k
    }

    results_n = {
        k: v.loc[
            [c for c in v.index if c not in ["average", "stdev", "minimum", "maximum"]],
            :,
        ]
        for k, v in results_n.items()
    }

    for k, v in results_n.items():
        v.columns.name = None

    # Get the detection limits
    results_detlim_n = {
        k: v.loc[["N"], :].T for k, v in results_detlim["wtdata"].items() if suffix in k
    }

    if len(results_detlim_n) == 0:
        # If there were no detection limits calculated via calczaf (eg for the base)
        # then just add nulls and handle this later
        for k, v in results_n.items():
            v["N detlim"] = np.NaN
    else:
        results_detlim_n = {
            k: v.loc[
                [
                    c
                    for c in v.index
                    if c not in ["average", "stdev", "minimum", "maximum"]
                ],
                :,
            ]
            for k, v in results_detlim_n.items()
        }

        for k, v in results_n.items():
            v["N detlim"] = results_detlim_n[k + "_detlim"]

    return results_n


def compile_n_summary(
    suffix_list, results, results_detlim, sampledata, datalist, summary_tables, samples
):
    stdev_string = [
        "original.kraw_stdev_pcnt",
        "montecarlo.kraw_stdev_pcnt",
        "montecarlo.kraw_stdev_apf_pcnt",
    ]

    summary_by_suffix = {}

    for i, suffix in enumerate(suffix_list):
        # Get n wt% and detlim (grouped by sample)
        results_n = get_n_and_detlim(results, results_detlim, suffix)

        # Get stdev per spot (not grouped; join with datalist to enable grouping)
        stdev_pct_per_spot = pd.concat(
            [datalist[["sample", "comment"]], summary_tables[0][[stdev_string[i]]]],
            axis=1,
        )

        stdev_pct_per_spot.rename(
            columns={stdev_string[i]: "N stdev pct (individual spots)"}, inplace=True
        )

        # Group the stdev per spot and join with n wt% and detlim
        stdev_by_sample = dict(list(stdev_pct_per_spot.groupby("sample")))

        N_by_sample = dict()

        for s in samples:
            N_by_sample[s] = pd.concat(
                [
                    stdev_by_sample[s].reset_index(drop=True),
                    results_n[f"{s}_{suffix}"].reset_index(drop=True),
                ],
                axis=1,
            )

            N_by_sample[s]["N detlim (orig)"] = (
                np.array(
                    [
                        spot.peak.loc[
                            spot.peak["element"] == "N", 
                            "dl_ppm"
                            ].values[0]  # take the first detection limit
                        for spot in sampledata[s]
                    ]
                )
                / 10000
            )

            N_by_sample[s]["N stdev abs (individual spots)"] = (
                results_n[f"{s}_{suffix}"]["N"].reset_index(drop=True)
                * stdev_by_sample[s]["N stdev pct (individual spots)"].reset_index(
                    drop=True
                )
                / 100
            )

            # Reorder the columns
            N_by_sample[s] = N_by_sample[s][
                [
                    "sample",
                    "comment",
                    "N",
                    "N stdev pct (individual spots)",
                    "N stdev abs (individual spots)",
                    "N detlim",
                    "N detlim (orig)",
                ]
            ]

            print(f"Sample: {s}, suffix: {suffix}")
            display(N_by_sample[s])

        # Summarise each sample by taking mean and stdev
        df = pd.concat(N_by_sample.values(), axis=0).reset_index()
        summary_by_suffix[suffix] = (
            df.rename(columns={"level_0": "sample"})
            .groupby("sample")[
                ["N", "N detlim", "N detlim (orig)", "N stdev abs (individual spots)"]
            ]
            .mean()
        ).rename(columns={"N": "N wt% (average)"})

        summary_by_suffix[suffix]["N stdev abs (multiple spots)"] = (
            df.rename(columns={"level_0": "sample"}).groupby("sample")[["N"]].std()
        )

        summary_by_suffix[suffix]

    compiled_summary = pd.concat(summary_by_suffix, axis=1)

    return (compiled_summary, N_by_sample)
