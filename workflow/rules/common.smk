def get_dwi_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="dwi",
                    atlas="{atlas}",
                    suffix="struc.conn.csv",
                    **config["subj_wildcards"],
                ),
                subject=config["subjects"][dataset],
                dataset=dataset,
                atlas=config["atlas"].keys(),
            )
        )
    return targets


def get_func_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="func",
                    desc="preproc",
                    space="{space}",
                    den="32k",
                    task="{task}",
                    denoise="{denoise}",
                    fwhm="{fwhm}",
                    atlas="{atlas}",
                    suffix="bold.pconn.nii",
                    **config["subj_wildcards"],
                ),
                subject=config["subjects"][dataset],
                dataset=dataset,
                task=config["func"]["task"],
                denoise=config["func"]["denoise"].keys(),
                space=config["func"]["space"],
                fwhm=config["func"]["fwhm"],
                atlas=config["atlas"].keys(),
            )
        )
    return targets
