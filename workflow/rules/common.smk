def get_dwi_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="dwi",
                    den="32k",
                    atlas="{atlas}",
                    suffix="struc.pconn.png",
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
                    den="32k",
                    task="{task}",
                    denoise="{denoise}",
                    fwhm="{fwhm}",
                    atlas="{atlas}",
                    suffix="bold.pconn.png",
                    **config["subj_wildcards"],
                ),
                subject=config["subjects"][dataset],
                dataset=dataset,
                task=config["func"]["task"],
                denoise=config["func"]["denoise"].keys(),
                fwhm=config["func"]["fwhm"],
                atlas=config["atlas"].keys(),
            )
        )
    return targets



def get_sfc_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="func",
                    desc="preproc",
                    den="32k",
                    task="{task}",
                    denoise="{denoise}",
                    fwhm="{fwhm}",
                    atlas="{atlas}",
                    suffix="sfc.pscalar.nii",
                    **config["subj_wildcards"],
                ),
                subject=config["subjects"][dataset],
                dataset=dataset,
                task=config["func"]["task"],
                denoise=config["func"]["denoise"].keys(),
                fwhm=config["func"]["fwhm"],
                atlas=config["atlas"].keys(),
            )
        )

    return targets
