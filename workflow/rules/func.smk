

wildcard_constraints:
    subject="[a-zA-Z0-9]+",
    task="[a-zA-Z0-9]+",
    denoise="[a-zA-Z0-9]+",
    atlas="[a-zA-Z0-9]+",
    fwhm="[a-zA-Z0-9]+",


rule map_bold_to_surface_fsLR:
    input:
        bold_preproc=lambda wildcards: config["input_path"]["bold_nii"][
            wildcards.dataset
        ],
        mid_surf=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            den="32k",
            suffix="midthickness.surf.gii",
            **config["subj_wildcards"],
        ),
        white_surf=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            den="32k",
            suffix="white.surf.gii",
            **config["subj_wildcards"],
        ),
        pial_surf=bids(
            root=root,
            datatype="anat",
            hemi="{hemi}",
            den="32k",
            suffix="pial.surf.gii",
            **config["subj_wildcards"],
        ),
    output:
        metric=bids(
            root=root,
            datatype="func",
            hemi="{hemi}",
            desc="preproc",
            den="32k",
            task="{task}",
            suffix="bold.dtseries.func.gii",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -volume-to-surface-mapping {input.bold_preproc} {input.mid_surf}  {output.metric} -ribbon-constrained {input.white_surf} {input.pial_surf}"


#--- this was added new:
rule resample_structure_volume_to_mni_boldref:
    """ this is the HCP grayordinates segmentation, with segmentation labels for
        the standard 19 structures. It is used for sampling the BOLD data, then
        the actual parcellation is applied later (with a cifti dlabel in the same
        grayordinates space)."""
    input:
        dseg=lambda wildcards: config["cifti"]["structure_label_volume"],
        ref=lambda wildcards: config["input_path"]["bold_mni_ref"][
            wildcards.dataset
        ],
    output:
        dseg=bids(
            root=root,
            datatype="func",
            desc="structurelabelvolume",
            space="boldref",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -volume-resample  {input.dseg}  {input.ref} ENCLOSING_VOXEL {output.dseg}"

rule resample_parc_volume_to_mni_boldref:
    input:
        dseg=lambda wildcards: config["atlas"][wildcards.atlas]["label_volume"],
        ref=lambda wildcards: config["input_path"]["bold_mni_ref"][
            wildcards.dataset
        ],
    output:
        dseg=bids(
            root=root,
            datatype="func",
            atlas="{atlas}",
            desc="labelvolume",
            space="boldref",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -volume-resample  {input.dseg}  {input.ref} ENCLOSING_VOXEL {output.dseg}"



rule create_cifti_label_boldref:
    input:
        structure_label_volume=rules.resample_structure_volume_to_mni_boldref.output.dseg,
        label_volume=rules.resample_parc_volume_to_mni_boldref.output.dseg,
        label_left="resources/atlas/atlas-{atlas}_space-fsLR_hemi-L_parc.label.gii", 
        label_right="resources/atlas/atlas-{atlas}_space-fsLR_hemi-R_parc.label.gii",
    output:
        cifti=bids(
            root=root,
            datatype="func",
            atlas="{atlas}",
            den="91k",
            space="boldref",
            suffix="parc.dlabel.nii",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -cifti-create-label {output} -volume {input.label_volume} {input.structure_label_volume} "
        " -left-label {input.label_left} -right-label {input.label_right} "

rule create_bold_cifti:
    input:
        left_metric=bids(
            root=root,
            datatype="func",
            hemi="L",
            desc="preproc",
            den="32k",
            task="{task}",
            suffix="bold.dtseries.func.gii",
            **config["subj_wildcards"],
        ),
        right_metric=bids(
            root=root,
            datatype="func",
            hemi="R",
            desc="preproc",
            den="32k",
            task="{task}",
            suffix="bold.dtseries.func.gii",
            **config["subj_wildcards"],
        ),
        bold_preproc=lambda wildcards: config["input_path"]["bold_mni"][
            wildcards.dataset
        ],
        dseg=bids(
            root=root,
            datatype="func",
            desc="structurelabelvolume",
            space="boldref",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    output:
        cifti=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            space="boldref",
            task="{task}",
            suffix="bold.dtseries.nii",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -cifti-create-dense-timeseries {output.cifti} -left-metric {input.left_metric} -right-metric {input.right_metric} "
        " -volume {input.bold_preproc} {input.dseg}"


rule denoise_cifti:
    input:
        cifti=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            space="boldref",
            task="{task}",
            suffix="bold.dtseries.nii",
            **config["subj_wildcards"],
        ),
        json=lambda wildcards: config["input_path"]["bold_json"][wildcards.dataset],
        confounds_tsv=lambda wildcards: config["input_path"]["bold_confounds"][
            wildcards.dataset
        ],
    params:
        denoise_params=lambda wildcards: config["func"]["denoise"][wildcards.denoise],
    output:
        cifti=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            space="boldref",
            task="{task}",
            denoise="{denoise}",
            suffix="bold.dtseries.nii",
            **config["subj_wildcards"],
        ),
    group:
        "grouped_subject"
    threads: 8
    resources:
        mem_mb=32000,
    script:
        "../scripts/denoise_cifti.py"


rule smooth_cifti:
    input:
        cifti=rules.denoise_cifti.output.cifti,
        left_surf=bids(
            root=root,
            datatype="anat",
            hemi="L",
            den="32k",
            suffix="midthickness.surf.gii",
            **config["subj_wildcards"],
        ),
        right_surf=bids(
            root=root,
            datatype="anat",
            hemi="R",
            den="32k",
            suffix="midthickness.surf.gii",
            **config["subj_wildcards"],
        ),
    params:
        fwhm=lambda wildcards: float(wildcards.fwhm),
    output:
        cifti=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            space="boldref",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            suffix="bold.dtseries.nii",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -cifti-smoothing {input.cifti} {params.fwhm} {params.fwhm} COLUMN "
        " {output.cifti} -fwhm -left-surface {input.left_surf} -right-surface {input.right_surf}"


rule parcellate_bold:
    input:
        cifti_dtseries=rules.smooth_cifti.output.cifti,
        cifti_dlabel=bids(
            root=root,
            datatype="func",
            atlas="{atlas}",
            den="91k",
            space="boldref",
            suffix="parc.dlabel.nii",
            **config["subj_wildcards"],
        ),

    params:
        exclude_opt=(
            "-exclude-outliers {nstdev} {nstdev}".format(
                    nstdev=config["func"]["parcellation"]["n_stdevs_exclude"]
                )
                if config["func"]["parcellation"]["do_exclude_outliers"]
            else ""
        ),
    output:
        cifti_ptseries=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="bold.ptseries.nii",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -cifti-parcellate {input.cifti_dtseries} {input.cifti_dlabel} "
        " COLUMN {output.cifti_ptseries} {params.exclude_opt}"


rule correlate_parcels:
    input:
        cifti=rules.parcellate_bold.output.cifti_ptseries,
    params:
        fisher_z="-fisher-z" if config["func"]["correlation"]["do_fisher_z"] else "",
    output:
        cifti=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="bold.pconn.nii",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "grouped_subject"
    shell:
        "wb_command -cifti-correlation {input.cifti} {output.cifti} {params.fisher_z} "

rule struc_conn_csv_to_pconn_cifti:
    input:
        #for reference pconn, pick the first task, denoise, fwhm values from config (any would do)
        ref_cifti_pconn=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            task=config["func"]["task"][0],
            denoise=next(iter(config["func"]["denoise"])),
            fwhm=config["func"]["fwhm"][0],
            atlas="{atlas}",
            suffix="bold.pconn.nii",
            **config["subj_wildcards"],
        ),
        conn_csv=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            suffix="{struc}.conn.csv",
            **config["subj_wildcards"],
        ),
    output:
        cifti_pconn=bids(
            root=root,
            datatype="dwi",
            den="91k",
            atlas="{atlas}",
            suffix="{struc,struc|strucFA}.pconn.nii",
            **config["subj_wildcards"],
        ),
    group:
        "grouped_subject"
    script:
        "../scripts/struc_conn_csv_to_pconn_cifti.py"


rule calc_degree:
    input:
        pconn="{prefix}_{suffix}.pconn.nii",
    output:
        pscalar="{prefix}_{suffix}degree.pscalar.nii",
    group:
        "grouped_subject"
    container:
        config["singularity"]["diffparc"]
    shell:
        "wb_command -cifti-reduce {input} SUM {output}"


rule calc_sfc:
    input:
        pconn_struc=bids(
            root=root,
            datatype="dwi",
            den="91k",
            atlas="{atlas}",
            suffix="struc.pconn.nii",
            **config["subj_wildcards"],
        ),
        pconn_func=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="bold.pconn.nii",
            **config["subj_wildcards"],
        ),
    output:
        pscalar_sfc=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="sfc.pscalar.nii",
            **config["subj_wildcards"],
        ),
    group:
        "grouped_subject"
    script:
        "../scripts/calc_sfc.py"

rule calc_network_fc:
    input:
        pconn=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="bold.pconn.nii",
            **config["subj_wildcards"],
        ),
        label_tsv="resources/atlas/atlas-{atlas}_dseg.tsv",
    output:
        pconn=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="netbold.pconn.nii",
            **config["subj_wildcards"],
        ),
    group:
        "grouped_subject"
    script:
        "../scripts/calc_network_conn.py"

rule calc_network_sc:
    input:
        pconn=bids(
            root=root,
            datatype="dwi",
            den="91k",
            atlas="{atlas}",
            suffix="{struc}.pconn.nii",
            **config["subj_wildcards"],
        ),
        label_tsv="resources/atlas/atlas-{atlas}_dseg.tsv",
    output:
        pconn=bids(
            root=root,
            datatype="dwi",
            den="91k",
            atlas="{atlas}",
            suffix="net{struc}.pconn.nii",
            **config["subj_wildcards"],
        ),
    group:
        "grouped_subject"
    script:
        "../scripts/calc_network_conn.py"


rule calc_network_sfc:
    input:
        pconn_sc=bids(
            root=root,
            datatype="dwi",
            den="91k",
            atlas="{atlas}",
            suffix="struc.pconn.nii",
            **config["subj_wildcards"],
        ),
        pconn_fc=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="bold.pconn.nii",
            **config["subj_wildcards"],
        ),
        label_tsv="resources/atlas/atlas-{atlas}_dseg.tsv",
    output:
        pconn=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="91k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="netsfc.pconn.nii",
            **config["subj_wildcards"],
        ),
    group:
        "grouped_subject"
    script:
        "../scripts/calc_network_sfc.py"


