# Prepare references and deconvolution inputs for the spatial module.
spatial_pp = config["spatial_pp"]
prepare_refs_config = spatial_pp["prepare_refs"]
dot_config = spatial_pp["dot"]


rule prepare_refs:
    input:
        **prepare_refs_config["input"]
    output:
        **prepare_refs_config["output"]
    params:
        liver_samples=spatial_pp["liver_samples"],
        **prepare_refs_config["params"]
    resources:
        mem_mb = 100000,
        slurm = "gres=gpu:1"
    conda:
        "../envs/ccc_cell2location.yaml"
    script:
        "../scripts/spatial_pp/prepare_refs.py"


rule dot:
    input:
        liver_reference=rules.prepare_refs.output.liver_reference,
        lung_reference=rules.prepare_refs.output.lung_reference,
        heart_reference=rules.prepare_refs.output.heart_reference,
        kidney_reference=rules.prepare_refs.output.kidney_reference,
        liver_spatial=rules.prepare_refs.output.liver_spatial,
        **dot_config["input"]
    output:
        **{name: directory(path) for name, path in dot_config["output"].items()}
    params:
        liver_samples=spatial_pp["liver_samples"],
        **dot_config["params"]
    resources:
        mem_mb = 200000,
        slurm = "gres=gpu:1"
    conda:
        "../envs/ccc_cell2location.yaml"
    script:
        "../scripts/spatial_pp/dot.py"
