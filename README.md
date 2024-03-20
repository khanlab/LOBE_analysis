# lobe_analysis


## before using, ensure you have pipx and poetry:


install pipx
```
python3 -m pip install --user pipx
python3 -m pipx ensurepath
```

install poetry
```
pipx install poetry
```


1. Clone repository


2. If using poetry/pipx:

```
poetry install
poetry shell
```

3. otherwise:

```
python3.11 -m venv venv
source venv/bin/activate
pip install .
```


### To QC fmriprep data from zip files:

Use the -l (list) and pipe to less to find the filepaths you want to unzip, e.g.:
```
unzip -l  /home/ROBARTS/alik/graham/LOBE/derivatives/snakebatch_LOBE_2024_03_19.zip  | less
```


Then, on /localscratch run:
```
unzip   /home/ROBARTS/alik/graham/LOBE/derivatives/snakebatch_LOBE_2024_03_19.zip  "LOBE/derivatives/fmriprep_23.1.0/sub-*.html" "LOBE/derivatives/fmriprep_23.1.0/sub-*/figures/*"
```

