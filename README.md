# rdf_compute

A command-line tool for computing per-atom radial distribution functions (RDFs) from molecular dynamics trajectories in XYZ format, powered by [OVITO](https://www.ovito.org/).

---

📄 Author: **Ouail Zakary**  
- 📧 Email: [Ouail.Zakary@oulu.fi](mailto:Ouail.Zakary@oulu.fi)  
- 🔗 ORCID: [0000-0002-7793-3306](https://orcid.org/0000-0002-7793-3306)  
- 🌐 Website: [Personal Webpage](https://cc.oulu.fi/~nmrwww/members/Ouail_Zakary.html)  
- 📁 Portfolio: [GitHub Portfolio](https://ozakary.github.io/)

---

Given a set of **source** atoms (e.g. all Xe atoms, or a specific subset by index), the script computes the RDF between each source atom and a set of **target** atoms (defaulting to all other atoms in the system). Each source atom produces an independent g(r) curve, all written to a single data file and plotted on a single figure.

## Features

- Works with any element — not limited to noble gases or any specific chemistry
- Handles single or multiple source atoms; each gets its own labeled g(r) curve
- Per-atom isolation: when multiple source atoms of the same element are present, other source atoms are masked so they do not contribute to each atom's individual RDF
- Same-element RDFs (e.g. Ru-Ru, P-P) supported by passing the same symbol to both `--source` and `--target`
- Frame striding for fast analysis of long trajectories
- Validated input with clear error messages (unknown elements, invalid indices)
- Progress bar per source atom showing frame-level progress
- Outputs a plain-text data file and a publication-ready PNG plot

## Requirements

- Python 3.9+
- [OVITO Python package](https://www.ovito.org/python/) (`ovito`)
- `numpy`
- `matplotlib`
- `tqdm`

Install dependencies:

```bash
pip install ovito numpy matplotlib tqdm
```

If you are using an Anaconda environment, prefer the OVITO conda package to avoid Qt conflicts:

```bash
conda install --channel ovito ovito
pip install numpy matplotlib tqdm
```

## Usage

```
python rdf_compute.py TRAJECTORY --source ELEMENT [options]
```

### Arguments

**Input**

| Argument | Default | Description |
|---|---|---|
| `TRAJECTORY` | required | Path to XYZ trajectory file (plain or extended XYZ) |
| `--stride N` | 1 | Use every Nth frame. Use larger values for long trajectories |

**Atom selection**

| Argument | Default | Description |
|---|---|---|
| `--source ELEMENT` | required | Element symbol of the source atoms (e.g. `Xe`, `Ru`, `Kr`) |
| `--source-indices IDX [IDX ...]` | all source atoms | Restrict to specific source atoms by 0-based global particle index |
| `--target ELEMENT [ELEMENT ...]` | all other elements | Target element(s) to correlate against. Pass the same symbol as `--source` for same-element RDFs (e.g. `--source Ru --target Ru`) |

**RDF parameters**

| Argument | Default | Description |
|---|---|---|
| `--cutoff ANGSTROM` | 10.0 | RDF cutoff distance in Å |
| `--bins N` | 200 | Number of histogram bins |

**Output**

| Argument | Default | Description |
|---|---|---|
| `--outdir DIR` | `rdf_output/` | Output directory (created automatically if absent) |
| `--prefix STR` | `rdf` | Filename stem; produces `<prefix>_data.dat` and `<prefix>_plot.png` |
| `--dpi N` | 150 | Plot resolution in DPI |

### Output files

`<prefix>_data.dat` — space-delimited text file with columns:

```
# r(Angstrom)  Xe1  Xe2  ...  Xe10
0.006000 0.000000 0.000000 ...
...
```

`<prefix>_plot.png` — all g(r) curves on a single figure, one curve per source atom.

## Examples

**Single source atom vs. all other atoms**

```bash
python rdf_compute.py traj.xyz --source Xe --outdir results/xe
```

**Multiple Xe atoms vs. all other atoms (per-atom RDFs)**

```bash
python rdf_compute.py traj.xyz --source Xe \
    --outdir results/xe_rest \
    --cutoff 12.0 --bins 1000 --stride 10
```

**Xe vs. specific target elements only**

```bash
python rdf_compute.py traj.xyz --source Xe --target O H \
    --outdir results/xe_water --prefix xe_water
```

**Restrict to specific source atoms by index**

```bash
python rdf_compute.py traj.xyz --source Xe --source-indices 0 4 7 \
    --outdir results/xe_selected
```

**Different element — Ru vs. P in a RuP system**

```bash
python rdf_compute.py traj.xyz --source Ru --target P \
    --outdir results/ru_p --prefix ru_p --cutoff 8.0
```

**Kr in a clathrate hydrate**

```bash
python rdf_compute.py traj.xyz --source Kr --target O H \
    --outdir results/kr_water --cutoff 12.0 --bins 300 --stride 5
```

**Same-element RDFs (e.g. Ru-Ru, P-P)**

Pass the same symbol to both `--source` and `--target`. The default target (all other elements) never includes the source element itself, so this must be requested explicitly.

```bash
# Ru-Ru
python rdf_compute.py traj.xyz --source Ru --target Ru \
    --outdir results/ru_ru --prefix ru_ru

# P-P
python rdf_compute.py traj.xyz --source P --target P \
    --outdir results/p_p --prefix p_p

# Xe-Xe (multiple Xe atoms — one shared curve, population-averaged)
python rdf_compute.py traj.xyz --source Xe --target Xe \
    --outdir results/xe_xe --prefix xe_xe
```

## Performance notes

**Stride** is the most effective lever for long trajectories. RDF frames are highly correlated in MD, so `--stride 10` on a 10 000-frame run gives 1000 effectively independent frames at 10x the speed with no meaningful loss of accuracy.

**Large systems (10 000+ atoms):** converting the trajectory from XYZ to LAMMPS binary dump format before analysis reduces I/O time by 5-10x:

```bash
python3 -c "
import warnings; warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
from ovito.io import import_file, export_file
p = import_file('traj.xyz', multiple_frames=True)
export_file(p, 'traj.dump', 'lammps/dump',
            columns=['Particle Type', 'Position.X', 'Position.Y', 'Position.Z'],
            multiple_frames=True)
"
python rdf_compute.py traj.dump --source Xe --cutoff 12.0 --bins 1000 --stride 10
```

**Multiple source atoms** are processed sequentially. Each atom requires one full pass over the selected frames, so total runtime scales linearly with the number of source atoms.

## Notes on element names

The `--source` and `--target` element symbols must match exactly how atom types are named in your XYZ file. If the script cannot find a requested element it exits immediately and prints the type names it detected in frame 0, so there is no ambiguity.

## License

MIT
