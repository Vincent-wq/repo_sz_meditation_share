# Meditation, clinical outcomes and effective connectivity

Analysis code and figure materials for our longitudinal meditation study in schizophrenia.

## Files

- **Analysis source:** the `1a`, `1b`, `3a–3c` and `4a–4b` notebooks and MATLAB scripts in `analysis/`. The notebooks retain their saved figures and statistical results; participant-level table previews are omitted.
- **Figure reproduction:** `analysis/00_reproduce_figure_panels.ipynb` plots results from the aggregate tables in `data/source_data/`. It can run without the participant-level data.

## Reproduce the figures

In a Python environment, run these commands from the repository root:

```text
python -m pip install -r requirements_public.txt
python -m jupyter lab
```

Open `analysis/00_reproduce_figure_panels.ipynb` and run the blocks in order. Set `SAVE_FIGURES = True` to save panels to `fig/public_reproduction/`; the current setting displays them in the notebook.

The tables cover clinical treatment estimates and trajectories, the mediation summary, local-network contrasts, CLPM coefficients, EC trajectories and contrasts, and neural–clinical coefficients and R². CLPM coefficients are shown as effect-size plots; manuscript diagrams were assembled manually. Glass-brain maps and participant scatter plots remain available as saved figures in the source notebooks. An optional omission-range plot summarizes the released influence ranges.

## Data access

Refitting the original analyses requires the participant-level inputs, analysis software and project helper functions. Underlying participant-level data are available from the corresponding authors upon reasonable request. This repository provides analysis source, aggregate figure data and a figure-reproduction notebook.
