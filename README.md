# FCNC NanoAODv15 2024

Required release: CMSSW_13_3_3.

- LeptonFR: lepton fake rate
- LeptonSF: lepton scale factors
- TriggerSF: trigger scale factors
- central: central analysis

From an empty CMSSW_13_3_3/src directory, run:

    git clone -b central --single-branch https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation.git .

## 2024 three-lepton analysis workflow

This analysis uses CMS Run 3 NanoAODv15 data and simulation from 2024. It selects events with three leptons (electrons or muons) and a Z boson candidate, compares data with Monte Carlo (MC) simulation in four channels, and estimates the nonprompt lepton (NPL) contribution.

The `central` branch includes electron and muon selection scale factors (SFs) and fake-rate (FR) maps. **To use the supplied corrections, follow Sections 1, 3, and 4.** To measure new corrections, complete Section 2 before producing the central ntuples.

### 1. Set up the environment and input file lists

The clone command above places `TZWi/` directly under `CMSSW_13_3_3/src/`. To repeat a correction measurement, use a separate empty CMSSW work area and replace `central` in the clone command with the measurement branch name. Build each checkout and enter its analysis directory.

```bash
# Run from the src directory of each CMSSW_13_3_3 work area.
cmsenv
scram b -j 4
cd "$CMSSW_BASE/src/TZWi/TopAnalysis/test/fcncTriLepton"
```

Batch production requires `nano_postproc.py` from NanoAODTools, HTCondor, and [`create-batch3_mine`][batch]. The measurement and analysis scripts use PyROOT and PyYAML. The central analysis also uses `correctionlib` to read its JSON corrections.

In the analysis directory, `NanoAOD/2024/` points to `TZWi/TopAnalysis/data/datasets/NanoAOD/2024/`. The YAML files define the datasets. The `.txt` files in the subdirectories list one input ROOT file per line. Check that the batch workers can open the supplied `/new_xrootd/...` and `file:///new_xrootd/...` paths. At another site, update the lists with accessible paths or XRootD URLs. [2024 dataset lists][datasets]

If the datasets change, update the corresponding 2024 YAML files on the branch you are using. For example, on `central`, a site with the same `/new_xrootd` mount can regenerate the file lists from DAS with the command below. The script uses a CMS grid proxy and rewrites both the YAML files and the text lists.

```bash
(
  cd "$CMSSW_BASE/src/TZWi/TopAnalysis/data/datasets/NanoAOD/2024" || exit 1
  python3 "$CMSSW_BASE/src/TZWi/TopAnalysis/scripts/tzwi-updatedataset__new_xrootd" \
    MC.RunIISummer24.central.yaml RD.Run2024.yaml
)
```

### 2. Measure new corrections

Each row below is an independent measurement. Run the submission command from `TZWi/TopAnalysis/test/fcncTriLepton/` on the indicated branch. Once its batch jobs finish, run the extraction command from the same directory. The `makentuple*` scripts extract the ROOT files from each job's `result_*.tgz` archive.

| Branch | Default samples | Submit jobs | Extract completed jobs |
| --- | --- | --- | --- |
| [`TriggerSF`][trig-submit] | `Run2024.JetMET` data and `WZto3LNu` MC | `python3 01.1_submit_trigsf_2024.py` | `(cd submit_2024_trigSF && bash "$CMSSW_BASE/src/TZWi/TopAnalysis/scripts/makentuple2024trigSF")` |
| [`LeptonSF`][sf-submit] | EGamma, Muon, and MuonEG data, `WZto3LNu` signal MC, and background MC | `python3 01.1_submit_leptonSF.py` | `(cd submit_2024_leptonSF && bash "$CMSSW_BASE/src/TZWi/TopAnalysis/scripts/makentuple2024leptonSF")` |
| [`LeptonFR`][fr-submit] | `Run2024.JetMET` data | `python3 01.1_submit_FR.py` | `(cd submit_2024_FR && bash "$CMSSW_BASE/src/TZWi/TopAnalysis/scripts/makentuple2024FR")` |

#### Trigger scale factors

Among events that pass the MET/MHT reference trigger, measure the fraction that also pass the analysis lepton triggers. The trigger SF is this efficiency in data divided by the corresponding efficiency in MC. The command below processes all four channels and removes duplicate data events. Keep `{channel}` in both quoted paths so the script can substitute each channel name. [Trigger SF measurement script][trig]

```bash
python3 make_trigger_sf_2024.py \
  --data-glob 'submit_2024_trigSF/ntuple_2024_trigSF/reco/{channel}/JetMET*Run2024*/*.root' \
  --mc-glob 'submit_2024_trigSF/ntuple_2024_trigSF/reco/{channel}/*WZto3LNu*/*.root' \
  --deduplicate-data
```

The default outputs are `triggerSF_2024.root`, `triggerSF_2024.csv`, and `triggerSF_2024.json`. Add `--make-2d` to produce two-dimensional efficiency and SF maps. The default MC event weight on `central` omits the trigger SF. To apply it, include the chosen per-channel value or two-dimensional map in the event-weight calculation.

#### Lepton selection scale factors

In `LepSFEl` and `LepSFMu`, a tight, opposite-sign, same-flavor lepton pair forms the Z candidate. The remaining lepton is the probe. Subtract background and nonprompt WZ probe contributions from both the passing and total probe counts in data. Divide the resulting data efficiency by the prompt-probe efficiency in WZ MC. The default promptMVA thresholds are `> 0.90` for electrons and `> 0.70` for muons. [Sample and selection configuration][sf-config] · [Lepton SF measurement script][sf-run]

```bash
LEPSF_BASE="$PWD/submit_2024_leptonSF/ntuple_2024_leptonSF" \
  bash run_make_lepton_sf_2024.sh
```

The script writes `leptonSF_electron_2024_signedEta.root` and `leptonSF_muon_2024_signedEta.root`, each with a matching CSV file. The map axes are probe transverse momentum ($p_T$) and signed pseudorapidity ($\eta$). Electrons use supercluster $\eta$.

If you change the input MC samples, recompute `sum_gen_weight` from `Runs.genEventSumw` in the **original NanoAOD files**, then run the measurement with the updated YAML configuration.

```bash
python3 fill_leptonSF_sum_gen_weight.py \
  --config leptonSF_samples_2024.yaml --dataset-base NanoAOD/2024 \
  --output leptonSF_samples_2024_filled.yaml

CONFIG=leptonSF_samples_2024_filled.yaml \
  LEPSF_BASE="$PWD/submit_2024_leptonSF/ntuple_2024_leptonSF" \
  bash run_make_lepton_sf_2024.sh
```

#### Lepton fake rates

The `FREl` and `FRMu` channels use JetMET data selected by the OR of the configured single-jet triggers from `PFJet200` through `PFJet550`. The default measurement region requires a selected trigger jet with $p_T > 250$ GeV and $p_T^{\mathrm{miss}} + m_T < 50$ GeV. For electrons and muons separately, measure the fraction of loose probes that pass the tight selection in bins of probe $p_T$ and signed $\eta$. [FR ntuple production script][fr-prod] · [FR map generation script][fr-run]

```bash
INPUT_BASE="$PWD/submit_2024_FR/ntuple_2024_FR/reco" \
  bash run_make_fake_rates_2024_signedEta.sh
```

The script writes ROOT and CSV files to `fakerate_2024/` for each flavor and each of three categories. `Inclusive` includes all b-jet multiplicities, `0B` requires `FR_nBJet == 0`, and `GE1B` requires `FR_nBJet >= 1`. The default maps use JetMET data without prompt-MC subtraction. The FR production script also leaves the data unfiltered by the Golden JSON. For a measurement restricted to certified runs and luminosity sections, apply the mask to the inputs or during ntuple production.

#### Use the new corrections in the central analysis

| Correction | Files from the measurement | Destination and use on `central` |
| --- | --- | --- |
| Electron and muon selection SFs | Two `leptonSF_{electron,muon}_2024_signedEta.root` files | Copy to `TZWi/TopAnalysis/data/leptonSF/2024/`. The [lepton SF module][sf-apply] reads them when processing MC. |
| Electron and muon FRs | Six `fakerate_{electron,muon}_promptMVA_2024_JetMET_SignedEta_{Inclusive,0B,GE1B}.root` files | Copy to `TZWi/TopAnalysis/data/leptonFR/2024/`. The [FR module][fr-apply] reads the selected maps when processing NPL channels. |

Copy the new files **before producing the central ntuples**. Ntuple production writes the SF and FR weights to per-event branches.

### 3. Produce the central analysis ntuples

Run the following commands from `TZWi/TopAnalysis/test/fcncTriLepton/` in the `central` work area. The [submission script][central-submit] reads [`config/grouping.yaml`][grouping] and the input file lists, then runs [`01_prod_ntuple.sh`][central-prod]. It processes `ElElEl`, `ElMuMu`, `MuElEl`, and `MuMuMu`, along with the corresponding nonprompt channels prefixed with `NPL`.

```bash
cd "$CMSSW_BASE/src/TZWi/TopAnalysis/test/fcncTriLepton"
DRY_RUN=1 python3 01.1_submit.py
python3 01.1_submit.py
```

Use the dry run to check the channels, datasets, and job counts before submitting. The default submission directory is `submit_2024/`. Set `NFILE` to change the number of input files per job (default 10).

The default NPL weights use the `Inclusive` FR maps. To produce a separate sample with maps split by b-jet multiplicity, run `FAKERATE_MAP_MODE=split SUBMIT_DIR=submit_2024_split python3 01.1_submit.py`. Extract each production in its own submission directory because the two sets of archives use the same internal ntuple paths.

The baseline selection requires exactly three electrons or muons in total. In descending $p_T$ order, the leptons must have $p_T$ of at least 25, 10, and 10 GeV, respectively. The Z candidate must have zero net charge and satisfy $|m_{\ell\ell}-91.2\,\mathrm{GeV}|\leq15\,\mathrm{GeV}$. Data production applies the 2024 Golden JSON. MC production applies electron reconstruction SFs, electron and muon selection SFs, pileup weights, and UParTAK4 b-tag weights. Additional selections on missing transverse momentum, W transverse mass, and b-jet multiplicity can be compared at the histogram stage. [Three-lepton selection code][cutflow]

After the default batch jobs finish, run the commands below. The failure checker looks for reported failures inside existing archives in `submit_2024/`. Also check that the submitted jobs produced the expected `result_*.tgz` archives and ROOT files. Jobs with an explicitly selected FR map mode have a `.FRMap_*` suffix that falls outside the checker's default directory pattern. Check those jobs separately in their submission directory.

```bash
bash 01.2_checkFailed_2024_parallel.sh
(cd submit_2024 && bash "$CMSSW_BASE/src/TZWi/TopAnalysis/scripts/makentuple2024")
```

The extracted files are under `submit_2024/ntuple_2024/reco/<channel>/<dataset>/`. The histogram scripts read `ntuple_2024/` from the analysis directory. If the name `ntuple_2024` is unused in that directory, create the symbolic link below.

```bash
ln -sT submit_2024/ntuple_2024 ntuple_2024
```

### 4. Produce histograms and data–MC comparison plots

Run the following scripts in order from the same `central` analysis directory. The first script deletes and rebuilds `raw_hist_2024/`. Move any previous histograms you want to compare to a separate directory before rerunning it.

```bash
python3 02_make_histograms_2024.py
python3 03_scalemerge_2024.py
python3 04_drawPlots_2024.py
```

| Stage | Configuration and output |
| --- | --- |
| Fill histograms | Use the variables and selections in [`config/histogramming.yaml`][hist-config] to write per-dataset ROOT files under `raw_hist_2024/<channel>/<process>/`. |
| Normalize and merge | [`03_scalemerge_2024.py`][merge] takes cross sections from `config/crosssection.yaml`. It estimates the nonprompt contribution by subtracting normalized prompt MC from FR-weighted data and writes the channel results to `hist_2024/<channel>.root`. |
| Draw plots | Use `config/plots.yaml` to produce `plots/<channel>/<selection step>/<variable>.png`. Open `plots/index.html` to browse the plots. |

The default MC event weight is `genWeight * Electron_SF * Electron_RecoSF * Muon_SF * puWeight * BtagWeight`. MC normalization uses the configured integrated luminosity of $109{,}000\,\mathrm{pb}^{-1}$ and the MC cross section. The [histogram tool][hist-tool] sums `Runs.genEventSumw` over the skim ROOT files for the normalization denominator. If this branch is missing, it falls back to the sum of `Events.genWeight` after selection. If the tool reports this fallback, check the sum of weights for the original generated sample.

For the Run 2 production and correction workflows, see [FCNC UL Run 2 variations](https://github.com/panthera718/FCNC_20ULRun2_Variation#readme) and [FCNC UL Run 2 Tau Fake Rate](https://github.com/panthera718/FCNC_20ULRun2_TauFakeRate#readme).

[batch]: https://github.com/panthera718/hep-tools__create-batch
[datasets]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/tree/central/TZWi/TopAnalysis/data/datasets/NanoAOD/2024
[trig]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/TriggerSF/TZWi/TopAnalysis/test/fcncTriLepton/make_trigger_sf_2024.py
[sf-config]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/LeptonSF/TZWi/TopAnalysis/test/fcncTriLepton/leptonSF_samples_2024.yaml
[sf-run]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/LeptonSF/TZWi/TopAnalysis/test/fcncTriLepton/run_make_lepton_sf_2024.sh
[fr-prod]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/LeptonFR/TZWi/TopAnalysis/test/fcncTriLepton/01_prod_ntuple_FR.sh
[fr-run]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/LeptonFR/TZWi/TopAnalysis/test/fcncTriLepton/run_make_fake_rates_2024_signedEta.sh
[sf-apply]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/central/TZWi/TopAnalysis/python/postprocessing/lepSFProducer.py
[fr-apply]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/central/TZWi/TopAnalysis/python/postprocessing/FakeRateProducer.py
[central-submit]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/central/TZWi/TopAnalysis/test/fcncTriLepton/01.1_submit.py
[central-prod]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/central/TZWi/TopAnalysis/test/fcncTriLepton/01_prod_ntuple.sh
[grouping]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/central/TZWi/TopAnalysis/test/fcncTriLepton/config/grouping.yaml
[cutflow]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/central/TZWi/TopAnalysis/python/postprocessing/fcncTriLeptonCutFlow.py
[hist-config]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/central/TZWi/TopAnalysis/test/fcncTriLepton/config/histogramming.yaml
[merge]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/central/TZWi/TopAnalysis/test/fcncTriLepton/03_scalemerge_2024.py
[hist-tool]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/central/TZWi/TopAnalysis/scripts/tzwi-makehistograms
[trig-submit]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/TriggerSF/TZWi/TopAnalysis/test/fcncTriLepton/01.1_submit_trigsf_2024.py
[sf-submit]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/LeptonSF/TZWi/TopAnalysis/test/fcncTriLepton/01.1_submit_leptonSF.py
[fr-submit]: https://github.com/panthera718/FCNC_NanoAODv15_2024_Variation/blob/LeptonFR/TZWi/TopAnalysis/test/fcncTriLepton/01.1_submit_FR.py
