#!/usr/bin/env python3
"""DeepDock: Iterative docking and selection of molecules using gnina and SPRINT."""

import os
import sys
# the default mimalloc is a hog and will get killed by slurm unless you reserve a lot more memory than you need
os.environ['ARROW_DEFAULT_MEMORY_POOL'] = 'jemalloc'

import pyarrow as pa
import rdkit
from rdkit.Chem import AllChem as Chem
import gzip, glob, re
import numpy as np
import matplotlib.pyplot as plt
import ultrafast
from ultrafast import embed
import torch
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader
from ultrafast.datamodules import EmbedDataset, embed_collate_fn
from ultrafast.model import DrugTargetCoembeddingLightning
from ultrafast.utils import get_featurizer
import torch.nn.functional as F
import os
import uuid
import shutil
import tempfile
import logging
from rdkit import DataStructs
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Descriptors import MolWt, ExactMolWt
import xgboost
from functools import partial
import multiprocessing
import dask.config as dc
import dask
import dask.bag as db
from distributed import WorkerPlugin
import dask.dataframe as dd
import shlex
import yaml
import pyarrow.parquet as pq
from dask.distributed import print as dask_print
from dask_jobqueue import *
from dask.distributed import Client, LocalCluster
from torch.nn import CosineSimilarity
from rdkit.Chem import AllChem as Chem
from io import StringIO
import math
import os
import tempfile
import subprocess
from pathlib import Path
from typing import Iterable, Union
import pickle
from dask.distributed import get_worker
import pyarrow.dataset as ds
from dask.distributed import Client, Semaphore, get_client, wait
from scipy.stats import gaussian_kde
import argparse
from Bio import SeqIO
import resource



def makesdf(smi, name, molweight_cutoff=1200):
    """Generate low energy 3D conformation and return as sdf string.
    Returns empty string on failure or if above molweight cutoff.
    @param smi: input smiles
    @param name: molecule name
    @param molweight_cutoff: maximum allowed molecular weight
    @return: sdf string or empty string on failure"""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return ""
        if ExactMolWt(mol) > molweight_cutoff:
            return ""
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol)
        mol.SetProp("_Name", name)
        params = Chem.ETKDG()
        params.timeout = 100
        cids = Chem.EmbedMultipleConfs(mol, 10, Chem.ETKDG())

        if not cids:
            return ""
        cenergy = []
        for conf in cids:
            converged = not Chem.UFFOptimizeMolecule(mol, confId=conf)
            cenergy.append(Chem.UFFGetMoleculeForceField(mol, confId=conf).CalcEnergy())

        mol = Chem.RemoveHs(mol)
        sortedcids = sorted(cids, key=lambda cid: cenergy[cid])
        sio = StringIO()
        sdwriter = Chem.SDWriter(sio)
        sdwriter.write(mol, sortedcids[0])
        sdwriter.close()
        return sio.getvalue()
    except:
        return ""


def dock_partition(
    df: pd.DataFrame,
    args: Union[str, Iterable[str]] = (),
    gnina_executable: str = "gnina",
    gnina_timeout: int = 900,
    molweight_cutoff=1200,
    partition_info=None,
) -> str:
    """Dock molecules in dataframe using gnina.
    
    Args:
        df: DataFrame with 'smiles' and 'name' columns
        args: Additional gnina command-line arguments
        gnina_executable: Path to gnina executable
        gnina_timeout: Timeout in seconds per molecule
        molweight_cutoff: Maximum allowed molecular weight
        partition_info: Dask partition info for logging (auto-provided by dask)
    
    Returns:
        DataFrame with columns: mol, smiles, minimizedAffinity, CNNscore, 
        CNNaffinity, CNN_VS. Empty DataFrame on failure.
    """
    in_path = out_path = None

    logger = logging.getLogger()
    level = logger.getEffectiveLevel()
    if partition_info:
        logger = logging.Logger("dock_partition")
        logger.setLevel(level)
        log_path = f"dask_logs/partition.{partition_info['number']}.log"
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        handler = logging.FileHandler(log_path)
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
        )
        logger.addHandler(handler)

    logger.info(f"docking partition {partition_info['number']}")
    try:
        # Normalize args to a list of strings
        if isinstance(args, str):
            extra = shlex.split(args)
        else:
            extra = list(args) if args else []

        sdfs = ""
        for r, row in df.iterrows():
            smi = row.smiles
            name = row["name"]
            sdf = makesdf(smi, f"{name} {smi}", molweight_cutoff=molweight_cutoff)
            if sdf:
                sdfs += sdf
            else:
                logger.warning(f"ERROR with {smi} {name}")

        in_fd = out_fd = None

        # Create input SDF temp file
        in_fd, in_path = tempfile.mkstemp(suffix=".sdf", prefix="gnina_in_")
        with os.fdopen(in_fd, mode="w", encoding="utf-8", newline="\n") as f:
            f.write(sdfs)
        # fd is now closed by context manager; avoid double-close in finally
        in_fd = None

        # Create output temp file path (don’t open; gnina will write it)
        out_fd, out_path = tempfile.mkstemp(suffix=".sdf", prefix="gnina_out_")
        os.close(out_fd)  # Close immediately; gnina will overwrite
        out_fd = None

        cmd = [
            gnina_executable,
            *extra,
            "-l",
            in_path,
            "-o",
            out_path,
            "--seed",
            "0",
            "--cpu",
            "1",
        ]
        timeout_seconds = gnina_timeout * len(df)
        try:
            # Run gnina, capture stderr for helpful errors
            proc = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,  # gnina usually writes to files; keep stdout just in case
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout_seconds,
            )
        except subprocess.TimeoutExpired as e:
            logger.warning(
                f"Command timed out after {timeout_seconds} seconds: {' '.join(cmd)}"
            )
            logger.warning(
                f"Stdout captured before timeout: {e.stdout.decode() if e.stdout else ''}"
            )
            logger.warning(
                f"Stderr captured before timeout: {e.stderr.decode() if e.stderr else ''}"
            )

        except subprocess.CalledProcessError as e:
            logger.error(f"Error running gnina: {' '.join(cmd)}")
            logger.error(proc.stdout)
            logger.error(proc.stderr)

        # Read result
        with open(out_path, "r", encoding="utf-8", errors="replace") as f:
            result = f.read()

        # take first output for each unique name
        seen = set()
        results = []
        for mol in result.split("$$$$\n"):
            m = re.search(
                r"^(\S+) (\S+)\n.*?> <minimizedAffinity>\n(\S+)\n\n> <CNNscore>\n(\S+)\n\n> <CNNaffinity>\n(\S+)\n\n> <CNN_VS>\n(\S+)\n.*",
                mol,
                re.DOTALL | re.MULTILINE,
            )
            if not m:
                logger.warning(f"Failed to match: {mol}")
                continue
            else:
                try:
                    name = m.group(1)
                    smiles = m.group(2)
                    if m and name not in seen:
                        seen.add(name)
                        vina = float(m.group(3))
                        cnnscore = float(m.group(4))
                        aff = float(m.group(5))
                        vs = float(m.group(6))
                        results.append((mol, smiles, vina, cnnscore, aff, vs))
                except Exception as e:
                    logger.error(
                        f"An unexpected error happened during sdf parsing: {e}\n{df}"
                    )
                    continue

        # make dataframe
        ret = pd.DataFrame(
            results,
            columns=[
                "mol",
                "smiles",
                "minimizedAffinity",
                "CNNscore",
                "CNNaffinity",
                "CNN_VS",
            ],
        )
        return ret

    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}\n{df}")
        return pd.DataFrame(
            columns=[
                "mol",
                "smiles",
                "minimizedAffinity",
                "CNNscore",
                "CNNaffinity",
                "CNN_VS",
            ]
        )
    finally:
        # Clean up temp files
        for p in (in_path, out_path):
            try:
                if p and Path(p).exists():
                    os.remove(p)
            except Exception:
                # Best-effort cleanup; ignore failures
                pass


def dock_batch(outdir, batch, gcmd, cluster):
    """Dock a batch of molecules using gnina.
    
    Overwrites existing results if already present. Requires molecules to have been 
    selected and saved via initial_select() or select_next_batch().
    
    Args:
        outdir: Output directory containing batch subdirectories
        batch: Batch number to dock
        gcmd: Gnina command-line arguments (additional flags to pass to gnina)
        cluster: Dask cluster for distributed docking
    
    Returns:
        None. Results are written to outdir/batch{batch}/docked.parquet
    """
    with Client(cluster) as client:
        # split and dock
        combined = dd.read_parquet(
            f"{outdir}/batch{batch}/selected.parquet", split_row_groups=True
        )
        meta = [
            ("mol", "object"),
            ("smiles", "object"),
            ("minimizedAffinity", "f4"),
            ("CNNscore", "f4"),
            ("CNNaffinity", "f4"),
            ("CNN_VS", "f4"),
        ]
        docked = combined.map_partitions(dock_partition, gcmd, meta=meta)
        docked.to_parquet(f"{outdir}/batch{batch}/docked.parquet", overwrite=True)


def initial_select(
    infile: str,
    prefix: str,
    protein_sequence: str,
    sprint_checkpoint: str,
    cluster,
    iolimit,
    device: str = 'cpu',
    N: int = 100_000,
):
    """Select top N molecules by SPRINT score and random sample from input file.
    
    Selects the highest-scoring molecules using the SPRINT model with the provided
    protein sequence, combined with a random sample of the same size. Results are
    written to the output directory.
    
    Args:
        infile: Input parquet file path with 'sprint' embeddings and 'name' columns
        prefix: Output directory prefix (results saved to prefix/batch0/)
        protein_sequence: Target protein sequence string
        sprint_checkpoint: Path to SPRINT model checkpoint file
        cluster: Dask cluster for distributed processing
        N: Number of molecules to select (default: 100,000)
    
    Returns:
        None. Results are written to {prefix}/batch0/selected.parquet
    """
    outdir = f"{prefix}/batch0"
    os.makedirs(outdir, exist_ok=True)


    if not sprint_checkpoint:
        model = None
        T = None
    else:
        torch.set_grad_enabled(False)
        model = DrugTargetCoembeddingLightning.load_from_checkpoint(sprint_checkpoint)
        model.eval().to(device)
        featurizer = get_featurizer(model.args.target_featurizer, batch_size=1)
        e = featurizer._transform_single(protein_sequence).to(device)

        # check featurizer output is not all zeros
        if e.abs().sum().item() == 0:
            raise ValueError("Error: Target featurizer produced all-zero features.")
        T = model.embed(e.unsqueeze(0), sample_type="target")

    def score(D):
        return F.cosine_similarity(T, torch.tensor(D), dim=0).item()

    with Client(cluster) as client:
        sem = Semaphore(max_leases=iolimit, name=f"Limiter{iolimit}")

        #calculate size for random sampling
        dataset = ds.dataset(infile)
        npart = len(list(dataset.get_fragments()))
        total_rows = sum(f.metadata.num_rows for f in dataset.get_fragments())
        fraction = N/total_rows

        def iter_parquet_batches(fname, batch_size):

            with sem:
                # Get Dask temporary directory (falls back to system temp if not set)
                worker = get_worker()
                tmpdir = worker.local_directory
                if tmpdir is None:
                    # Absolute fallback: system temp directory
                    tmpdir = tempfile.gettempdir()

                # Construct destination path
                local_name = os.path.basename(fname)
                tmp_path = os.path.join(tmpdir, f"{uuid.uuid4()}_{local_name}")

                # Copy file if not already present (or you may force overwrite)
                try:
                    shutil.copy2(fname, tmp_path)
                except Exception:
                    dask_print(f'Failed copy {fname}')
                    tmp_path = fname
                
            table = pq.ParquetFile(fname)
            it = table.iter_batches(batch_size=batch_size, use_threads=False, columns=["smiles", "name"])
            i = 0
            while True:
                try:
                    i += 1
                    batch = next(it)          # I/O happens here
                except StopIteration:
                    table.close()
                    try:
                        os.remove(tmp_path)
                    except Exception:
                        pass
                    return
                df = batch.to_pandas()
                yield df

        def make_df(fname):
            dfs = []
            for df_chunk in iter_parquet_batches(fname, batch_size=1000):
                if "sprint" not in df_chunk.columns:
                    df_chunk = ensure_sprint_column(df_chunk, model, device=device)
                df_chunk["score"] = df_chunk["sprint"].map(score)
                dfs.append(df_chunk[["name", "smiles", "score"]])
            df = pd.concat(dfs, ignore_index=True)
            sample = df.sample(frac=fraction).copy().assign(score=0)
            return df.nlargest(N, "score"), sample

        def reduce_df(df1, df2):
            df1_top, df1_sample = df1
            df2_top, df2_sample = df2
            return (
                pd.concat([df1_top, df2_top], ignore_index=True).nlargest(N, "score"),
                pd.concat([df1_sample, df2_sample], ignore_index=True),
            )

        if model is None:  # no sprinting, just random sample
            df = dd.read_parquet(infile, columns=["name", "smiles"], memory_map=True, pre_buffer=False)
            combined = df[["name", "smiles"]].sample(N / len(df)).assign(score=0).compute()
        else:
            # get top sprint scoring and random sample and checkpoint to disk
            files = glob.glob(f"{infile}/*.parquet")
            b = db.from_sequence(files, npartitions=len(files))
            topdf, rand = b.map(make_df).fold(reduce_df).compute()
            combined = pd.concat([topdf, rand])

        combined = combined.drop_duplicates("smiles")

        npart = min(2 * N // 10, npart)
        combined = dd.from_pandas(combined, npart)

        combined.to_parquet(
            f"{outdir}/selected.parquet",
            overwrite=True,
            row_group_size=10,
            engine="pyarrow",
        )


def smiles_to_fp(smi, fpgen):
    mol = Chem.MolFromSmiles(smi)
    fp = fpgen.GetFingerprint(mol)
    arr = np.zeros((fp.GetNumBits(),), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def smiles_to_fp_binary(smi, fpgen):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = fpgen.GetFingerprint(mol)
    return fp.ToBinary()


def ensure_fp_column(pdf: pd.DataFrame) -> pd.DataFrame:
    if "fp" not in pdf.columns:
        fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        pdf["fp"] = pdf["smiles"].map(lambda smi: smiles_to_fp_binary(smi, fpgen))
        pdf = pdf[pdf["fp"].notna()].copy()
    return pdf


def ensure_sprint_column(pdf: pd.DataFrame, model, device: str = "cpu") -> pd.DataFrame:
    if "sprint" not in pdf.columns:
        pdf = ensure_fp_column(pdf)
        if len(pdf) == 0:
            return pdf
        fps = []
        for row in pdf.itertuples(index=False):
            fp = DataStructs.ExplicitBitVect(row.fp)
            arr = np.zeros((fp.GetNumBits(),), dtype=np.float32)
            DataStructs.ConvertToNumpyArray(fp, arr)
            fps.append(arr)
        with torch.no_grad():
            tensor = torch.from_numpy(np.stack(fps)).to(device)
            emb = model.embed(tensor)
        pdf["sprint"] = [x.astype(np.float32).tolist() for x in emb.cpu().numpy()]
    return pdf


def fp_partition(s):
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    return np.stack(s["smiles"].apply(lambda x: smiles_to_fp(x, fpgen)).values)


def predict_part(pdf: pd.DataFrame, n: int, model, seenpath: str, partition_info=None) -> pd.DataFrame:
    '''Predict molecular binding scores for a partition of compounds using a trained model.

    This function processes a dataframe partition containing molecular fingerprints,
    generates predictions using the provided XGBoost model, and returns the top-scoring
    compounds.

    Args:
        pdf (pd.DataFrame): Input dataframe partition containing columns:
            - 'smiles': SMILES string representation of molecules
            - 'fp': Molecular fingerprint data
            - 'name': Compound identifier
        n (int): Number of top-scoring results to return, sorted by score descending
        model: Trained XGBoost model for molecular property prediction
        seenpath: path to pickle of set of seen smiles
        partition_info (dict, optional): Metadata about the Dask partition provided
            automatically by Dask, containing 'number' key for partition identifier.
            If provided, enables partition-specific logging. Defaults to None.

    Returns:
        pd.DataFrame: Dataframe containing the top n compounds with columns:
            - 'name': Compound identifier
            - 'smiles': SMILES string representation
            - 'score': Predicted binding score
    '''

    logger = logging.getLogger()
    level = logger.getEffectiveLevel()
    number = 0
    if partition_info:
        logger = logging.Logger("dock_partition")
        logger.setLevel(level)
        log_path = f"dask_logs/partition.{partition_info['number']}.log"
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        handler = logging.FileHandler(log_path)
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
        )
        logger.addHandler(handler)
        number = partition_info["number"]

    with open(seenpath,'rb') as inset:
        seen = pickle.load(inset)
    logger.info(f"About to read {number}")
    pdf = ensure_fp_column(pdf)
    if len(pdf) == 0:
        return pd.DataFrame(columns=["name", "smiles", "score"])
    fpsize = 2048
    if len(pdf):
        fpsize = DataStructs.ExplicitBitVect(pdf.iloc[0]["fp"]).GetNumBits()
    fps = np.zeros((len(pdf), fpsize), dtype=np.float32)
    mask = np.zeros(len(pdf),dtype=bool)
    for r, row in enumerate(pdf.itertuples(index=False)):
        mask[r] = row.smiles not in seen
        fp = DataStructs.ExplicitBitVect(row.fp)
        DataStructs.ConvertToNumpyArray(fp, fps[r])

    logger.info(f"About to predict {number} with shape {fps.shape}")
    pdf["score"] = model.predict(fps).astype(np.float32)

    del fps  # try to trigger garbage collection
    pdf = pdf.drop(
        columns=[c for c in pdf.columns if c not in ["name", "smiles", "score"]]
    )

    #wait until here to mask out to reduce memory, assuming most rows are kept
    pdf = pdf[mask]
    logger.info(f"Finished predicting {number} ({len(pdf)} unique)")

    return pdf.nlargest(n, "score")


def select_next_batch(
    infile,
    outdir,
    curbatch,
    cluster,
    score="CNN_VS",
    N=100_000,
    recompute_model=False,
    iolimit=50,
    local_workers=16
):
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    sz = fpgen.GetOptions().fpSize

    modelpath = f"{outdir}/batch{curbatch}/model.json"
    seenpath = f"{outdir}/batch{curbatch}/seen.pkl"
    seen = set()
    if os.path.exists(modelpath) and os.path.exists(seenpath) and not recompute_model:
        model = xgboost.XGBRegressor()
        model.load_model(modelpath)
        with open(seenpath,'rb') as inset:
            seen = pickle.load(inset)
    else:
        # do model training locally
        npart = multiprocessing.cpu_count()  # TODO, better way to set this
        with LocalCluster(threads_per_worker=1,n_workers=local_workers, memory_limit="40GB") as lcluster, Client(lcluster) as client:
            dfs = []
            for i in range(0, curbatch + 1):
                df = dd.read_parquet(
                    f"{outdir}/batch{i}/docked.parquet/", columns=["smiles", score]
                ).repartition(npartitions=npart)
                dfs.append(df)
            docked = dd.concat(dfs)
            seen = set(docked.smiles.compute())
            print(f"Training model with {len(seen)} unique SMILES")
            with open(seenpath,'wb') as out:
                pickle.dump(seen,out)
            X = docked.map_partitions(
                fp_partition, meta=np.empty((0, sz), dtype=np.float32)
            )
            y = docked[score].to_dask_array()

            X = X.compute()
            y = y.compute()

            model = xgboost.XGBRegressor(max_depth=5,n_jobs=-1)
            model.fit(X, y)
            train_score = model.score(X, y)
            print(f"Model training R² score: {train_score:.4f}")
            model.save_model(modelpath)

    with Client(cluster) as client:
        sem = Semaphore(max_leases=iolimit, name=f"Limiter{iolimit}")

        # dask insists on trying to load the whole parquet before processing,
        # so manually implement the map and reduction
        def make_df(fname):
            with sem:
                try:
                    df = pd.read_parquet(fname, columns=["smiles", "name", "fp"])
                except (KeyError, ValueError, pa.lib.ArrowInvalid):
                    df = pd.read_parquet(fname, columns=["smiles", "name"])
            return predict_part(df, n=N, model=model, seenpath=seenpath)

        def reduce_df(df1, df2):
            return pd.concat([df1, df2], ignore_index=True).nlargest(N, "score")

        files = glob.glob(f"{infile}/*.parquet")

        b = db.from_sequence(files, npartitions=len(files))
        topdf = b.map(make_df).fold(reduce_df).compute()  # parallel reduce

        topdf = topdf.drop_duplicates("smiles")

        topdf = topdf[~topdf["smiles"].isin(seen)]
        print(f"Selected {len(topdf)} for batch {curbatch+1}")

        maxpart = 2 * N // 10
        topdf = dd.from_pandas(topdf, maxpart)
        os.makedirs(f"{outdir}/batch{curbatch+1}/", exist_ok=True)
        topdf.to_parquet(
            f"{outdir}/batch{curbatch+1}/selected.parquet",
            overwrite=True,
            row_group_size=10,
            engine="pyarrow",
        )


def get_batch_dirs(bench_dir):
    if not os.path.isdir(bench_dir):
        raise FileNotFoundError(f"Missing bench directory: {bench_dir}")

    # Find batchN directories in numeric order starting at N=0
    batch_dirs = []
    for name in os.listdir(bench_dir):
        m = re.fullmatch(r"batch(\d+)", name)
        if m:
            batch_dirs.append((int(m.group(1)), name))
    batch_dirs.sort(key=lambda t: t[0])

    if not batch_dirs:
        raise FileNotFoundError(f"No batch<N> directories found under {bench_dir}")

    return batch_dirs


def next_batch(benchdir):
    """Get next batch index. This is the lowest numbered batch
    that has not been docked."""
    batch_dirs = get_batch_dirs(benchdir)
    lastb = batch_dirs[-1][0]
    nextb = lastb + 1
    if not os.path.exists(os.path.join(benchdir, f"batch{lastb}/docked.parquet")):
        nextb = lastb
    return nextb

def plot_dists_by_batch(
    bench_dir,
    features=("CNN_VS", "CNNscore", "CNNaffinity","minimizedAffinity"),
    bins=50,
    figsize=(16, 4),
    style="seaborn-v0_8-colorblind",
    limits=None,
    local_workers=16):
    """Plot molecular property distributions across docking batches.
    
    Crawls all batch<N> directories under bench_dir and plots selected features
    as separate lines, with batches distinguished by color/label. Useful for
    visualizing the evolution of molecular scores across iterative selection rounds.
    
    Args:
        bench_dir: Root directory containing batch0, batch1, ... subdirectories
        features: Tuple of column names to plot (default: ("CNN_VS", "minimizedAffinity"))
        bins: Number of histogram bins when kde=False (default: 50)
        figsize: Figure size as (width, height) tuple (default: (12, 4))
        style: Matplotlib style name (default: "seaborn-v0_8-colorblind")
        limits: Dict of (min,max) tuples to set x-axis limits per feature
    
    Returns:
        Tuple of (fig, axes) from matplotlib
    
    Raises:
        FileNotFoundError: If bench_dir does not exist or contains no batch<N> directories
    """

    plt.style.use(style)

    bench = os.path.basename(bench_dir.rstrip('/'))
    if not os.path.isdir(bench_dir):
        raise FileNotFoundError(f"Missing bench directory: {bench_dir}")

    # Find batchN directories in numeric order starting at N=0
    batch_dirs = get_batch_dirs(bench_dir)

    # Collect data per feature per batch
    data_per_feature = {feat: [] for feat in features}
    labels_per_feature = {feat: [] for feat in features}

    def add_vals(mdf,label):
        # Append per-feature arrays
        for feat in features:
            vals = mdf[feat].to_numpy(dtype=float, copy=False)
            vals = vals[~np.isnan(vals)]
            if feat == "minimizedAffinity":
                # flip sign so larger = better (if that's your convention)
                vals = -vals
            data_per_feature[feat].append(vals)
            labels_per_feature[feat].append(f'{label} ({vals.mean():.2f})')
            
    with LocalCluster(threads_per_worker=1,n_workers=local_workers,memory_limit="auto") as lcluster, Client(lcluster) as client: 
        # Iterate batches
        for n, bname in batch_dirs:
            bpath = os.path.join(bench_dir, bname)
            docked_path = os.path.join(bpath, "docked.parquet")
            selected_path = os.path.join(bpath, "selected.parquet")

            if not (os.path.isdir(docked_path) and os.path.isdir(selected_path)):
                # Skip incomplete batches
                continue

            docked = dd.read_parquet(docked_path,columns=[*features,'smiles']).compute()
            selected = dd.read_parquet(selected_path,columns=['smiles','score']).compute()

            # Success rate per batch
            success = 100.0 * (len(docked) / max(1, len(selected)))
            print(
                f"{bname}: Docking success rate = {success:.2f}% "
                f"(docked={len(docked)}, selected={len(selected)})"
            )

            # Merge and mask
            mdf = selected.merge(docked, on="smiles", how="inner")

            if n == 0:
                add_vals(mdf[mdf['score'] != 0], "SPRINT")
                add_vals(mdf[mdf['score'] == 0], "Random")
            else:
                add_vals(mdf, f'Batch{n}')


    for kde in (True, False):
        # Plot
        fig, axes = plt.subplots(1, len(features), figsize=figsize, constrained_layout=True)
        if len(features) == 1:
            axes = [axes]

        for ax, feat in zip(axes, features):
            series_list = data_per_feature[feat]
            labels = labels_per_feature[feat]

            # Skip features with no data
            if not series_list or all(len(s) == 0 for s in series_list):
                ax.set_title(f"{feat} (no data)")
                ax.set_xlabel(feat if feat != "minimizedAffinity" else f"-{feat}")
                ax.set_ylabel("Density")
                continue

            # Common x-range for KDE if needed
            if kde:
                all_vals = np.concatenate([s for s in series_list if len(s)])
                xmin, xmax = np.min(all_vals), np.max(all_vals)
                xmin = max(xmin, 0.0)  
                if xmin == xmax:
                    xmin -= 1.0
                    xmax += 1.0
                x = np.linspace(xmin, xmax, 256)

            # Plot each batch as a line
            for vals, lab in zip(series_list, labels):
                lab = f'{lab} ({len(vals)})'
                if len(vals) == 0:
                    continue
                vals = vals[vals > 0]
                if kde:
                    # Safe KDE
                    kde_obj = (
                        gaussian_kde(vals)
                        if (len(vals) >= 2 and np.std(vals) > 0)
                        else None
                    )
                    if kde_obj is None:
                        # Fallback to step hist if KDE degenerate
                        counts, edges = np.histogram(vals, bins=bins, density=True)
                        centers = 0.5 * (edges[:-1] + edges[1:])
                        ax.plot(centers, counts, label=lab)
                    else:
                        ax.plot(x, kde_obj(x), label=lab)
                else:
                    # Simple line histogram (step outline)
                    ax.hist(vals, bins=bins, density=True, histtype="step", label=lab)
            if limits and feat in limits:
                ax.set_xlim(*limits[feat])
            ax.set_ylim(0, None)

            xlab = feat if feat != "minimizedAffinity" else f"-{feat}"

            ax.set_xlabel(xlab)
            ax.set_ylabel("Density")
            ax.legend(loc="best", fontsize="small")
            ax.set_title(feat)

        plt.suptitle(bench)

        out = f'{bench_dir}/{bench}_batch_dists_{"kde" if kde else "hist"}.pdf'
        plt.savefig(out, bbox_inches="tight")
        print(f"Saved: {out}")




def process_smiles_line(line, dbname, fpgen=None, model=None, reduced=False):
    """Process a whitespace-separated line containing a SMILES string and a name.

    Parameters
    ----------
    line : str
        Input line where the first token is expected to be a SMILES string and the last token is the compound name.
    dbname : str
        Identifier for the source database or collection the entry came from.
    fpgen : object, optional
        Fingerprint generator providing a GetFingerprint(mol) method that returns an RDKit fingerprint-like object
        (supports GetNumBits(), ToBinary(), etc.). Required when reduced=False.
    model : object, optional
        Embedding model exposing an embed(torch.Tensor) -> torch.Tensor method; the function will pass a single-batch
        tensor to model.embed under torch.no_grad(). Required when reduced=False.
    reduced : bool
        If True, only return db, name, and canonical SMILES; skip fp and sprint generation.

    Returns
    -------
    tuple or None
        If the SMILES parses successfully, returns a tuple:
            (dbname, name, canonical_smiles) if reduced=True,
            otherwise (dbname, name, canonical_smiles, fingerprint_binary, embedding_array).
        Returns None if RDKit cannot parse the input SMILES.
    """
    vals = line.split()
    if not vals:
        return None
    smi = vals[0]
    name = vals[-1]
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        return None

    smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    if reduced:
        return (dbname, name, smiles)

    fp = fpgen.GetFingerprint(mol)
    features = np.zeros((fp.GetNumBits(),), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, features)
    b = torch.from_numpy(features).unsqueeze(0)
    with torch.no_grad():
        e = model.embed(b)

    return (dbname, name, smiles, fp.ToBinary(), e.squeeze().numpy())


def process_partition(partition_iter, dbname, smiles_only=False):
    """Process a Dask partition of SMILES representations.

    Initializes a Morgan fingerprint generator and optionally loads a pre-trained drug-target
    coembedding model. Iterates through each SMILES line in the partition and applies
    process_smiles_line to generate results.

    Args:
        partition_iter: An iterable of SMILES representation strings to be processed.
        dbname: Database identifier to associate with the processed results.
        smiles_only: If True, generate only db/name/smiles and skip fp/sprint creation.

    Yields:
        Processed results from process_smiles_line for each input SMILES string.
    """
    if smiles_only:
        for line in partition_iter:
            yield process_smiles_line(line, dbname, reduced=True)
        return

    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    model = DrugTargetCoembeddingLightning.load_from_checkpoint(
        "/net/galaxy/home/koes/dkoes/git/panspecies-dti/checkpoints/sprint.ckpt"
    )
    model.eval()

    for line in partition_iter:
        yield process_smiles_line(line, dbname, fpgen=fpgen, model=model)


def process_db(inputfile, dbname, prefix, cluster, blocksize="100KB", repartition=0, smiles_only=False):
    """Process a database of chemical compounds and generate embeddings.
    
    Reads SMILES strings from a text file, canonicalizes them, optionally generates Morgan
    fingerprints and SPRINT embeddings, then writes results to Parquet format.
    
    Args:
        inputfile: Path to the input file containing SMILES strings
        dbname: Database identifier 
        prefix: Output directory path 
        cluster: Cluster object
        blocksize: Block size for reading input file (default: '100KB')
        repartition: Number of partitions for repartitioning (0 = no repartition, default: 0)
        smiles_only: If True, write reduced parquet with only db/name/smiles.
    
    Returns:
        None. Results written to {prefix}/{dbname}.parquet with columns:
        - db: Database identifier
        - name: Compound name
        - smiles: Canonical SMILES string
        - fp: Morgan fingerprint (binary) [optional]
        - sprint: SPRINT embedding (list of float32) [optional]
    """
    with Client(cluster) as client:

        if smiles_only:
            schema = pa.schema(
                [
                    pa.field("db", pa.string()),
                    pa.field("name", pa.string()),
                    pa.field("smiles", pa.string()),
                ]
            )
            columns = ["db", "name", "smiles"]
        else:
            schema = pa.schema(
                [
                    pa.field("db", pa.string()),
                    pa.field("name", pa.string()),
                    pa.field("smiles", pa.string()),
                    pa.field("fp", pa.binary()),  # or pa.large_binary()
                    pa.field("sprint", pa.list_(pa.float32())),  # or pa.float64()
                ]
            )
            columns = ["db", "name", "smiles", "fp", "sprint"]

        bag = db.read_text(inputfile, linedelimiter="\n", blocksize=blocksize)
        if repartition:
            bag = bag.repartition(repartition)
        parsed = bag.map_partitions(partial(process_partition, dbname=dbname, smiles_only=smiles_only)).filter(
            lambda x: x is not None
        )
        df = parsed.to_dataframe(columns=columns)
        df.to_parquet(
            f"{prefix}/{dbname}.parquet",
            engine="pyarrow",
            schema=schema,
            compression="zstd",
        )

def cluster_from_yaml(yaml_path):
    """Initialize a Dask cluster from a YAML configuration file.
    
    Args:
        yaml_path: Path to the YAML configuration file.
        """

    clusters = {
        "slurm": SLURMCluster,
        "pbs": PBSCluster,
        "sge": SGECluster,
        "lsf": LSFCluster,
        "htcondor": HTCondorCluster,
        "oar": OARCluster,
        "moab": MoabCluster,
    }        

    with open(yaml_path, "r") as f:
        config = yaml.safe_load(f)

    dask.config.update(dask.config.config,config)
    jc = config.get("jobqueue")
    if not jc:
        print('Creating LocalCluster')
        return LocalCluster(memory_limit="auto")
    else:
        for n in jc.keys():
            if n in clusters:
                return clusters[n](scheduler_options=jc[n]['scheduler-options'])
            
    raise ValueError(f"Unsupported cluster type in YAML: {yaml_path}")
        
def dask_blocksize_for_file(fname, min_block_kb=100, max_blocks=100_000):
    """
    Compute a blocksize (in bytes) for Dask such that:
      - The file is split into at most `max_blocks` blocks
      - Each block is at least `min_block_kb` kilobytes
    """
    file_size = os.path.getsize(fname)          # in bytes
    min_block = min_block_kb * 1024             # 100 KB → 102400 bytes

    # Smallest blocksize that keeps number of blocks <= max_blocks
    blocksize_for_block_limit = math.ceil(file_size / max_blocks) or 1

    # Enforce both constraints
    blocksize = max(min_block, blocksize_for_block_limit)

    return blocksize

def get_topn(infile, outdir, topn, target_metric="CNN_VS"):
    """Extract top N molecules from all docked batches.
    
    Args:
        infile: Input parquet file path with 'sprint' embeddings and 'name' columns
        outdir: Output directory prefix (results saved to outdir/topn.parquet)
        topn: Number of top molecules to extract by target_metric"""
    
    with LocalCluster(threads_per_worker=1) as lcluster, Client(lcluster) as client:
     
        batch_dirs = get_batch_dirs(outdir)
        dfs = []
        for n, bname in batch_dirs:
            bpath = os.path.join(outdir, bname)
            docked_path = os.path.join(bpath, "docked.parquet")

            if not os.path.isdir(docked_path):
                # Skip incomplete batches
                continue

            docked = dd.read_parquet(docked_path,columns=['mol',target_metric])
            dfs.append(docked)

        if not dfs:
            print("No completed batches found; skipping topn extraction.")
            return

        all_df = dd.concat(dfs, ignore_index=True)
        topdf = all_df.nlargest(topn, target_metric).compute()
        outpath = os.path.join(outdir, "topn.sdf.gz")
        with gzip.open(outpath, "wt", encoding="utf-8") as f:
            for mol in topdf["mol"]:
                f.write(mol)
                if not mol.endswith("$$$$\n"):
                    f.write("$$$$\n")
   
    print(f"Saved top {topn} molecules to: {outpath}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="DeepDock: molecular docking pipeline")
    parser.add_argument( "command", choices=[
            "prepare",
            "initial_batch",
            "next_batch",
            "dock_batch",
            "select_batch",
            "analyze",
            "topn",
            "all", ],
        help="Command to execute. All will skipped already completed steps. All other commands will overwrite existing ouputs.",
    )
    parser.add_argument("--input",  help="Input file (SMILES or parquet)")
    parser.add_argument("--dir", help="Output directory. Must be writeable by all workers. Will be created if not exists.", required=True)
    parser.add_argument("--cluster", help="Path to dask cluster YAML configuration. If not provided will use LocalCluster.")
    parser.add_argument("--num_batches", type=int, default=5, help="Number of batches for 'all' command")
    parser.add_argument("--batch", type=int, default=-1,help="Batch number for dock_batch and select_batch")
    parser.add_argument("--iolimit", type=int, default=50, help="I/O semaphore limit")
    parser.add_argument("--max_workers", type=int, default=500, help="Maximum number of workers")
    parser.add_argument("--target_sequence", help="Target protein sequence")
    parser.add_argument("--sprint_checkpoint", help="Path to SPRINT model checkpoint")
    parser.add_argument("--smiles-only", action="store_true", dest="smiles_only", help="Write reduced parquet with only db,name,smiles and omit fp/sprint generation.")
    parser.add_argument("-N","--batch_size",type=int, default=100_000, help="Number of molecules to select for each batch")
    parser.add_argument("--target_metric", default="CNN_VS", help="Target metric for scoring")
    parser.add_argument("--molweight_cutoff", type=int, default=1200, help="Maximum molecular weight")
    parser.add_argument("--receptor", "-r", help="Receptor file path. Must be accessible by all workers.",required=True)
    parser.add_argument("--autobox_ligand", help="Ligand for autobox generation. Must be accessible by all workers.",required=True)
    parser.add_argument("--gnina_executable", default="gnina", help="Path to gnina executable. Must be accessible by all workers.")
    parser.add_argument("--gnina_timeout", type=int, default=900, help="Timeout per molecule in seconds")
    parser.add_argument("--topn", type=int, default=1000, help="Number of top molecules to extract for 'topn' command")
    parser.add_argument("--unfix_plot_limits", action="store_true", help="Do not fix x-axis limits in analysis plots")
    parser.add_argument("--local_workers", type=int, default=16, help="Number of local workers for LocalCluster")
    parser.add_argument("--local_memory", type=str, default="40GB", help="Memory per worker for LocalCluster")

    args = parser.parse_args()
    os.makedirs(args.dir, exist_ok=True)
    # initialize dask cluster from YAML (tolerance is a virtue)
    dask.config.set({
        "distributed.comm.timeouts.connect": "120s",
        "distributed.comm.timeouts.tcp": "120s",
        "distributed.scheduler.worker-ttl": "120s",
    })
    cluster = cluster_from_yaml(args.cluster) 
    cluster.adapt(minimum=0, maximum=args.max_workers)
    print('Cluster Dashboard:', cluster.dashboard_link,flush=True)
    name = os.path.basename(args.dir.rstrip("/"))

    # need to open lots of files
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (hard, hard))

    if args.command == "prepare" or args.command == "all":
        if os.path.exists(f"{args.dir}/{name}.parquet"):
            print(f"Found existing parquet file: {args.dir}/{name}.parquet; skipping processing.")
        else:
            if args.input is None:
                parser.error("The --input argument is required for 'prepare' and 'all' commands.")
            if args.input.rstrip('/').endswith(".parquet"):
                print("Input appears to be parquet; skipping processing.")
            else:
                process_db(
                    args.input,
                    name,
                    prefix=args.dir,
                    cluster=cluster
                    blocksize=dask_blocksize_for_file(args.input),
                    smiles_only=args.smiles_only,
                )

    if args.input is None:
        args.input = f"{args.dir}/{name}.parquet"

    if not os.path.exists(args.input):
        parser.error(f"Input file does not exist: {args.input}")

    if args.command == "initial_batch" or args.command == "all":
        if args.command == "all" and os.path.exists(f"{args.dir}/batch0/selected.parquet"):
            print(f"Found existing batch0 selected.parquet; skipping initial selection.")
        else:
            # Read target sequence. If a FASTA file path is provided, prefer Biopython SeqIO to read the first sequence.
            seq = ""
            if args.target_sequence:
                if os.path.exists(args.target_sequence):                    
                    with open(args.target_sequence, "r") as fh:
                        rec = next(SeqIO.parse(fh, "fasta"), None)
                        seq = str(rec.seq) if rec is not None else ""
                else:
                    # Treat provided value as the raw sequence string
                    seq = args.target_sequence
            
            nseq = re.sub(r'[^A-Za-z]', '',seq)
            if seq != nseq:
                print("Error: Non-alphabetic characters in target sequence:")
                print(seq);
                sys.exit(-1)

            if args.iolimit < args.max_workers:
                cluster.adapt(minimum=0,maximum=args.iolimit)
            initial_select(
                infile=args.input,
                prefix=args.dir,
                protein_sequence=seq,
                sprint_checkpoint=args.sprint_checkpoint,
                cluster=cluster,
                iolimit=args.iolimit,
                N=args.batch_size,
            )
            cluster.adapt(minimum=0, maximum=args.max_workers)


    if args.command == 'select_batch' or args.command == 'next_batch':
        if args.batch < 0:
            args.batch = next_batch(args.dir)
        if args.batch == 0:
            parser.error("Use initial_batch command for batch 0 selection.")

        select_next_batch(
            infile=args.input,
            outdir=args.dir,
            curbatch=args.batch - 1,
            cluster=cluster,
            score=args.target_metric,
            N=args.batch_size,
            iolimit=args.iolimit,
            local_workers=args.local_workers
        )

    if args.command == "dock_batch" or args.command == 'next_batch' or args.command == 'all':
        if args.batch < 0:
            args.batch = next_batch(args.dir)
        if os.path.exists(f"{args.dir}/batch{args.batch}/selected.parquet"):
            dock_batch(
            outdir=args.dir,
            batch=args.batch,
            gcmd=[
                "-r",
                args.receptor,
                "--autobox_ligand",
                args.autobox_ligand,
            ],
            cluster=cluster,
            )
        else:
            print(f"{args.dir}/batch{args.batch}/selected.parquet does not exist, not docking.")

    if args.command == "all":
        for b in range(1, args.num_batches):

            if os.path.exists(f"{args.dir}/batch{b}/selected.parquet"):
                print(f"Found existing batch {b} selected.parquet; skipping selection.")
            else:
                print(f"Selecting batch {b}...")
                select_next_batch(
                    infile=args.input,
                    outdir=args.dir,
                    curbatch=b - 1,
                    cluster=cluster,
                    score=args.target_metric,
                    N=args.batch_size,
                    iolimit=args.iolimit,
                    local_workers=args.local_workers
                )

            if os.path.exists(f"{args.dir}/batch{b}/docked.parquet"):
                print(f"Found existing batch {b} docked.parquet; skipping docking.")
            else:
                print(f"Docking batch {b}...")
                dock_batch(
                    outdir=args.dir,
                    batch=b,
                    gcmd=[
                    "-r",
                    args.receptor,
                    "--autobox_ligand",
                    args.autobox_ligand,
                ],
                cluster=cluster,
                )

    if args.command == "analyze" or args.command == "all":
        plot_dists_by_batch(
            bench_dir=args.dir,
            features=(
                "CNN_VS",
                "CNNscore",
                "CNNaffinity",
                "minimizedAffinity",
            ),
            bins=50,
            figsize=(16, 4),
            limits=None if args.unfix_plot_limits else {'CNN_VS':(0,10),'CNNscore':(0,1),'CNNaffinity':(2,12),'minimizedAffinity':(0,15)},
            local_workers=args.local_workers
        )

    if args.command == "topn" or args.command == "all":
        get_topn(
            infile=args.input,
            outdir=args.dir,
            topn=args.topn,
            target_metric=args.target_metric,
        )
