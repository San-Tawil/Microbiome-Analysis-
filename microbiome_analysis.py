#!/usr/bin/env python3
"""
Complete Microbiome Analysis Pipeline
======================================
Taxonomic profiling, alpha diversity, beta diversity, PERMANOVA & PERMDISP
for MetaPhlAn relative-abundance tables.

Dependencies: numpy, pandas, matplotlib, seaborn, scipy, scikit-learn
(NO scikit-bio required — all diversity metrics implemented directly)
"""

# ---------------------------------------------------------------------------
# 0. Imports & global config
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use('Agg')

import os
import sys
import warnings
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import seaborn as sns

from scipy.spatial.distance import pdist, squareform, braycurtis, jaccard
from scipy.stats import kruskal, mannwhitneyu
from sklearn.manifold import MDS

warnings.filterwarnings('ignore')
sns.set_style("whitegrid")

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

TABLE_FILES = {
    'phylum':  os.path.join(SCRIPT_DIR, '8_rel-phyla-table.tsv'),
    'family':  os.path.join(SCRIPT_DIR, '8_rel-family-table.tsv'),
    'genus':   os.path.join(SCRIPT_DIR, '8_rel-genus-table.tsv'),
    'species': os.path.join(SCRIPT_DIR, '8_rel-species-table.tsv'),
}
METADATA_FILE = os.path.join(SCRIPT_DIR, 'metadata.tsv')
RESULTS_DIR   = os.path.join(SCRIPT_DIR, 'results')

GROUP_COL = 'State'
TOP_N_VALUES = [10, 15]

# ---------------------------------------------------------------------------
# Consistent colour palette (fixed dictionary for major phyla, then auto)
# ---------------------------------------------------------------------------
FIXED_COLORS = {
    # Phylum-level (matching reference image style)
    'Bacteroidota':       '#6BAED6',   # light blue
    'Chloroflexi':        '#3182BD',   # medium blue
    'Firmicutes':         '#A1D99B',   # light green
    'Actinobacteriota':   '#31A354',   # green
    'Synergistota':       '#FC9272',   # salmon/pink
    'Proteobacteria':     '#DE2D26',   # red
    'Verrucomicrobiota':  '#FEC44F',   # gold/yellow
    'Patescibacteria':    '#E31A1C',   # darker red
    'Desulfobacterota':   '#9E9AC8',   # lavender
    'Cyanobacteria':      '#3F007D',   # dark purple
    'Bdellovibrionota':   '#FFFFB3',   # pale yellow
    'TA06':               '#B15928',   # brown-red
    'Armatimonadota':     '#80B1D3',   # pastel blue
    'Spirochaetota':      '#1F78B4',   # blue
    'Thermotogota':       '#33A02C',   # green
    'Euryarchaeota':      '#006D2C',   # dark green
    'WPS-2':              '#FB6A4A',   # light red
    'Halanaerobiaeota':   '#C51B7D',   # magenta
    'Fusobacteriota':     '#FD8D3C',   # orange
    'Dependentiae':       '#E6550D',   # dark orange
    'Other':              '#D9D9D9',   # light grey
    'Unassigned':         '#BDBDBD',   # grey
}

_TAB20  = plt.cm.tab20(np.linspace(0, 1, 20))
_TAB20B = plt.cm.tab20b(np.linspace(0, 1, 20))
_SET3   = plt.cm.Set3(np.linspace(0, 1, 12))
_ALL_COLORS = list(_TAB20) + list(_TAB20B) + list(_SET3)

def _get_color(name, idx):
    """Return colour from fixed dict or fall back to auto palette."""
    if name in FIXED_COLORS:
        return FIXED_COLORS[name]
    return _ALL_COLORS[idx % len(_ALL_COLORS)]


# ===================================================================
# 1. DATA LOADING & VALIDATION
# ===================================================================

def load_metadata(path):
    meta = pd.read_csv(path, sep='\t', comment=None)
    first_col = meta.columns[0]
    if first_col.startswith('#'):
        meta = meta.rename(columns={first_col: first_col.lstrip('#')})
    id_col = meta.columns[0]
    meta = meta.set_index(id_col)
    meta.index = meta.index.astype(str).str.strip()
    log.info(f"Metadata: {meta.shape[0]} samples, columns={list(meta.columns)}")
    return meta


def _extract_taxon_label(otu_id, level):
    parts = otu_id.split(';')
    prefix_map = {'phylum': 'p__', 'family': 'f__', 'genus': 'g__', 'species': 's__'}
    target_prefix = prefix_map.get(level, None)

    if target_prefix:
        for p in parts:
            p = p.strip()
            if p.startswith(target_prefix):
                name = p[len(target_prefix):]
                if name and name != '__' and name.strip():
                    return name

    for p in reversed(parts):
        p = p.strip()
        if any(p.startswith(pfx) for pfx in ['d__', 'k__']):
            continue
        for pfx in ['p__', 'c__', 'o__', 'f__', 'g__', 's__']:
            if p.startswith(pfx):
                p = p[len(pfx):]
                break
        if p and p != '__' and p.strip():
            return p
    return 'Unassigned'


def load_abundance_table(path, level):
    df = pd.read_csv(path, sep='\t', skiprows=1)
    otu_col = df.columns[0]
    sample_cols = [c for c in df.columns if c != otu_col]
    rename = {c: c.strip() for c in sample_cols}
    df = df.rename(columns=rename)
    sample_cols = [c.strip() for c in sample_cols]

    labels = df[otu_col].apply(lambda x: _extract_taxon_label(str(x), level))
    data = df[sample_cols].T
    data.columns = labels
    data = data.T.groupby(level=0).sum().T
    log.info(f"[{level}] {data.shape[0]} samples × {data.shape[1]} taxa")
    return data


def validate(tables, meta):
    for level, df in tables.items():
        common = set(meta.index) & set(df.index)
        log.info(f"[{level}] {len(common)} samples matched")


# ===================================================================
# 2. TAXONOMIC PROFILING — Grouped horizontal stacked bar (reference style)
# ===================================================================

def taxonomic_barplots(tables, meta, out_dir):
    """
    Grouped horizontal stacked barplots matching reference style:
    - Y-axis = Group (Patient / Control)
    - X-axis = Relative Abundance (%) 0–100, gridlines every 5%
    - Taxa sorted by mean abundance
    - Unassigned / Other excluded from vis, remaining renormalized to 100%
    - Legend on the right
    """
    bar_dir = os.path.join(out_dir, 'taxa_barplots')
    os.makedirs(bar_dir, exist_ok=True)

    EXCLUDE_VIS = {'Unassigned', 'Other', 'unassigned', 'other', '__'}

    for level, df in tables.items():
        common = sorted(set(df.index) & set(meta.index))
        df_c = df.loc[common].copy()
        groups = meta.loc[common, GROUP_COL]

        # Drop Unassigned / Other columns for visualization
        vis_cols = [c for c in df_c.columns if c not in EXCLUDE_VIS]
        df_vis = df_c[vis_cols].copy()

        # Renormalize each sample row to sum to 1 (then ×100 for %)
        row_sums = df_vis.sum(axis=1)
        row_sums = row_sums.replace(0, 1)  # avoid division by zero
        df_vis = df_vis.div(row_sums, axis=0)

        for top_n in TOP_N_VALUES:
            # Mean abundance per taxon (after renorm)
            means = df_vis.mean().sort_values(ascending=False)
            top_taxa = means.head(top_n).index.tolist()

            plot_df = df_vis[top_taxa].copy()

            # Convert to %
            plot_df_pct = plot_df * 100
            plot_df_pct['Group'] = groups.values
            group_means = plot_df_pct.groupby('Group').mean()

            all_taxa = top_taxa

            # Rename groups for display
            display_names = {'FMF': 'Patient', 'Control': 'Control'}
            group_means.index = [display_names.get(g, g) for g in group_means.index]
            group_order = ['Patient', 'Control']
            group_means = group_means.reindex([g for g in group_order if g in group_means.index])

            # Colours
            colors = [_get_color(t, i) for i, t in enumerate(all_taxa)]

            # --- Grouped barplot ---
            fig, ax = plt.subplots(figsize=(12, 3.5 + 0.6 * len(group_means)))
            y_pos = np.arange(len(group_means))
            bar_height = 0.6
            left = np.zeros(len(group_means))

            for idx, taxon in enumerate(all_taxa):
                vals = group_means[taxon].values
                ax.barh(y_pos, vals, left=left, height=bar_height,
                        color=colors[idx], label=taxon,
                        edgecolor='white', linewidth=0.4)
                left += vals

            ax.set_yticks(y_pos)
            ax.set_yticklabels(group_means.index, fontsize=12, weight='bold')
            ax.set_xlabel('Relative Abundance (%)', fontsize=11)
            ax.set_title(f'{level.capitalize()}', fontsize=14, weight='bold')
            ax.set_xlim(0, 100)
            ax.set_xticks(np.arange(0, 101, 5))
            ax.tick_params(axis='x', labelsize=9)
            ax.grid(axis='x', linestyle='-', alpha=0.4, color='grey')
            ax.set_axisbelow(True)
            ax.invert_yaxis()
            ax.set_ylabel('Group', fontsize=11)

            ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                      fontsize=8, title=level.capitalize(), title_fontsize=9,
                      ncol=1, frameon=True, edgecolor='grey')

            plt.tight_layout()
            fname = f'{level}_top{top_n}_barplot.png'
            fig.savefig(os.path.join(bar_dir, fname), dpi=300, bbox_inches='tight')
            plt.close(fig)

            # --- Per-sample barplot ---
            fig2, ax2 = plt.subplots(figsize=(14, max(8, len(common) * 0.35)))
            order = groups.sort_values().index
            df_plot_sample = plot_df_pct.drop(columns='Group').loc[order]
            y_pos2 = np.arange(len(df_plot_sample))
            left2 = np.zeros(len(df_plot_sample))

            for idx, taxon in enumerate(all_taxa):
                vals = df_plot_sample[taxon].values
                ax2.barh(y_pos2, vals, left=left2, height=0.8,
                         color=colors[idx], label=taxon,
                         edgecolor='white', linewidth=0.3)
                left2 += vals

            ax2.set_yticks(y_pos2)
            ax2.set_yticklabels(df_plot_sample.index, fontsize=7)
            ax2.set_xlabel('Relative Abundance (%)', fontsize=11)
            ax2.set_title(f'{level.capitalize()} — Per Sample (Top {top_n})', fontsize=13, weight='bold')
            ax2.set_xlim(0, 100)
            ax2.set_xticks(np.arange(0, 101, 5))
            ax2.tick_params(axis='x', labelsize=9)
            ax2.grid(axis='x', linestyle='-', alpha=0.4, color='grey')
            ax2.set_axisbelow(True)
            ax2.invert_yaxis()

            sample_groups = groups.loc[order]
            gcm = {'FMF': '#E74C3C', 'Control': '#3498DB'}
            for i, sample in enumerate(df_plot_sample.index):
                g = sample_groups.loc[sample]
                ax2.barh(i, 0.6, left=-0.8, height=0.8,
                         color=gcm.get(g, 'grey'), clip_on=False)

            g_patches = [mpatches.Patch(color=c, label=g) for g, c in gcm.items()]
            handles, labels_leg = ax2.get_legend_handles_labels()
            ax2.legend(handles + g_patches, labels_leg + list(gcm.keys()),
                       bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7,
                       title='Taxa / Group')

            plt.tight_layout()
            fname2 = f'{level}_top{top_n}_per_sample_barplot.png'
            fig2.savefig(os.path.join(bar_dir, fname2), dpi=300, bbox_inches='tight')
            plt.close(fig2)

            # CSV exports
            group_means.to_csv(os.path.join(bar_dir, f'{level}_top{top_n}_group_means.csv'))
            plot_df_pct.drop(columns='Group').to_csv(
                os.path.join(bar_dir, f'{level}_top{top_n}_per_sample_data.csv'))
            log.info(f"  Saved {fname}, {fname2}")


# ===================================================================
# 3. ALPHA DIVERSITY
# ===================================================================

def _shannon(p):
    p = p[p > 0]
    return -np.sum(p * np.log(p))

def _simpson(p):
    return np.sum(p ** 2)

def _inverse_simpson(p):
    s = _simpson(p)
    return 1.0 / s if s > 0 else 0.0


def compute_alpha(tables, meta, out_dir):
    plot_dir = os.path.join(out_dir, 'alpha_diversity', 'plots')
    val_dir  = os.path.join(out_dir, 'alpha_diversity', 'values')
    stat_dir = os.path.join(out_dir, 'alpha_diversity', 'stats')
    for d in [plot_dir, val_dir, stat_dir]:
        os.makedirs(d, exist_ok=True)

    metrics = {'Shannon': _shannon, 'Simpson': _simpson, 'Inverse_Simpson': _inverse_simpson}
    all_stats = []

    for level, df in tables.items():
        common = sorted(set(df.index) & set(meta.index))
        df_c = df.loc[common]
        groups = meta.loc[common, GROUP_COL]

        records = []
        for sample in common:
            row = df_c.loc[sample].values.astype(float)
            rec = {'Sample': sample, 'Group': groups.loc[sample]}
            for mname, mfunc in metrics.items():
                rec[mname] = mfunc(row)
            records.append(rec)

        div_df = pd.DataFrame(records)
        div_df.to_csv(os.path.join(val_dir, f'{level}_alpha_diversity.csv'), index=False)

        # Boxplots
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        for ax, mname in zip(axes, metrics):
            sns.boxplot(data=div_df, x='Group', y=mname, ax=ax, palette='Set2',
                        width=0.5, fliersize=3)
            sns.stripplot(data=div_df, x='Group', y=mname, ax=ax,
                          color='black', size=4, alpha=0.6, jitter=True)
            ax.set_title(mname, fontsize=12, weight='bold')
            ax.set_xlabel('')
        fig.suptitle(f'Alpha Diversity — {level.capitalize()}', fontsize=14, weight='bold', y=1.02)
        plt.tight_layout()
        fig.savefig(os.path.join(plot_dir, f'{level}_alpha_boxplots.png'),
                    dpi=300, bbox_inches='tight')
        plt.close(fig)

        # Stats
        unique_groups = div_df['Group'].unique()
        for mname in metrics:
            gv = [div_df.loc[div_df['Group'] == g, mname].values for g in unique_groups]
            if len(gv) >= 2 and all(len(v) > 0 for v in gv):
                stat_kw, p_kw = kruskal(*gv)
            else:
                stat_kw, p_kw = np.nan, np.nan
            all_stats.append({'Level': level, 'Metric': mname, 'Test': 'Kruskal-Wallis',
                              'Statistic': stat_kw, 'p-value': p_kw,
                              'Groups': ' vs '.join(unique_groups)})
            if len(unique_groups) == 2:
                u_stat, p_mw = mannwhitneyu(gv[0], gv[1], alternative='two-sided')
                all_stats.append({'Level': level, 'Metric': mname, 'Test': 'Mann-Whitney U',
                                  'Statistic': u_stat, 'p-value': p_mw,
                                  'Groups': ' vs '.join(unique_groups)})

        log.info(f"  Alpha diversity done for {level}")

    pd.DataFrame(all_stats).to_csv(os.path.join(stat_dir, 'alpha_diversity_stats.csv'), index=False)


# ===================================================================
# 4. BETA DIVERSITY (no scikit-bio — pure scipy/sklearn)
# ===================================================================

def _clr_transform(df, pseudocount=1e-6):
    data = df.values + pseudocount
    log_data = np.log(data)
    gm = log_data.mean(axis=1, keepdims=True)
    return pd.DataFrame(log_data - gm, index=df.index, columns=df.columns)


def _compute_dm(df, metric='braycurtis'):
    """Compute square distance matrix using scipy."""
    ids = list(df.index)
    n = len(ids)
    vals = df.values
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if metric == 'braycurtis':
                d = braycurtis(vals[i], vals[j])
            elif metric == 'jaccard':
                a = (vals[i] > 0).astype(float)
                b = (vals[j] > 0).astype(float)
                d = jaccard(a, b)
            elif metric == 'euclidean':
                d = np.sqrt(np.sum((vals[i] - vals[j]) ** 2))
            else:
                d = 0
            mat[i, j] = d
            mat[j, i] = d
    return mat, ids


def _pcoa(dm_matrix):
    """Classical MDS (PCoA) from a distance matrix."""
    n = dm_matrix.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (dm_matrix ** 2) @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Keep positive eigenvalues
    pos = eigvals > 0
    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]

    coords = eigvecs_pos * np.sqrt(eigvals_pos)
    total_var = eigvals_pos.sum()
    explained = eigvals_pos / total_var * 100 if total_var > 0 else eigvals_pos * 0
    return coords, explained


def _draw_confidence_ellipse(ax, x, y, color, n_std=1.96, **kwargs):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    if np.any(np.isnan(cov)) or np.any(np.isinf(cov)):
        return
    try:
        vals, vecs = np.linalg.eigh(cov)
    except np.linalg.LinAlgError:
        return
    order = vals.argsort()[::-1]
    vals = np.abs(vals[order])
    vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)
    ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                  width=width, height=height, angle=angle,
                  edgecolor=color, facecolor=color, alpha=0.15,
                  linewidth=1.5, linestyle='--', **kwargs)
    ax.add_patch(ell)


def _ordination_plot(coords, ids, groups_series, method, metric, out_path,
                     explained=None, col_prefix='PC'):
    group_colors = {'FMF': '#E74C3C', 'Control': '#3498DB'}
    fig, ax = plt.subplots(figsize=(8, 7))

    for g in groups_series.unique():
        mask = np.array([groups_series.loc[i] == g for i in ids])
        x = coords[mask, 0]
        y = coords[mask, 1]
        c = group_colors.get(g, 'grey')
        ax.scatter(x, y, c=c, label=g, s=60, edgecolors='black', linewidth=0.5, alpha=0.85, zorder=3)
        _draw_confidence_ellipse(ax, x, y, color=c)

    xl = f'{col_prefix}1'
    yl = f'{col_prefix}2'
    if explained is not None and len(explained) >= 2:
        xl += f' ({explained[0]:.1f}%)'
        yl += f' ({explained[1]:.1f}%)'
    ax.set_xlabel(xl, fontsize=11)
    ax.set_ylabel(yl, fontsize=11)
    ax.set_title(f'{method} — {metric}', fontsize=13, weight='bold')
    ax.legend(fontsize=10)
    ax.axhline(0, color='grey', lw=0.5, ls='--', zorder=1)
    ax.axvline(0, color='grey', lw=0.5, ls='--', zorder=1)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def _distance_boxplot(dm_mat, ids, groups_series, metric, out_path):
    records = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            gi = groups_series.loc[ids[i]]
            gj = groups_series.loc[ids[j]]
            dtype = 'Within' if gi == gj else 'Between'
            records.append({'Distance': dm_mat[i, j], 'Type': dtype})
    bp_df = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.boxplot(data=bp_df, x='Type', y='Distance', palette='Set2', ax=ax, width=0.5)
    sns.stripplot(data=bp_df, x='Type', y='Distance', color='black',
                  size=3, alpha=0.4, jitter=True, ax=ax)
    ax.set_title(f'Distance Distribution — {metric}', fontsize=12, weight='bold')
    ax.set_ylabel(metric)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    bp_df.to_csv(out_path.replace('.png', '_data.csv'), index=False)


def _permanova(dm_mat, groups_arr, n_perms=999):
    """PERMANOVA (one-way) — Anderson 2001."""
    unique_grps = np.unique(groups_arr)
    n = len(groups_arr)

    def _pseudo_f(labels):
        # SS_T: total sum of squared distances / N
        ss_t = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                ss_t += dm_mat[i, j] ** 2
        ss_t /= n

        # SS_W: within-group sum of squared distances / n_g
        ss_w = 0.0
        for g in unique_grps:
            idx_g = np.where(labels == g)[0]
            ng = len(idx_g)
            if ng > 1:
                for ii in range(ng):
                    for jj in range(ii + 1, ng):
                        ss_w += dm_mat[idx_g[ii], idx_g[jj]] ** 2 / ng

        ss_a = ss_t - ss_w
        k = len(unique_grps)
        df_a = k - 1
        df_w = n - k

        if df_w <= 0 or ss_w == 0:
            return 0.0
        return (ss_a / df_a) / (ss_w / df_w)

    observed_f = _pseudo_f(groups_arr)

    count = 0
    for _ in range(n_perms):
        perm = np.random.permutation(groups_arr)
        if _pseudo_f(perm) >= observed_f:
            count += 1

    p_value = (count + 1) / (n_perms + 1)
    return observed_f, p_value


def _permdisp(dm_mat, groups_arr, n_perms=999):
    """PERMDISP: test homogeneity of dispersions."""
    unique_grps = np.unique(groups_arr)
    n = len(groups_arr)

    # Compute centroid distance for each sample (distance to group centroid in PCoA space)
    coords, _ = _pcoa(dm_mat)
    coords2 = coords[:, :min(2, coords.shape[1])]

    dispersions = np.zeros(n)
    for g in unique_grps:
        mask = groups_arr == g
        centroid = coords2[mask].mean(axis=0)
        dists = np.sqrt(np.sum((coords2[mask] - centroid) ** 2, axis=1))
        dispersions[mask] = dists

    # F-statistic comparing dispersions between groups
    group_disps = [dispersions[groups_arr == g] for g in unique_grps]
    if len(group_disps) >= 2 and all(len(v) > 1 for v in group_disps):
        # One-way F-test on dispersions
        overall_mean = dispersions.mean()
        ss_b = sum(len(gd) * (gd.mean() - overall_mean) ** 2 for gd in group_disps)
        ss_w = sum(np.sum((gd - gd.mean()) ** 2) for gd in group_disps)
        k = len(unique_grps)
        df_b = k - 1
        df_w = n - k
        f_obs = (ss_b / df_b) / (ss_w / df_w) if ss_w > 0 and df_w > 0 else 0

        count = 0
        for _ in range(n_perms):
            perm = np.random.permutation(groups_arr)
            gd_p = [dispersions[perm == g] for g in unique_grps]
            om = dispersions.mean()
            ss_b_p = sum(len(gd) * (gd.mean() - om) ** 2 for gd in gd_p)
            ss_w_p = sum(np.sum((gd - gd.mean()) ** 2) for gd in gd_p)
            f_p = (ss_b_p / df_b) / (ss_w_p / df_w) if ss_w_p > 0 else 0
            if f_p >= f_obs:
                count += 1
        p_value = (count + 1) / (n_perms + 1)
    else:
        f_obs, p_value = np.nan, np.nan

    return f_obs, p_value


def compute_beta(tables, meta, out_dir):
    dm_dir   = os.path.join(out_dir, 'beta_diversity', 'distance_matrices')
    pcoa_dir = os.path.join(out_dir, 'beta_diversity', 'pcoa')
    nmds_dir = os.path.join(out_dir, 'beta_diversity', 'nmds')
    box_dir  = os.path.join(out_dir, 'beta_diversity', 'boxplots')
    perm_dir = os.path.join(out_dir, 'beta_diversity', 'permanova')
    disp_dir = os.path.join(out_dir, 'beta_diversity', 'permdisp')
    for d in [dm_dir, pcoa_dir, nmds_dir, box_dir, perm_dir, disp_dir]:
        os.makedirs(d, exist_ok=True)

    df = tables['genus']
    common = sorted(set(df.index) & set(meta.index))
    df_c = df.loc[common]
    groups = meta.loc[common, GROUP_COL]
    groups_arr = groups.values

    df_relabund = df_c.copy()
    df_clr = _clr_transform(df_c)

    perm_results = []

    configs = {
        'Bray-Curtis': ('braycurtis', df_relabund),
        'Jaccard':     ('jaccard',    df_relabund),
        'Aitchison':   ('euclidean',  df_clr),
    }

    for metric_name, (metric_key, data_df) in configs.items():
        log.info(f"  Computing {metric_name} ...")
        dm_mat, ids = _compute_dm(data_df, metric=metric_key)
        metric_safe = metric_name.replace('-', '_')

        # Save distance matrix
        dm_df = pd.DataFrame(dm_mat, index=ids, columns=ids)
        dm_df.to_csv(os.path.join(dm_dir, f'{metric_safe}_distance_matrix.csv'))

        # PCoA
        log.info(f"  PCoA for {metric_name} ...")
        coords_pcoa, explained = _pcoa(dm_mat)
        for top_n in TOP_N_VALUES:
            _ordination_plot(coords_pcoa, ids, groups, f'PCoA (top {top_n})',
                             metric_name,
                             os.path.join(pcoa_dir, f'PCoA_{metric_safe}_top{top_n}.png'),
                             explained=explained, col_prefix='PC')

        pd.DataFrame(coords_pcoa[:, :min(5, coords_pcoa.shape[1])],
                     index=ids,
                     columns=[f'PC{i+1}' for i in range(min(5, coords_pcoa.shape[1]))]).to_csv(
            os.path.join(pcoa_dir, f'PCoA_{metric_safe}_coords.csv'))

        # NMDS
        log.info(f"  NMDS for {metric_name} ...")
        mds = MDS(n_components=2, dissimilarity='precomputed',
                  random_state=RANDOM_SEED, max_iter=300, n_init=10,
                  normalized_stress='auto')
        coords_nmds = mds.fit_transform(dm_mat)
        for top_n in TOP_N_VALUES:
            _ordination_plot(coords_nmds, ids, groups, f'NMDS (top {top_n})',
                             metric_name,
                             os.path.join(nmds_dir, f'NMDS_{metric_safe}_top{top_n}.png'),
                             col_prefix='NMDS')

        pd.DataFrame(coords_nmds, index=ids, columns=['NMDS1', 'NMDS2']).to_csv(
            os.path.join(nmds_dir, f'NMDS_{metric_safe}_coords.csv'))

        # Distance boxplot
        _distance_boxplot(dm_mat, ids, groups, metric_name,
                          os.path.join(box_dir, f'{metric_safe}_boxplot.png'))

        # PERMANOVA
        log.info(f"  PERMANOVA for {metric_name} ...")
        try:
            f_stat, p_val = _permanova(dm_mat, groups_arr, n_perms=999)
            entry = {'Metric': metric_name, 'Test': 'PERMANOVA',
                     'F_statistic': f_stat, 'p_value': p_val, 'permutations': 999}
            perm_results.append(entry)
            pd.DataFrame([entry]).to_csv(
                os.path.join(perm_dir, f'PERMANOVA_{metric_safe}.csv'), index=False)
            log.info(f"    F={f_stat:.4f}, p={p_val:.4f}")
        except Exception as e:
            log.error(f"  PERMANOVA failed: {e}")

        # PERMDISP
        log.info(f"  PERMDISP for {metric_name} ...")
        try:
            f_disp, p_disp = _permdisp(dm_mat, groups_arr, n_perms=999)
            entry = {'Metric': metric_name, 'Test': 'PERMDISP',
                     'F_statistic': f_disp, 'p_value': p_disp, 'permutations': 999}
            perm_results.append(entry)
            pd.DataFrame([entry]).to_csv(
                os.path.join(disp_dir, f'PERMDISP_{metric_safe}.csv'), index=False)
            log.info(f"    F={f_disp:.4f}, p={p_disp:.4f}")
        except Exception as e:
            log.error(f"  PERMDISP failed: {e}")

    if perm_results:
        pd.DataFrame(perm_results).to_csv(
            os.path.join(out_dir, 'beta_diversity', 'beta_stats_summary.csv'), index=False)
    log.info("  Beta diversity complete")


# ===================================================================
# MAIN
# ===================================================================

def main():
    log.info("=" * 60)
    log.info("MICROBIOME ANALYSIS PIPELINE")
    log.info("=" * 60)

    meta = load_metadata(METADATA_FILE)
    tables = {}
    for level, path in TABLE_FILES.items():
        tables[level] = load_abundance_table(path, level)
    validate(tables, meta)

    os.makedirs(RESULTS_DIR, exist_ok=True)

    log.info("STEP 1: Taxonomic profiling")
    taxonomic_barplots(tables, meta, RESULTS_DIR)

    log.info("STEP 2: Alpha diversity")
    compute_alpha(tables, meta, RESULTS_DIR)

    log.info("STEP 3: Beta diversity & statistics")
    compute_beta(tables, meta, RESULTS_DIR)

    log.info("=" * 60)
    log.info("COMPLETE — results in: %s", RESULTS_DIR)
    log.info("=" * 60)


if __name__ == '__main__':
    import traceback
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
