#!/usr/bin/env python3
"""Check that our tiles are the input Granite was trained on, using the model itself.

1. Configuration. Compares tile size, patch size, channel count, band names and
   order, and the per-band normalisation statistics in our code with the model's
   own config.yaml from the Hugging Face snapshot.
2. Positional embedding. Compares the checkpoint's positional embedding with the
   fixed 2D sin/cos embedding the architecture defines for a 42x42 input.
3. Input statistics. For a sample of tiles: per-channel mean and s.d. after
   normalisation (should be of order 0 and 1 if units and scaling are right), and
   how much of each tile is flagged, filled or missing SST.
4. Reconstruction. Granite is a masked autoencoder: shown a tile with 75% of its
   patches hidden, its decoder predicts the hidden pixels. If our tiles are in the
   form it was trained on, it reconstructs them far better than a constant guess.
   The test is repeated for alternative ways of filling unusable pixels, scoring
   only clear pixels, so the conventions can be compared on equal terms.
5. Embedding sensitivity. How much the embeddings move under those alternatives.

Writes data/results/granite_input_check.json and prints a summary.

Usage (CPU is enough; no internet once the model is cached):
    python scripts/check_granite_inputs.py --n-tiles 300
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import granite, tiles  # noqa: E402
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, from_repo, read_csv  # noqa: E402

OUR_BANDS = [f"{band}_reflectance" for band in tiles.OLCI_BANDS] + ["SST"]
OUR_MEAN = np.array([*tiles.OLCI_MEANS, tiles.SST_MEAN_K], dtype=np.float64)
OUR_STD = np.array([*tiles.OLCI_STDS, tiles.SST_STD_K], dtype=np.float64)


def main() -> None:
    args = parse_args()
    torch.manual_seed(args.seed)
    report: dict = {}
    snapshot = snapshot_dir(args.model_name, args.local_files_only)

    report["config"] = check_config(snapshot)
    checkpoint = torch.load(snapshot / "checkpoint.pt", map_location="cpu")
    state = granite._state_dict(checkpoint)
    report["positional_embedding"] = check_pos_embed(state)

    paths = sample_tiles(args)
    batch = [load(path) for path in paths]
    report["input_statistics"] = input_statistics(batch)

    mae = build_mae(state)
    variants = {name: torch.from_numpy(np.stack([fill(t, name) for t in batch])) for name in VARIANTS}
    scored = torch.from_numpy(np.stack([t["scored"] for t in batch]))
    report["reconstruction"] = reconstruction(mae, variants, scored, args)

    encoder, _ = granite.build_encoder(args.model_name, local_files_only=True)
    report["embedding_sensitivity"] = embedding_sensitivity(encoder.eval(), variants)

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(report, indent=2))
    summarise(report)
    print(f"\nWrote {args.out_json}")


# ----------------------------------------------------------------------------
def snapshot_dir(model_name: str, local_files_only: bool) -> Path:
    from huggingface_hub import snapshot_download

    return Path(snapshot_download(model_name, local_files_only=local_files_only))


def check_config(snapshot: Path) -> dict:
    cfg = json.loads((snapshot / "config.yaml").read_text())["pretrained_cfg"]  # the file is JSON despite its name
    mean, std = np.array(cfg["mean"]), np.array(cfg["std"])
    out = {
        "model_img_size": cfg["img_size"],
        "our_tile_size": tiles.TILE_SIZE,
        "img_size_matches": cfg["img_size"][-2:] == [tiles.TILE_SIZE, tiles.TILE_SIZE],
        "patch_size_matches": cfg["patch_size"] == [1, 2, 2],
        "in_chans_matches": cfg["in_chans"] == len(OUR_BANDS),
        "bands_and_order_match": cfg["bands"] == OUR_BANDS,
        "mean_max_relative_difference": float(np.max(np.abs(mean - OUR_MEAN) / np.abs(mean))),
        "std_max_relative_difference": float(np.max(np.abs(std - OUR_STD) / np.abs(std))),
    }
    out["all_match"] = bool(
        out["img_size_matches"] and out["patch_size_matches"] and out["in_chans_matches"]
        and out["bands_and_order_match"] and out["mean_max_relative_difference"] < 1e-6
        and out["std_max_relative_difference"] < 1e-6
    )  # fmt: skip
    return out


def check_pos_embed(state: dict) -> dict:
    from terratorch.models.backbones.prithvi_mae import get_3d_sincos_pos_embed

    key = next(k for k in state if k.endswith("pos_embed") and "decoder" not in k)
    stored = state[key].float().numpy()[0]
    expected = get_3d_sincos_pos_embed(stored.shape[-1], (1, tiles.TILE_SIZE // 2, tiles.TILE_SIZE // 2), add_cls_token=True)
    return {"key": key, "shape": list(stored.shape), "max_abs_difference": float(np.abs(stored - expected).max())}


# ----------------------------------------------------------------------------
def sample_tiles(args: argparse.Namespace) -> list[Path]:
    rows = read_csv(args.manifest_csv)
    paths = sorted({row["tile_path"] for row in rows if from_repo(row["tile_path"]).exists()})
    random.Random(args.seed).shuffle(paths)
    return [from_repo(path) for path in paths[: args.n_tiles]]


def load(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as tile:
        data = tile["data"].astype(np.float32)
        clear = tile["valid_mask"].astype(bool)
    sst_ok = np.isfinite(data[-1])
    # Pixels a fair comparison can score: clear water, with SST present.
    return {"data": data, "clear": clear, "sst_ok": sst_ok, "scored": np.broadcast_to(clear & sst_ok, data.shape).copy()}


def input_statistics(batch: list[dict]) -> dict:
    data = np.stack([t["data"] for t in batch])  # (n, 17, 42, 42)
    clear = np.stack([t["clear"] for t in batch])
    per_channel = {}
    for c, name in enumerate(OUR_BANDS):
        values = data[:, c][clear & np.isfinite(data[:, c])]
        per_channel[name] = {"mean": float(values.mean()), "std": float(values.std())}
    return {
        "n_tiles": len(batch),
        "clear_fraction_mean": float(clear.mean()),
        "sst_missing_fraction_mean": float(np.mean([~t["sst_ok"] for t in batch])),
        "tiles_with_any_missing_sst": int(sum((~t["sst_ok"]).any() for t in batch)),
        "olci_exact_zero_fraction": float((data[:, :16] == 0).mean()),  # fill pixels are stored as 0
        "normalised_clear_pixels": per_channel,
    }


# ----------------------------------------------------------------------------
# Ways to fill pixels that are not clear water. All in normalised units.
VARIANTS = {
    "as_embedded": "flagged pixels keep their values; missing SST -> 0 (channel mean). What embed_tiles.py did.",
    "flagged_to_mean": "flagged pixels and missing SST -> 0 (channel mean) in every channel",
    "flagged_to_raw_zero": "flagged pixels and missing SST -> physical value 0, then normalised "
    "(TerraTorch's no_data_replace=0 convention, as in IBM's fine-tuning config)",
}


def fill(tile: dict, variant: str) -> np.ndarray:
    data = tile["data"].copy()
    unusable = ~tile["clear"]
    if variant == "as_embedded":
        data = np.nan_to_num(data, nan=0.0)
    elif variant == "flagged_to_mean":
        data[:, unusable] = 0.0
        data = np.nan_to_num(data, nan=0.0)
    elif variant == "flagged_to_raw_zero":
        raw_zero = (-OUR_MEAN / OUR_STD).astype(np.float32)
        data[:, unusable] = raw_zero[:, None]
        sst = data[-1]
        sst[~np.isfinite(sst)] = raw_zero[-1]
    return data.astype(np.float32)


def build_mae(state: dict):
    from terratorch.models.backbones.prithvi_mae import PrithviMAE

    mae = PrithviMAE(
        img_size=tiles.TILE_SIZE, patch_size=(1, 2, 2), num_frames=1, in_chans=17, embed_dim=granite.EMBED_DIM,
        depth=12, num_heads=8, decoder_embed_dim=256, decoder_depth=4, decoder_num_heads=8, mlp_ratio=4,
        norm_pix_loss=False, mask_ratio=0.75,
    )  # fmt: skip
    target = mae.state_dict()
    mapped = {}
    for key, value in state.items():
        key = key.replace("_timm_module.", "")
        if key.startswith(("encoder.", "decoder.")):
            new = key
        elif "decoder" in key or "_dec" in key or key == "mask_token":
            new = "decoder." + key
        else:
            new = "encoder." + key
        if new in target and target[new].shape == value.shape:
            mapped[new] = value
    result = mae.load_state_dict(mapped, strict=False)
    if result.missing_keys:
        raise RuntimeError(f"MAE weights missing after loading: {result.missing_keys[:10]}")
    return mae.eval()


def reconstruction(mae, variants: dict, scored: torch.Tensor, args: argparse.Namespace) -> dict:
    """RMSE (normalised units) on hidden patches, clear pixels only, per input convention."""
    results = {}
    reference = variants["as_embedded"].unsqueeze(2)  # what the clear pixels truly are
    weight = scored.unsqueeze(2).float()
    for name, inputs in variants.items():
        squared, constant, local, count, spread_in, spread_out, n_spread = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0
        per_channel_sq, per_channel_n, per_channel_bias = torch.zeros(17), torch.zeros(17), torch.zeros(17)
        for repeat in range(args.repeats):
            torch.manual_seed(args.seed + repeat)  # same hidden patches for every variant
            with torch.no_grad():
                _, pred, mask = mae(inputs.unsqueeze(2), mask_ratio=0.75)
            image, hidden = as_images(mae, pred, mask)
            w = weight * hidden
            error = (image - reference) ** 2 * w
            squared += float(error.sum())
            constant += float((reference**2 * w).sum())  # guessing the training mean (0) everywhere
            # Guessing each hidden pixel as the mean of the tile's visible clear pixels, per channel.
            seen = weight * (1 - hidden)
            visible_mean = (reference * seen).sum(dim=(2, 3, 4), keepdim=True) / seen.sum(dim=(2, 3, 4), keepdim=True).clamp(min=1)
            local += float(((visible_mean - reference) ** 2 * w).sum())
            count += float(w.sum())
            per_channel_sq += error.sum(dim=(0, 2, 3, 4))
            per_channel_bias += ((image - reference) * w).sum(dim=(0, 2, 3, 4))
            per_channel_n += w.sum(dim=(0, 2, 3, 4))
            # Spread across the 17 channels at each hidden clear pixel: does the decoder reproduce spectral shape?
            pixels = (w[:, 0, 0] > 0)
            spread_in += float(reference[:, :, 0].permute(0, 2, 3, 1)[pixels].std(1).sum())
            spread_out += float(image[:, :, 0].permute(0, 2, 3, 1)[pixels].std(1).sum())
            n_spread += int(pixels.sum())
        results[name] = {
            "description": VARIANTS[name],
            "rmse_hidden_clear_pixels": (squared / count) ** 0.5,
            "rmse_constant_guess": (constant / count) ** 0.5,
            "rmse_visible_mean_guess": (local / count) ** 0.5,
            "channel_spread_input": spread_in / max(n_spread, 1),
            "channel_spread_decoder_output": spread_out / max(n_spread, 1),
            "bias_per_channel": {b: float(per_channel_bias[i] / per_channel_n[i]) for i, b in enumerate(OUR_BANDS)},
            "rmse_per_channel": {b: float((per_channel_sq[i] / per_channel_n[i]) ** 0.5) for i, b in enumerate(OUR_BANDS)},
        }
    return results


def as_images(mae, pred: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    """Prediction (n, 17, 1, 42, 42) and hidden-pixel mask (n, 1, 1, 42, 42), whatever form TerraTorch returns.

    Older TerraTorch versions return patch sequences; newer ones return images.
    """
    size = (tiles.TILE_SIZE, tiles.TILE_SIZE)
    if pred.ndim == 3:  # (n, patches, values per patch)
        if mask.ndim == 2:
            mask = mask.unsqueeze(-1).repeat(1, 1, pred.shape[-1])
        pred, mask = mae.unpatchify(pred, image_size=size), mae.unpatchify(mask, image_size=size)[:, :1]
    image = pred if pred.ndim == 5 else pred.unsqueeze(2)
    hidden = mask.float()
    while hidden.ndim < 5:
        hidden = hidden.unsqueeze(1)
    return image, hidden[:, :1]


def embedding_sensitivity(encoder, variants: dict) -> dict:
    with torch.no_grad():
        embedded = {name: granite.encode(encoder, inputs).numpy() for name, inputs in variants.items()}
    base = embedded["as_embedded"]
    out = {}
    for name, other in embedded.items():
        cosine = (base * other).sum(1) / (np.linalg.norm(base, axis=1) * np.linalg.norm(other, axis=1))
        out[name] = {"cosine_to_as_embedded_median": float(np.median(cosine)), "cosine_min": float(cosine.min())}
    return out


# ----------------------------------------------------------------------------
def summarise(report: dict) -> None:
    c = report["config"]
    print("1. Configuration vs the model's config.yaml")
    print(f"   tile {c['our_tile_size']}x{c['our_tile_size']} vs model {c['model_img_size']}: {'OK' if c['img_size_matches'] else 'MISMATCH'}")
    print(f"   patch size: {'OK' if c['patch_size_matches'] else 'MISMATCH'} | channels: {'OK' if c['in_chans_matches'] else 'MISMATCH'}"
          f" | band names and order: {'OK' if c['bands_and_order_match'] else 'MISMATCH'}")  # fmt: skip
    print(f"   normalisation: mean max rel. diff {c['mean_max_relative_difference']:.1e}, s.d. {c['std_max_relative_difference']:.1e}")
    p = report["positional_embedding"]
    print(f"2. Positional embedding {p['shape']}: max |checkpoint - sin/cos| = {p['max_abs_difference']:.2e}")
    s = report["input_statistics"]
    print(f"3. Inputs ({s['n_tiles']} tiles): clear fraction {s['clear_fraction_mean']:.3f}; SST missing on "
          f"{s['sst_missing_fraction_mean']:.3f} of pixels ({s['tiles_with_any_missing_sst']} tiles with some); "
          f"OLCI fill pixels {s['olci_exact_zero_fraction']:.4f}")  # fmt: skip
    print("   normalised clear pixels, mean / s.d. per channel:")
    for name, stats in s["normalised_clear_pixels"].items():
        print(f"     {name:18s} {stats['mean']:7.2f} / {stats['std']:5.2f}")
    print("4. Reconstruction of hidden patches (clear pixels only; normalised units; lower is better)")
    for name, r in report["reconstruction"].items():
        print(f"   {name:20s} model {r['rmse_hidden_clear_pixels']:.3f} | constant guess {r['rmse_constant_guess']:.3f}"
              f" | tile's own visible-pixel mean {r['rmse_visible_mean_guess']:.3f}")  # fmt: skip
    r = report["reconstruction"]["as_embedded"]
    print(f"   spread across the 17 channels at a pixel: input {r['channel_spread_input']:.3f}, "
          f"decoder output {r['channel_spread_decoder_output']:.3f}")  # fmt: skip
    bias = report["reconstruction"]["as_embedded"]["bias_per_channel"]
    print("   mean error (prediction - truth) per channel, as embedded: "
          + ", ".join(f"{name.split('_')[0]} {value:+.2f}" for name, value in bias.items()))  # fmt: skip
    print("5. Embedding change vs as embedded (cosine similarity, median / min)")
    for name, r in report["embedding_sensitivity"].items():
        print(f"   {name:20s} {r['cosine_to_as_embedded_median']:.4f} / {r['cosine_min']:.4f}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-csv", type=Path, default=MANIFEST_DIR / "biogoship_tile_manifest.csv")
    parser.add_argument("--out-json", type=Path, default=REPO_ROOT / "data" / "results" / "granite_input_check.json")
    parser.add_argument("--n-tiles", type=int, default=300)
    parser.add_argument("--repeats", type=int, default=3, help="Random masks per tile.")
    parser.add_argument("--model-name", default=granite.MODEL_NAME)
    parser.add_argument("--local-files-only", action="store_true")
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


if __name__ == "__main__":
    main()
