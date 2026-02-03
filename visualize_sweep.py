#!/usr/bin/env python3
"""Visualize sweeping turbulence frames as an animated GIF."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import glob
from PIL import Image
import io
from tqdm import tqdm

def read_bin(filename):
    """Read binary file with (nx, ny) header, k_current, followed by float64 data."""
    with open(filename, 'rb') as f:
        nx, ny = np.fromfile(f, dtype=np.int32, count=2)
        k_current = np.fromfile(f, dtype=np.float64, count=1)[0]
        data = np.fromfile(f, dtype=np.float64).reshape((ny, nx))
    return data, k_current

def main():
    frames = sorted(glob.glob('frame_*.bin'))
    print(f"Found {len(frames)} frames")
    
    if not frames:
        print("No frames found!")
        return
    
    # Skip first frame
    frames = frames[1:]
    
    # Read all frames to get global min/max for consistent colormap
    all_data = []
    all_k = []
    for f in tqdm(frames, desc='Reading frames'):
        data, k_current = read_bin(f)
        all_data.append(data)  # Already |F|² from Fortran
        all_k.append(k_current)
    
    # Get range for log scale (need positive values)
    all_positive = [d[d > 0] for d in all_data]
    vmin = min(p.min() for p in all_positive if len(p) > 0)
    vmax = max(d.max() for d in all_data)
    print(f"Global |F|² range: [{vmin:.2e}, {vmax:.2e}]")
    print(f"k range: [{min(all_k):.1f}, {max(all_k):.1f}]")
    norm = LogNorm(vmin=vmin, vmax=vmax)
    
    # Create images
    images = []
    time_axis = np.arange(len(frames))
    for i, (f, data) in enumerate(tqdm(zip(frames, all_data), total=len(frames), desc='Rendering')):
        fig, (ax_main, ax_k) = plt.subplots(2, 1, figsize=(8, 10), 
                                             gridspec_kw={'height_ratios': [4, 1]})
        
        # Main turbulence plot (|F|²) with log scale
        im = ax_main.imshow(data, cmap='inferno', norm=norm, origin='lower')
        ax_main.set_title(f'Frame {i+1}/{len(frames)}  k={all_k[i]:.1f}')
        ax_main.axis('off')
        plt.colorbar(im, ax=ax_main, fraction=0.046, pad=0.04)
        
        # k vs time plot
        ax_k.plot(time_axis, all_k, 'b-', linewidth=1)
        ax_k.axvline(x=i, color='r', linewidth=2)
        ax_k.set_xlim(0, len(frames)-1)
        ax_k.set_ylim(min(all_k) - 2, max(all_k) + 2)
        ax_k.set_xlabel('Frame')
        ax_k.set_ylabel('k')
        ax_k.set_title('Driven wavenumber')
        
        plt.tight_layout()
        
        # Save to buffer
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        images.append(Image.open(buf).copy())
        plt.close()
    
    # Save as GIF
    print("Creating animated GIF...")
    images[0].save('sweep_animation.gif', save_all=True, append_images=images[1:],
                   duration=100, loop=0)
    print("Saved sweep_animation.gif")
    
    # Also save first, middle, last frames as static images
    for idx, name in [(0, 'sweep_first.png'), (len(all_data)//2, 'sweep_middle.png'), 
                       (-1, 'sweep_last.png')]:
        actual_idx = idx if idx >= 0 else len(frames) + idx
        fig, (ax_main, ax_k) = plt.subplots(2, 1, figsize=(10, 12),
                                             gridspec_kw={'height_ratios': [4, 1]})
        im = ax_main.imshow(all_data[idx], cmap='inferno', norm=norm, origin='lower')
        ax_main.set_title(f'Frame {actual_idx+1}  k={all_k[idx]:.1f}')
        ax_main.axis('off')
        plt.colorbar(im, ax=ax_main, fraction=0.046, pad=0.04)
        
        ax_k.plot(time_axis, all_k, 'b-', linewidth=1)
        ax_k.axvline(x=actual_idx, color='r', linewidth=2)
        ax_k.set_xlim(0, len(frames)-1)
        ax_k.set_ylim(min(all_k) - 2, max(all_k) + 2)
        ax_k.set_xlabel('Frame')
        ax_k.set_ylabel('k')
        ax_k.set_title('Driven wavenumber')
        
        plt.tight_layout()
        plt.savefig(name, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved {name}")

if __name__ == '__main__':
    main()
