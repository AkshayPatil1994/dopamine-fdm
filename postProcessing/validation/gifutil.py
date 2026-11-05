# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Frame-to-GIF helper shared by the validation animations (ffmpeg palette pipeline, fixed width and fps)."""
import subprocess


def frames_to_gif(pattern, out, fps=10, width=720):
    """Encode frames matching an ffmpeg pattern such as 'frames/f_%04d.png' into an optimised looping GIF."""
    vf = f"fps={fps},scale={width}:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128:stats_mode=diff[p];" \
         "[b][p]paletteuse=dither=bayer:bayer_scale=5"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", pattern, "-vf", vf,
                    "-loop", "0", out], check=True)
