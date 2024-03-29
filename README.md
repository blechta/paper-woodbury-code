# Supporting code for paper Efficient solution of parameter identification problems with $H^1$ regularization #

<a href="https://doi.org/10.5281/zenodo.6855783">
<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.6521320.svg" alt="DOI">
</a>
<a href="https://doi.org/10.48550/arXiv.2209.02815">
<img src="https://zenodo.org/badge/DOI/10.48550/arXiv.2209.02815.svg" alt="DOI">
</a>


## Requirements ##

* MATLAB R2020b with Symbolic Math Toolbox,
* [HSL_MI20 2.0.0](https://www.hsl.rl.ac.uk/packages/hsl_mi20.html),
* [Gmsh 4.5.6](https://gmsh.info),
* [matlab2tikz f299888](https://github.com/matlab2tikz/matlab2tikz/tree/f299888fc381a6976009678a7dc00e6fe6872bd2).


## Steps to reproduce the paper results ##

1. Start Matlab and `cd` to the repository root
   (where this README resides).
2. Run
   ```
   app_dc.checkerboard_example.short_benchmark;
   ```
   This step produces some `checkerboard-*.{fig,xdmf,bin,mat}`
   files in the current directory.
   This step can take around 900 minutes and needs 15 GB of RAM.

   To get the full results of the paper, use
   ```
   app_dc.checkerboard_example.full_benchmark;
   ```
   instead. This can take around 75 hours (depending on machine)
   and needs around 250 GB of RAM.
3. Run
   ```
   app_dc.checkerboard_example.postprocess;
   app_dc.checkerboard_example.postprocess2;
   app_dc.checkerboard_example.postprocess3;
   app_dc.checkerboard_example.postprocess4;
   app_dc.checkerboard_example.plot_configs;
   ```
   to produce additional
   `checkerboard-*.{pdf,tex}` files in
   the current directory.


## Steps to produce 3d plots ##

The 3d plots are produced by the Python script `plot-3d.py`
from `.xdmf` and `.bin` files produced by the aforementioned
steps. The script uses `pyvista` and a `VTK` build with `XDMF3`
support, which is not featured on PyPI. One can circumvent this
by using Docker.
```
docker build -t pyvista dockerfiles/
docker run --rm -v "$PWD":/tmp pyvista python3 -u plot-3d.py
```
This produces `checkerboard-3dfw.svg` and
`checkerboard-3dfw-beta.svg`.


## Authors ##

* Jan Blechta <jan.blechta@math.tu-chemnitz.de>
* Mathias Scheunert <mathias.scheunert@mailserver.tu-freiberg.de>


## License ##

MIT License
```
Copyright (C) 2018-2022 Jan Blechta and Mathias Scheunert

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
