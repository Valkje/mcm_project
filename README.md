# Monte Carlo Methods Project

This repo contains: 

* Two report files, an HTML and a PDF version. Personally I think the HTML version is the prettiest, but you can use either to get a coherent overview of the code and the results it produces.
* A `project.jl` file, which just contains the code and is easiest to run.
* Some toml files that will help you install the necessary Julia packages (further explained below).

## Getting started

To run the code, first install Julia.

Visit https://julialang.org/downloads/ and select the installer for your system (I used 1.8.3, but probably 1.9.0 will work just as well). If you're on Linux, run from a location of your choosing:

```bash
wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.3-linux-x86_64.tar.gz
tar zxvf julia-1.8.3-linux-x86_64.tar.gz
```

On Linux, that's all you have to do. You can choose to make the `julia` executable more accessible by, e.g., opening `~/.bashrc` and adding the line:

```bash
export PATH="$PATH:/path/to/<Julia directory>/bin"
```

Make sure to start a new terminal or run `source ~/.bashrc` afterwards.

Consequently, you want to install the relevant packages.

Move one level above the current directory (the current directory is the one that contains this file) and activate Julia:

```bash
cd /path/to/mcm_project/..
julia
```

Inside Julia, change to the Pkg mode by typing `]`. In there, you want to activate the project toml file and instantiate the project. This will install all necessary packages and precompile them. Make sure that `mcm_project` reflects the actual folder name that you cloned from github.

```julia
]
activate mcm_project
instantiate
```

## Running the code

If you have Visual Studio Code installed, you could open the `project.ipynb` file, select Julia as your notebook kernel, and then hopefully be able to run everything cell by cell. Slightly simpler would be to run the more barebones `project.jl`, which contains some additional code to create an `images` directory and save the generated images as png files. To use the last option, run:

```bash
cd mcm_project
julia --project=. project.jl
```