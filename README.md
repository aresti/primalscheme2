# primalscheme2

### Why a new repo?
primalscheme2 is a greenfield codebase and uses a different algorithm to the original [primalscheme](https://github.com/aresti/primalscheme). We wanted to make the project public even in its unfinished state, since it includes useful new features.

**The test suite has not yet been ported so please consider this a beta, and subject to breaking changes.**

#### New features

 - Entirely new algorithm that yeilds superior results, particularly on 'difficult' genomes.
 - Significantly improved interaction checker.
 - 'Repair' mode, to check an existing scheme against a new N-masked input, replacing primers where required.
 - Limited ambiguity code support (alpha).
 
 #### Roadmap
 - 'Panel' inputs
 - Further ambiguity code support
 - 'Jackhammer' strategy/mode

### Usage notes
The original [primalscheme](https://github.com/aresti/primalscheme) accepted multi-fasta inputs, and handled alignment internally. For primalscheme2, you should now create your n-masked reference upfront using [Chris Kent's multi-fasta-mask pipeline](https://github.com/ChrisgKent/multi_fasta_mask), using the output as your primalscheme input.

That being said, you can still use a single genome input to primalscheme without any preliminary step.

### Installation
primalscheme2 is not yet available via. Conda or pypy. For now, please install from source.

#### Requirements

 - Python >= 3.9
 - [Poetry](https://python-poetry.org/)
 - [Chris Kent's multi-fasta-mask pipeline](https://github.com/ChrisgKent/multi_fasta_mask)

#### Install from source

 1. Install [Poetry](https://python-poetry.org/docs/#installation)
 2. `git clone https://github.com/aresti/primalscheme2.git`
 3. `cd primalscheme2`
 4. `poetry install` (by default, Poetry will handle your venv for you. If you'd prefer to manage that yourself, read the Poetry docs for further details)
 5. Activate your venv via. `poetry shell`
 6. Verify your installation: `primalscheme --help`

### Usage

Having first activated your venv via. `poetry shell`:

#### Basic usage
`primalscheme my_input.fa`

#### Options

 - `--output DIRECTORY`
 - `--amplicon-size-min INTEGER RANGE [100<=x<=2000]`
 - `--amplicon-size-max INTEGER RANGE [100<=x<=2000]`
 -  `--min-overlap INTEGER RANGE [x>=0]`
 - `--force / --no-force` (overwrite existing output dir)
 - `--debug / --no-debug`
 - `--prefix <str> Prefix name for your outputs.  [default: scheme]`
 - `--repair FILE`
 - `--repair-interactions / --no-repair-interactions`

