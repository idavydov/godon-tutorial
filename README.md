# Godon tutorial

## Installing godon

If you are using GNU/Linux, go to the
[website](https://bitbucket.org/Davydov/godon/) and grab a precompiled
binary from the Downloads section.

```
$ wget https://bitbucket.org/Davydov/godon/downloads/godon-master-linux-gnu-x86_64 -O godon
```

Now make this binary executable:

```
$ chmod +x godon
```

You could always check find documentation by calling `godon --help`:

```
$ ./godon --help
usage: godon [<flags>] <command> [<args> ...]

codon models optmizer and sampler

Flags:
  -h, --help                     Show context-sensitive help (also try --help-long and --help-man).
  -v, --version                  Show application version.
      --gcode=1                  NCBI genetic code id, standard by default
      --fg-branch=-1             foreground branch number
      --max-branch-length=100    maximum branch length
  -n, --no-branch-length         don't optimize branch lengths
      --codon-frequency="F3X4"   codon frequecny (F0 or F3X4)
      --codon-frequency-file=CODON-FREQUENCY-FILE  
                                 codon frequencies file (overrides --codon-frequency)
      --ncat-site-rate=1         number of categories for the site rate variation (no variation by default)
      --ncat-codon-rate=1        number of categories for the codon rate variation (no variation by default)
      --proportional             use three rates and three proportions instead of gamma distribution
      --ncat-beta=4              number of the categories for the beta distribution (models M7&M8)
  -s, --start=START              read start position from the trajectory or JSON file
      --randomize-start          use uniformly distributed random starting point; by default random starting point is distributed around realistic parameter values
      --iter=10000               number of iterations
      --report=1                 report every N iterations
  -m, --method="lbfgsb"          optimization method to use (lbfgsb: limited-memory Broyden–Fletcher–Goldfarb–Shanno with bounding constraints, simplex: downhill simplex, annealing:
                                 simullated annealing, mh: Metropolis-Hastings, n_lbfgs: LBFGS from nlopt, n_simplex: downhill simplex from nlopt, n_cobyla: COBYLA from nlopt,
                                 n_bobyqa: BOBYQA from nlopt, n_sqp: SQP from nlopt, n_mlsl: MLSL from nlopt (BOBYQA local optimizer), none: just compute likelihood, no optimization)
      --final                    perform final extra computations, i.e. NEB and BEB site posterior (default on, use --no-final to disable)
      --neb                      perform naive empirical bayes of positive selection (default on, use --no-neb to disable)
      --beb                      perform bayes empirical bayes of positive selection (default on, use --no-beb to disable)
      --codon-rates              perform NEB analysis of codon rates
      --site-rates               perform NEB analysis of site rates
      --codon-omega              perform NEB analysis of codon omega
      --report-acceptance=200    report acceptance rate every N iterations
      --adaptive                 use adaptive MCMC or sumulated annealing
      --skip-adaptive=-1         number of iterations to skip for adaptive mcmc (5% by default)
      --maximum-adaptive=-1      stop adapting after iteration (20% by default)
      --aggregate=none           state aggregation mode: observed (all positions, keep observed states), observed_new (new implementation of observed), fixed (absolutely conserved
                                 positions, keep observed), random (like observed, but non-aggregated states are shuffled between the positions)
  -p, --procs=PROCS              number of threads to use
  -S, --seed=-1                  random generator seed, default time based
      --cpu-profile=CPU-PROFILE  write cpu profile to file
  -o, --out=OUT                  write log to a file
  -t, --trajectory=TRAJECTORY    write optimization trajectory to a file
  -l, --log-level=notice         set loglevel (critical, error, warning, notice, info, debug)
  -j, --json=JSON                write json output to a file

Commands:
  help [<command>...]
    Show help.

  optimize* [<flags>] <model> <alignment> <tree>
    Run model optimization or sampling

  test [<flags>] <model> <alignment> <tree>
    Run test for positive selection
```

## Gene sequence and a tree

We will use Drosohpila gene
[VPS13](http://www.flybase.org/reports/FBgn0033194.html) and its
orthologs from 12 Drosophila species. Vacuolar protein sorting 13
(Vps13) encodes a protein located to endosomal membrane that is
involved in protein homeostasis.

They are already available in this repository. When creating your own
files for analysis make sure to have identical sequence identifiers
used in the tree and in the sequence alignment.

There are a couple of caveats when analyzing your own sequence. First,
make sure to provide a good-quality codon-based sequence
alignment. Misalignment could lead to false positive selection events
identified. One way to create a high-quality sequence alignment is to
use [GUIDANCE2](https://dx.doi.org/10.1093/nar/gkv318) in codon mode
and mask (replace with `X`) residues with scores below certain
threshold. In
[Selectome](https://selectome.unil.ch/cgi-bin/methods.cgi) the
threshold is set to 0.93.

Second, mask the selenocysteine-encoding codons, as they are not
supported by majority of codon models. Finally, make sure to create a
high-quality phylogenetic tree. E.g., using software such are
[PhyML](http://www.atgc-montpellier.fr/phyml/) or
[RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html),
or software based on Bayesian approach. If branch support values are
low, consider removing some sequences or trying multiple topologies,
and choosing the most conservative estimates of positive selection.

Also, it is important to keep in mind that branch-length estimates
provided by software such as PhyML is expressed in the number of
nucleotide substitutions per unit of pseudotime. In codon models the
branch lengths are describing the number of codon substitutions. So it
is important to allow Godon or other codon model software to estimate
branch lengths.

## The simplest model

We first use Godon with the classical codon model, also known as
[Goldman Yang
1994](https://doi.org/10.1093/oxfordjournals.molbev.a040153) or M0.

By using this model we assume there is a single dN/dS ratio over all
the sites of the alignment and all the branches of the tree. This
assumption is almost never satisfied. However, this model allows to
detect an extremely strong positive selection affecting majority of
the positions over long period of time. This model has a low
computational cost, therefore it is frequently used to estimate branch
lengths for the codon alignment. These branch lengths could be then
used without further optimization for more complex codon models, to
improve computational performance without compromising statistical
properties too much (see supplementary materials of [this
paper](https://doi.org/10.1093/molbev/msz048)).

Now let's run the model:

```
$ ./godon M0 --out-tree M0_tree.nwk EMGT00050000000025.Drosophila.001.fst EMGT00050000000025.Drosophila.001.nwk
Tree is rooted. Will unroot.
Starting lnL=-55216.574
Maximum likelihood: -48511.34704457653
omega=0.08243971396522087
kappa=1.7393328157675232
Running time: 44.50772981s
```

By default Godon will try to use all the available CPUs, you can
control it with `--procs N` argument. We use `--out-tree` to save the
tree for further use.

