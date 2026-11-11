//! CRUSH — Compartment Refinement for Ultraprecise Stratification of Hi-C.
//!
//! Standalone Rust implementation. Reads .hic natively, holds the pipeline in
//! memory, and parallelises over chromosomes with rayon. Command-line flags
//! mirror the bash implementation so existing invocations port over unchanged.

mod bed;
mod diagnostic;
mod eigen;
mod gc;
mod hic;
mod score;
mod stats;
mod v1;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process;
use std::sync::Mutex;

use rayon::prelude::*;

use bed::{BedSet, ChromSizes};
use hic::{HicFile, ValueKind};
use score::{ChromScoreInput, Moments};

const VERSION: &str = env!("CARGO_PKG_VERSION");

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ShifterMode {
    /// Subtract the local rolling mean only (bash v1.0 behaviour).
    Highpass,
    /// Subtract (rolling mean - coarse score), preserving regional structure.
    Coarse,
}

/// Which of the walk's per-resolution tracks to keep alongside the merge.
///
/// CRUSH scores every resolution on the ladder on its way down, and the final
/// round re-scores all of them seeded from the finest result. Those tracks are
/// discarded once merged, but they are the natural thing to look at when asking
/// where in the walk a call was decided, so they can be written out on request.
#[derive(Debug, Clone, PartialEq, Eq)]
enum KeepRes {
    /// Write only the merged output (default).
    None,
    /// Write every resolution the ladder walked.
    All,
    /// Write just these resolutions. Validated against the ladder at run time,
    /// since the ladder is not known until the .hic is opened.
    Only(Vec<u64>),
}

#[derive(Clone)]
struct Config {
    hic: PathBuf,
    sizes: PathBuf,
    initial_a: PathBuf,
    initial_b: PathBuf,
    res: u64,
    coarsest_res: u64,
    eigen_res: u64,
    cpu: usize,
    outpre: String,
    outdir: PathBuf,
    norm: String,
    threshold: f64,
    window_kb: f64,
    min_distance: u64,
    max_distance: u64,
    qthresh: f64,
    pcalculation: bool,
    shifter_mode: ShifterMode,
    eigenfile: Option<PathBuf>,
    eigen_init: bool,
    /// Resolutions >= this are mean-centred per chromosome. `None` disables
    /// it, which is the default: bash v1 initialised endZ=0 and then rewrote it
    /// to the finest resolution before use, so the centring always ran — an
    /// "optional" parameter whose default meant "always on".
    end_z: Option<u64>,
    /// Drop bins CRUSH's own row marginal reports as having no coverage.
    vc_filter: bool,
    /// Convert outputs to bigWig (needs UCSC bedGraphToBigWig on PATH).
    bigwig: bool,
    /// Keep the bedGraph alongside the bigWig.
    keep_bedgraph: bool,
    /// Prefix bigWig chromosome names with "chr" when they lack it.
    chr_prefix: bool,
    /// Split re-derived seeds at the per-chromosome median rather than at zero.
    seed_split_median: bool,
    /// Fraction of bins used as seeds, taken from the extremes of the score
    /// distribution. 1.0 seeds every bin (historical); 0.5 seeds the top and
    /// bottom quartiles and leaves the ambiguous middle out.
    seed_frac: f64,
    /// Resolution at which seeds are maintained, independent of the pass being
    /// scored. 0 uses the pass's own resolution.
    seed_res: u64,
    /// Interpolate the coarse track between bin centres instead of scattering
    /// it as a step function (removes sign flips at coarse bin boundaries).
    coarse_interp: bool,
    /// Score the round's target resolution first, then re-run every resolution
    /// seeded from THAT result (v1's reprocess_resolutions_with_shifter).
    reverse_pass: bool,
    /// Restrict seeds to the positions the initial annotation supplied (gene
    /// bodies for A, low-GC bins for B), relabelling them from the previous
    /// round rather than seeding from every bin.
    seed_anchor: bool,
    /// How much of the previous round's seed state to carry forward, 0..1.
    /// 0 replaces the seeds outright each round (current); higher values blend,
    /// slowing how fast the walk can move and limiting how much a small data
    /// difference compounds across rounds.
    seed_damping: f64,
    /// Clamp on the FINAL output scores, applied after the mean-|score|=100
    /// rescale. A handful of bins land orders of magnitude off the bulk and
    /// there is nothing to read in their exact value, but IGV's autoscale is
    /// driven by the extremes, so one such bin flattens the whole track to a
    /// line. 0 disables.
    clamp: f64,
    /// Which per-resolution tracks from the final round to write out.
    keep_res: KeepRes,
    verbose: bool,
    /// Run the line-for-line port of the published bash v1 pipeline (v1.rs).
    v1: bool,
    /// Walk from both annotation and eigenvector seeds, keep the stronger per
    /// chromosome (--seed-select both). False when one seed is requested.
    seed_both: bool,
    /// Flip chromosomes whose track anti-correlates with gene/GC density.
    orient: bool,
    /// --v1 only: bash v1's Final Shifter and 3-point smoothing passes.
    reshift: bool,
    smoothing: bool,
    max_matrix_bytes: usize,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            hic: PathBuf::new(),
            sizes: PathBuf::new(),
            initial_a: PathBuf::new(),
            initial_b: PathBuf::new(),
            res: 1_000,
            coarsest_res: 1_000_000,
            eigen_res: 100_000,
            cpu: 24,
            outpre: String::new(),
            outdir: PathBuf::from("CRUSHtracks"),
            norm: "NONE".into(),
            threshold: 10.0,
            window_kb: 5.0,
            min_distance: 0,
            max_distance: 0,
            qthresh: 0.05,
            pcalculation: false,
            shifter_mode: ShifterMode::Coarse,
            eigenfile: None,
            eigen_init: true,
            end_z: Some(1_000),
            vc_filter: true,
            bigwig: true,
            keep_bedgraph: true,
            chr_prefix: true,
            seed_split_median: false,
            seed_frac: 1.0,
            seed_res: 500,
            seed_anchor: false,
            reverse_pass: true,
            coarse_interp: true,
            seed_damping: 0.0,
            clamp: 500.0,
            keep_res: KeepRes::None,
            verbose: false,
            v1: false,
            seed_both: true,
            orient: true,
            reshift: false,
            smoothing: false,
            max_matrix_bytes: 8 * 1024 * 1024 * 1024,
        }
    }
}

fn usage() -> String {
    format!(
        "crush-rs {VERSION} — kilobase-scale A/B compartment calling from Hi-C

USAGE:
  crush-rs -i FILE.hic -g SIZES -a INITIAL_A.bed -b INITIAL_B.fa [OPTIONS]

REQUIRED:
  -i, --hic FILE          Input .hic file
  -g, --genomesize FILE   Chromosome sizes (2 columns: name, length)
  -a, --initialA FILE     BED seeding the A compartment (gene annotations)
  -b, --initialB FILE     BED seeding the B compartment, or a genome FASTA to
                          derive low-GC B seeds from

COMMON OPTIONS:
  -r, --res BP            Finest resolution to analyse            [default: 1000]
  -m, --maxres BP         Coarsest resolution to walk from     [default: 1000000]
  -c, --cpu N             Worker threads                            [default: 24]
  -o, --outpre PREFIX     Output filename prefix
  -O, --outdir DIR        Output directory                 [default: CRUSHtracks]
  -N, --norm TYPE         NONE | VC | VC_SQRT | KR | SCALE        [default: NONE]
  -v, --verbose           Verbose progress reporting
  -h, --help              This message
  -H, --help-advanced     Every option, including the algorithm internals

OUTPUT:
  {{outpre}}mergedCrush_{{res}}.bedgraph / .bw   the compartment track
  Scores are rescaled so mean |score| = 100, then clamped to ±500 so a few
  extreme bins cannot flatten IGV's autoscale (--clamp).

Run with -H for p-values, per-resolution output, seeding and shifter controls."
    )
}

fn usage_advanced() -> String {
    format!(
        "crush-rs {VERSION} — kilobase-scale A/B compartment calling from Hi-C
ADVANCED OPTIONS. For the short list, run -h.

USAGE:
  crush-rs -i FILE.hic -g SIZES -a INITIAL_A.bed -b INITIAL_B.fa [OPTIONS]

REQUIRED:
  -i, --hic FILE          Input .hic file
  -g, --genomesize FILE   Chromosome sizes (2 columns: name, length)
  -a, --initialA FILE     BED seeding the A compartment (gene annotations)
  -b, --initialB FILE     BED seeding the B compartment (low-GC bins), or a
                          genome FASTA to derive them from

RESOLUTION LADDER:
  -r, --res BP            Finest resolution to analyse            [default: 1000]
  -m, --maxres BP         Coarsest resolution to walk from     [default: 1000000]
  -Z, --eigenres BP       Resolution for eigenvector init        [default: 100000]

OUTPUT:
  -o, --outpre PREFIX     Output filename prefix
  -O, --outdir DIR        Output directory                 [default: CRUSHtracks]
  -B, --bigwig 0|1        Write bigWig output (needs UCSC bedGraphToBigWig
                          on PATH); 0 writes bedGraph only           [default: 1]
      --keep-bedgraph 0|1 Keep the bedGraph next to the bigWig       [default: 1]
      --chr-prefix 0|1    Prefix bigWig chromosome names with 'chr' when they
                          lack it, so the file loads against IGV hg19/hg38.
                          Affects the bigWig only, not the bedGraph  [default: 1]
      --clamp V           Clamp final scores to ±V, applied after the
                          mean-|score|=100 rescale. A handful of bins sit orders
                          of magnitude off the bulk; their exact value carries
                          nothing, but IGV autoscales to them and flattens the
                          rest of the track to a line. 0 disables the clamp and
                          writes the raw rescaled values           [default: 500]
      --keepres SPEC      Also write the final round's per-resolution tracks,
                          which are otherwise discarded at the merge. 'all'
                          keeps every resolution on the ladder, or give a
                          comma-separated list of specific ones:
                            --keepres all
                            --keepres 100000,10000,5000
                          Spaces after the commas are accepted. Each track is
                          rescaled and clamped independently and lands in
                          {{outpre}}Crush_{{res}}.bedgraph. Listing a resolution
                          the ladder never visited is an error, which names the
                          ladder that was walked                    [default: off]

ALGORITHM:
  -N, --norm TYPE         NONE | VC | VC_SQRT | KR | SCALE        [default: NONE]
  -T, --threshold V       Clamp on distance-normalised values        [default: 10]
  -w, --window KB         Sliding-window width in kb                  [default: 5]
  -d, --distance BP       Minimum distance from the diagonal          [default: 0]
  -u, --upperlim BP       Maximum distance from the diagonal    [default: 0 = all]
  -q, --qvalue V          FDR threshold (0 disables filtering)      [default: 0.05]
  -p, --pcalculation 0|1  Compute p-values and write the q-value track and the
                          q-filtered track. Off by default: it costs a pass at
                          the finest resolution and nothing in the paper uses
                          the filtered output                        [default: 0]
  -E, --endZ BP           Mean-centre each chromosome at resolutions >= this.
                          Omit or pass 0 to disable                [default: 1000]
  -V, --vc-filter 0|1     Drop bins with no coverage, by CRUSH's own row
                          marginal (centromeres, telomeres, unmappable
                          regions). 0 reproduces pre-filter output   [default: 1]
  -M, --shifter-mode M    coarse | highpass                     [default: coarse]
      --coarse-interp 0|1 Interpolate the coarse track between bin centres when
                          re-centring, instead of holding it piecewise-constant.
                          A step function puts a discontinuity at every coarse
                          bin edge which enters the output at full magnitude and
                          compounds down the ladder                  [default: 1]
      --reverse-pass 0|1  Within each round, score the target resolution first
                          and then re-run every resolution seeded from THAT
                          result, rather than from the previous round's merge.
                          This is v1's reprocess_resolutions_with_shifter; it is
                          what makes the walk 'forward-and-reverse'  [default: 1]

INITIALISATION:
      --seed-select S     both | annotation | eigen. 'both' walks from the -a/-b
                          annotation seeds and from the eigenvector seeds, and
                          keeps, per chromosome, the result with the stronger
                          compartments. On shallow maps a seed can settle into
                          a weaker alternative partition on single chromosomes;
                          across 16 tested maps the stronger result never
                          needed an orientation flip and was never weaker than
                          the result from random seeds.
                          Doubles the run time                    [default: both]
      --orient 0|1        After the walk, flip any chromosome whose track
                          anti-correlates with gene density minus low-GC density,
                          so A is the gene-rich, GC-rich side      [default: 1]
      --eigen-init 0|1    Compute a principal component to seed round 0. 0 seeds
                          from the -a/-b annotations alone, which is what CRUSH
                          v1.81 did (it had no --eigenres). Required for maps too
                          sparse for a correlation matrix -- every Dip-C single
                          cell returns an undefined PC on every chromosome
                                                                     [default: 1]
  -e, --eigenfile FILE    Use this bedGraph to seed round 0 instead of
                          computing eigenvectors internally

SEEDING:
      --seed-res BP       Resolution at which seeds are maintained, independent
                          of the pass being scored; 0 uses the pass's own
                          resolution. Re-binning seeds to a coarse pass
                          collapses each coarse bin to a single label and throws
                          away the sub-structure the finer rounds established
                          -- v1 kept its seeds at the 500 bp GC granularity
                          throughout and only relabelled them      [default: 500]
      --seed-damping D    Fraction of the previous round's seed state carried
                          forward, 0..1. 0 replaces seeds outright each round;
                          higher values blend, so a small data difference
                          compounds more slowly across the walk      [default: 0]
      --seed-anchor 0|1   Keep the seed SET fixed to the positions the initial
                          annotation supplied (gene bodies, low-GC bins) and
                          only relabel them A/B from the previous round, rather
                          than seeding from every bin. v1 did this by accident:
                          B23 made its seed union collapse to the original
                          annotation every time                      [default: 0]
      --seed-frac F       Fraction of bins used as seeds each round, taken from
                          the extremes of the score distribution. 1.0 seeds
                          every bin; 0.5 seeds the top and bottom quartiles and
                          drops the ambiguous middle, which damps the
                          round-to-round noise the walk otherwise
                          amplifies                                [default: 1.0]
      --seed-split S      sign | median. Where to split re-derived seeds into A
                          and B each round. 'sign' uses a FIXED threshold of
                          zero, identical for every sample. 'median' uses a
                          per-sample threshold, so bins near it can land in
                          different seed classes in two replicates and the walk
                          amplifies that: measured on a Wi38 replicate pair at
                          25 kb, sign gives sqrtMSD 12.54 against median's
                          31.55 (the original pipeline gives 10.80). Median is
                          retained only as the degeneracy fallback, which fires
                          per chromosome when a seed class collapses below 5%
                                                                  [default: sign]

COMPATIBILITY:
      --v1                Run a line-for-line port of the bash v1 pipeline
                          that produced the published tracks, bugs included,
                          instead of the corrected algorithm. Defaults to that
                          pipeline's settings (-m 2500000 -w 5 -T 10); flags
                          given explicitly still override them. Add
                          --eigen-init 0 for runs that predate eigenvector
                          seeding.                              [default: off]
      --reshift 0|1       --v1 only: bash v1's Final Shifter, re-centring
                          each rung on the previous one         [default: 0]
      --smoothing 0|1     --v1 only: bash v1's 3-point smoothing after the
                          Final Shifter                         [default: 0]

MISC:
  -c, --cpu N             Worker threads                            [default: 24]
  -v, --verbose           Verbose progress reporting
  -h, --help              Short help
  -H, --help-advanced     This message

SHIFTER MODES:
  highpass  Subtract the local rolling mean. Reproduces CRUSH bash v1.0, which
            declared a coarse-resolution input but never read it, so this was
            the only behaviour ever actually available. Removes slowly varying
            regional signal.
  coarse    Subtract (rolling mean - coarse score), re-centring the finer track
            on the coarser one and preserving regional structure."
    )
}

/// Parse a --keepres spec: "all", or a comma-separated list of resolutions.
fn parse_keepres(spec: &str) -> Result<KeepRes, String> {
    let t = spec.trim();
    if t.eq_ignore_ascii_case("all") {
        return Ok(KeepRes::All);
    }
    if t.eq_ignore_ascii_case("none") || t.is_empty() {
        return Ok(KeepRes::None);
    }
    let mut out: Vec<u64> = Vec::new();
    for part in t.split(',') {
        let p = part.trim();
        if p.is_empty() {
            continue;
        }
        let v: u64 = p.parse().map_err(|_| {
            format!("--keepres takes 'all' or a comma-separated list of resolutions, got '{p}'")
        })?;
        if v == 0 {
            return Err("--keepres resolutions must be > 0".into());
        }
        if !out.contains(&v) {
            out.push(v);
        }
    }
    if out.is_empty() {
        return Err("--keepres was given no resolutions".into());
    }
    out.sort_unstable_by(|a, b| b.cmp(a));
    Ok(KeepRes::Only(out))
}

fn parse_args() -> Result<Config, String> {
    let mut cfg = Config::default();
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.iter().any(|a| a == "-H" || a == "--help-advanced") {
        println!("{}", usage_advanced());
        process::exit(0);
    }
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        println!("{}", usage());
        process::exit(0);
    }

    // --v1 is a preset, applied before any other flag is read so that an
    // explicit flag anywhere on the line still overrides it.
    if args.iter().any(|a| a == "--v1") {
        apply_v1_preset(&mut cfg);
    }

    // Every flag except -v takes exactly one value. Unknown flags are a hard
    // error rather than being silently discarded (the bash v1 behaviour).
    let mut i = 0;
    while i < args.len() {
        let flag = args[i].clone();
        let mut need = |i: usize| -> Result<String, String> {
            args.get(i + 1)
                .cloned()
                .ok_or_else(|| format!("option '{flag}' requires a value"))
        };
        match args[i].as_str() {
            "-v" | "--verbose" => {
                cfg.verbose = true;
                i += 1;
                continue;
            }
            "--v1" => {
                cfg.v1 = true;
                i += 1;
                continue;
            }
            "-i" | "--hic" => cfg.hic = PathBuf::from(need(i)?),
            "-g" | "--genomesize" => cfg.sizes = PathBuf::from(need(i)?),
            "-a" | "--initialA" => cfg.initial_a = PathBuf::from(need(i)?),
            "-b" | "--initialB" => cfg.initial_b = PathBuf::from(need(i)?),
            "-o" | "--outpre" => cfg.outpre = need(i)?,
            "-O" | "--outdir" => cfg.outdir = PathBuf::from(need(i)?),
            "-N" | "--norm" => cfg.norm = need(i)?,
            "-e" | "--eigenfile" => cfg.eigenfile = Some(PathBuf::from(need(i)?)),
            "--eigen-init" => {
                cfg.eigen_init = need(i)? != "0";
                cfg.seed_both = false;
            }
            "--seed-select" => {
                let v = need(i)?;
                match v.as_str() {
                    "both" => cfg.seed_both = true,
                    "annotation" => { cfg.seed_both = false; cfg.eigen_init = false; }
                    "eigen" | "eigenvector" => { cfg.seed_both = false; cfg.eigen_init = true; }
                    other => return Err(format!("--seed-select must be both, annotation or eigen, got '{other}'")),
                }
            }
            "--orient" => cfg.orient = need(i)? != "0",
            "-E" | "--endZ" => {
                let v = parse_num(&need(i)?, "-E")?;
                cfg.end_z = if v == 0 { None } else { Some(v) };
            }
            "-V" | "--vc-filter" => cfg.vc_filter = need(i)? != "0",
            "-B" | "--bigwig" => cfg.bigwig = need(i)? != "0",
            "--keep-bedgraph" => cfg.keep_bedgraph = need(i)? != "0",
            "--chr-prefix" => cfg.chr_prefix = need(i)? != "0",
            "--seed-frac" => cfg.seed_frac = parse_f64(&need(i)?, "--seed-frac")?,
            "--seed-res" => cfg.seed_res = parse_num(&need(i)?, "--seed-res")?,
            "--seed-anchor" => cfg.seed_anchor = need(i)? != "0",
            "--reverse-pass" => cfg.reverse_pass = need(i)? != "0",
            "--coarse-interp" => cfg.coarse_interp = need(i)? != "0",
            "--seed-damping" => cfg.seed_damping = parse_f64(&need(i)?, "--seed-damping")?,
            "--clamp" => {
                let v = parse_f64(&need(i)?, "--clamp")?;
                if v < 0.0 {
                    return Err("--clamp must be >= 0 (0 disables the clamp)".into());
                }
                cfg.clamp = v;
            }
            "--keepres" => {
                // A list typed with spaces after the commas arrives as several
                // argv entries. Glue back any that are purely digits and commas
                // so `--keepres 100000, 10000, 5000` means what it looks like,
                // without quoting. Anything else ends the list, so a following
                // flag is still parsed as a flag.
                let mut spec = need(i)?;
                let mut extra = 0usize;
                while spec.trim_end().ends_with(',') {
                    match args.get(i + 2 + extra) {
                        Some(n)
                            if !n.is_empty()
                                && n.chars().all(|c| c.is_ascii_digit() || c == ',') =>
                        {
                            spec.push_str(n);
                            extra += 1;
                        }
                        _ => break,
                    }
                }
                cfg.keep_res = parse_keepres(&spec)?;
                i += extra;
            }
            "--seed-split" => {
                let v = need(i)?;
                cfg.seed_split_median = match v.as_str() {
                    "median" => true,
                    "sign" => false,
                    other => {
                        return Err(format!(
                            "--seed-split must be 'median' or 'sign', got '{other}'"
                        ))
                    }
                };
            }
            "-M" | "--shifter-mode" => {
                let v = need(i)?;
                cfg.shifter_mode = match v.as_str() {
                    "highpass" => ShifterMode::Highpass,
                    "coarse" => ShifterMode::Coarse,
                    other => {
                        return Err(format!(
                            "--shifter-mode must be 'highpass' or 'coarse', got '{other}'"
                        ))
                    }
                };
            }
            "-r" | "--res" => cfg.res = parse_num(&need(i)?, "-r")?,
            "-m" | "--maxres" => cfg.coarsest_res = parse_num(&need(i)?, "-m")?,
            "-Z" | "--eigenres" => cfg.eigen_res = parse_num(&need(i)?, "-Z")?,
            "-d" | "--distance" => cfg.min_distance = parse_num(&need(i)?, "-d")?,
            "-u" | "--upperlim" => cfg.max_distance = parse_num(&need(i)?, "-u")?,
            "-c" | "--cpu" => cfg.cpu = parse_num(&need(i)?, "-c")? as usize,
            "-T" | "--threshold" => cfg.threshold = parse_f64(&need(i)?, "-T")?,
            "-w" | "--window" => cfg.window_kb = parse_f64(&need(i)?, "-w")?,
            "-q" | "--qvalue" => cfg.qthresh = parse_f64(&need(i)?, "-q")?,
            "-p" | "--pcalculation" => cfg.pcalculation = need(i)? != "0",
            "--reshift" => cfg.reshift = need(i)? != "0",
            "--smoothing" => cfg.smoothing = need(i)? != "0",
            other => return Err(format!("unrecognised option '{other}' (see --help)")),
        }
        i += 2;
    }

    if cfg.hic.as_os_str().is_empty() {
        return Err("missing required option -i/--hic".into());
    }
    for (p, what) in [
        (&cfg.sizes, "-g/--genomesize"),
        (&cfg.initial_a, "-a/--initialA"),
        (&cfg.initial_b, "-b/--initialB"),
    ] {
        if p.as_os_str().is_empty() {
            return Err(format!("missing required option {what}"));
        }
        if !p.exists() {
            return Err(format!("{what}: file not found: {}", p.display()));
        }
    }
    if !cfg.hic.exists() {
        return Err(format!("-i/--hic: file not found: {}", cfg.hic.display()));
    }
    if cfg.cpu == 0 {
        cfg.cpu = 1;
    }
    Ok(cfg)
}

/// Defaults of the published v1 runs (their CRUSHparamters.txt dumps:
/// --maxres 2500000 --window 5 --threshold 10). Applied before the other
/// flags are read, so any of them given explicitly still wins.
fn apply_v1_preset(cfg: &mut Config) {
    cfg.coarsest_res = 2_500_000;
    cfg.window_kb = 5.0;
    cfg.threshold = 10.0;
}

fn parse_num(s: &str, flag: &str) -> Result<u64, String> {
    s.parse()
        .map_err(|_| format!("option '{flag}' expects a non-negative integer, got '{s}'"))
}
fn parse_f64(s: &str, flag: &str) -> Result<f64, String> {
    s.parse()
        .map_err(|_| format!("option '{flag}' expects a number, got '{s}'"))
}

// ---------------------------------------------------------------------------
// Pipeline
// ---------------------------------------------------------------------------

/// A whole-genome bedGraph track: chrom -> [(start, end, value)].
type Track = HashMap<String, Vec<(u64, u64, f64)>>;

/// A/B seed bins per chromosome, as (span_of_each_seed, seed_start_positions).
/// The span matters: seeds are intervals, and a coarse seed must mark every
/// finer bin it covers when the pipeline walks down to a smaller resolution.
type Seeds = HashMap<String, (u64, Vec<u64>)>;

fn main() {
    let cfg = match parse_args() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("crush-rs: {e}");
            process::exit(2);
        }
    };
    if let Err(e) = run(&cfg) {
        eprintln!("crush-rs: {e}");
        process::exit(1);
    }
}

fn run(cfg_in: &Config) -> Result<(), String> {
    let mut cfg_owned: Config = cfg_in.clone();
    let cfg: &mut Config = &mut cfg_owned;
    rayon::ThreadPoolBuilder::new()
        .num_threads(cfg.cpu)
        .build_global()
        .ok();

    let sizes = ChromSizes::load(&cfg.sizes).map_err(|e| e.to_string())?;
    let genes = BedSet::load(&cfg.initial_a).map_err(|e| e.to_string())?;
    // -b accepts either a Bbins BED or a genome FASTA. bash decides by looking
    // for a '>' on the first line; do the same, and compute GC-based B seeds
    // rather than silently parsing a FASTA as BED (which yields 0 intervals and
    // then a confusing failure three stages later).
    let bbins = if gc::looks_like_fasta(&cfg.initial_b) {
        let bin = if cfg.res >= 500 { 500 } else { cfg.res };
        println!("      -b is a FASTA: deriving low-GC B seeds at {bin} bp");
        gc::bbins_from_fasta(&cfg.initial_b, &sizes, bin).map_err(|e| e.to_string())?
    } else {
        BedSet::load(&cfg.initial_b).map_err(|e| e.to_string())?
    };

    println!("crush-rs {VERSION}");
    println!("  hic        : {}", cfg.hic.display());
    println!("  chromosomes: {}", sizes.len());
    println!("  initial A  : {} intervals", genes.total());
    println!("  initial B  : {} intervals", bbins.total());

    if cfg.v1 {
        println!("  mode       : --v1 (port of the published bash v1 pipeline)");
        let v1cfg = v1::V1Config {
            hic: cfg.hic.clone(),
            sizes: sizes.clone(),
            genes: v1::load_ivs(&cfg.initial_a)?,
            bbins: if gc::looks_like_fasta(&cfg.initial_b) {
                let (b, t) = v1::bbins_from_fasta(&cfg.initial_b, &sizes)?;
                println!(
                    "      v1 GC seeding: threshold {t}, {} B bins",
                    b.values().map(|v| v.len()).sum::<usize>()
                );
                b
            } else {
                v1::load_ivs(&cfg.initial_b)?
            },
            res: cfg.res,
            coarsest_res: cfg.coarsest_res,
            eigen_res: cfg.eigen_res,
            eigen_init: cfg.eigen_init && cfg.eigenfile.is_none(),
            window_kb: cfg.window_kb,
            threshold: cfg.threshold,
            reshift: cfg.reshift,
            smoothing: cfg.smoothing,
            outdir: cfg.outdir.clone(),
            outpre: cfg.outpre.clone(),
            verbose: cfg.verbose,
        };
        let written = v1::run_v1(&v1cfg, cfg.cpu)?;
        if let Some(last) = written.last() {
            println!("      {}", last.display());
            if cfg.bigwig {
                bedgraph_to_bigwig(last, &sizes, cfg.keep_bedgraph, cfg.chr_prefix);
            }
        }
        println!("\nFinished.");
        return Ok(());
    }

    let available = HicFile::open(&cfg.hic)
        .map_err(|e| e.to_string())?
        .bp_resolutions
        .clone();
    let mut ladder: Vec<u64> = available
        .iter()
        .map(|&r| r as u64)
        .filter(|&r| r >= cfg.res && r <= cfg.coarsest_res)
        .collect();
    ladder.sort_unstable_by(|a, b| b.cmp(a)); // coarse -> fine

    if ladder.is_empty() {
        return Err(format!(
            "no resolutions stored in {} between -r {} and -m {} (file has {:?})",
            cfg.hic.display(),
            cfg.res,
            cfg.coarsest_res,
            available
        ));
    }
    println!("  ladder     : {ladder:?}");
    if let KeepRes::Only(want) = &cfg.keep_res {
        let missing: Vec<u64> = want.iter().copied().filter(|r| !ladder.contains(r)).collect();
        if !missing.is_empty() {
            return Err(format!(
                "--keepres asked for {missing:?}, which the walk never visits. \
                 The ladder for -r {} -m {} on this file is {ladder:?}",
                cfg.res, cfg.coarsest_res
            ));
        }
    }

    match cfg.end_z {
        Some(e) => println!("  mean-centre: resolutions >= {e} bp"),
        None => println!("  mean-centre: off"),
    }
    let cfg: &Config = cfg;

    // Walk from one or two seeds; see walk() and select_by_strength().
    let use_both = cfg.seed_both && cfg.eigenfile.is_none();
    let wo = if use_both {
        println!("\n  seeding    : annotation AND eigenvector; each chromosome keeps the walk");
        println!("               with the stronger compartments (--seed-select)");
        let ann = walk(cfg, &sizes, &genes, &bbins, &ladder, false)?;
        let eig = walk(cfg, &sizes, &genes, &bbins, &ladder, true)?;
        select_by_strength(cfg, &sizes, &available, ann, eig)
    } else {
        walk(cfg, &sizes, &genes, &bbins, &ladder, cfg.eigen_init)?
    };
    let WalkOut { mut merged, mut merged_tracks, pvals_at_finest, mask_cache } = wo;
    if cfg.orient {
        orient_by_prior(&mut merged, &mut merged_tracks, &sizes, &genes, &bbins);
    }
    let finest = *ladder.last().unwrap();

    println!("\n[4/4] Writing output");
    std::fs::create_dir_all(&cfg.outdir).map_err(|e| e.to_string())?;

    // Rescale so mean |score| is 100, matching the bash GIave normalisation.
    let giave = giave_of(&merged);
    let clamp = if cfg.clamp > 0.0 { Some(cfg.clamp) } else { None };
    match clamp {
        Some(c) => println!("      clamping output scores to ±{c:.0}"),
        None => println!("      output score clamp disabled"),
    }

    let crush_path = cfg
        .outdir
        .join(format!("{}mergedCrush_{}.bedgraph", cfg.outpre, finest));
    write_track(&crush_path, &merged, &sizes, Some(giave / 100.0), clamp)?;
    let mut written: Vec<PathBuf> = vec![crush_path.clone()];
    println!("      {}", crush_path.display());

    if let Some(pv) = pvals_at_finest {
        let mut keys: Vec<(String, u64, u64)> = Vec::new();
        let mut ps: Vec<f64> = Vec::new();
        let mut chroms: Vec<&String> = pv.keys().collect();
        chroms.sort();
        for c in chroms {
            for &(s, e, p) in &pv[c] {
                keys.push((c.clone(), s, e));
                ps.push(p);
            }
        }
        let qs = stats::benjamini_hochberg(&ps);

        let mut qtrack: Track = HashMap::new();
        for ((c, s, e), q) in keys.iter().zip(&qs) {
            qtrack.entry(c.clone()).or_default().push((*s, *e, *q));
        }
        let q_path = cfg
            .outdir
            .join(format!("{}mergedqvalue_{}.bedgraph", cfg.outpre, finest));
        write_track(&q_path, &qtrack, &sizes, None, None)?;
        written.push(q_path.clone());
        println!("      {}", q_path.display());

        if cfg.qthresh > 0.0 {
            let mut lookup: HashMap<(&str, u64), f64> = HashMap::new();
            for ((c, s, _), q) in keys.iter().zip(&qs) {
                lookup.insert((c.as_str(), *s), *q);
            }
            let mut filtered: Track = HashMap::new();
            for (c, v) in &merged {
                for &(s, e, x) in v {
                    if lookup
                        .get(&(c.as_str(), s))
                        .map_or(false, |&q| q <= cfg.qthresh)
                    {
                        filtered.entry(c.clone()).or_default().push((s, e, x));
                    }
                }
            }
            let f_path = cfg.outdir.join(format!(
                "{}mergedCrush_{}_qfiltered.bedgraph",
                cfg.outpre, finest
            ));
            write_track(&f_path, &filtered, &sizes, Some(giave / 100.0), clamp)?;
            written.push(f_path.clone());
            println!("      {}", f_path.display());
        }
    }

    // --keepres: the walk's per-resolution tracks, which the merge above
    // otherwise consumes. Each is coverage-filtered and rescaled on its own
    // terms so it stands alone in a browser next to the merged track.
    let keep: Vec<u64> = match &cfg.keep_res {
        KeepRes::None => Vec::new(),
        KeepRes::All => ladder.clone(),
        KeepRes::Only(v) => v.clone(),
    };
    if !keep.is_empty() {
        for &res in &keep {
            let Some((_, track)) = merged_tracks.iter().find(|(r, _)| *r == res) else {
                eprintln!(
                    "      NOTE: --keepres {res} skipped — the final round produced                      no track at that resolution (its pass found no usable seeds)"
                );
                continue;
            };
            let mut t = track.clone();
            if cfg.vc_filter {
                if let Some(m) = mask_cache.get(&res) {
                    t = filter_by_coverage(&t, m, res);
                }
            }
            let g = giave_of(&t);
            let path = cfg
                .outdir
                .join(format!("{}Crush_{}.bedgraph", cfg.outpre, res));
            write_track(&path, &t, &sizes, Some(g / 100.0), clamp)?;
            written.push(path.clone());
            println!("      {}", path.display());
        }
    }

    if cfg.bigwig {
        for bg in &written {
            if let Some(bw) =
                bedgraph_to_bigwig(bg, &sizes, cfg.keep_bedgraph, cfg.chr_prefix)
            {
                println!("      {}", bw.display());
            }
        }
    }

    // Does the map actually support the track we just wrote? CRUSH rescales its
    // output to a fixed mean magnitude, so it emits a compartment-shaped track
    // for any input; this is the check that the contacts agree. Advisory only —
    // a low score never fails the run, because the answer "this map has no
    // compartments" is a legitimate result, not an error.
    println!("\n  compartment check");
    match diagnostic::evaluate(&cfg.hic, &cfg.norm, &sizes, &merged, &available) {
        Some(rep) => {
            println!(
                "      (AA+BB)/2AB over quartiles of the merged track, scored at {} kb,",
                rep.scale / 1000
            );
            println!("      contacts >= 1 Mb apart; 1.0 = the track shifted along itself");
            for c in &rep.per_chrom {
                println!(
                    "        {:<6} {:>6.3}   null {:>5.3}   z {:>7.1}   {:>7} bins",
                    c.chrom, c.strength, c.null, c.z, c.bins
                );
            }
            println!(
                "      median over {} chromosomes: {:.3}  [{}]",
                rep.per_chrom.len(),
                rep.median,
                rep.verdict()
            );
            if rep.median < diagnostic::WARN_BELOW {
                eprintln!(
                    "\n  WARNING: compartment strength {:.3} is at chance (< {:.2}).",
                    rep.median,
                    diagnostic::WARN_BELOW
                );
                eprintln!(
                    "      The contact map does not support the A/B track CRUSH just wrote.\n\
                     \x20     CRUSH always produces a signed track, so a plausible-looking\n\
                     \x20     bedGraph is not evidence that compartments exist. Likely causes:\n\
                     \x20       - this map genuinely has no compartments (mitotic, or similar)\n\
                     \x20       - too few contacts to resolve them\n\
                     \x20       - the wrong normalisation (-N); an unbalanced matrix can bury them\n\
                     \x20     Treat the track as unvalidated before interpreting it."
                );
            } else if rep.median < diagnostic::PASS_AT {
                eprintln!(
                    "\n  NOTE: compartment strength {:.3} is between {:.2} and {:.2} — weak or\n\
                     \x20     unusual compartmentalisation. Real but short-range compartments\n\
                     \x20     (e.g. pachytene) land here. Worth checking before interpreting.",
                    rep.median,
                    diagnostic::WARN_BELOW,
                    diagnostic::PASS_AT
                );
            }
        }
        None => println!("      skipped (no chromosome had enough coverage to score)"),
    }

    println!("\nFinished.");
    Ok(())
}


/// Everything the walk produces that the output stage needs.
struct WalkOut {
    merged: Track,
    merged_tracks: Vec<(u64, Track)>,
    pvals_at_finest: Option<Track>,
    mask_cache: HashMap<u64, HashMap<String, Vec<bool>>>,
}

/// Seeding, the forward-and-reverse walk, and the merge, from one seed:
/// the eigenvector (`eigen_init`) or the -a/-b annotations.
fn walk(
    cfg: &Config,
    sizes: &ChromSizes,
    genes: &BedSet,
    bbins: &BedSet,
    ladder: &Vec<u64>,
    eigen_init: bool,
) -> Result<WalkOut, String> {
    let (sizes, genes, bbins) = (sizes.clone(), genes, bbins);
    let ladder = ladder.clone();
    if eigen_init || cfg.eigenfile.is_some() {
        println!("\n[1/4] Eigenvector initialisation at {} bp", cfg.eigen_res);
    } else {
        println!("\n[1/4] Annotation initialisation at {} bp", ladder[0]);
    }

    // Round 0 is seeded from a compartment track: either one supplied with -e
    // (useful for isolating the walk from the eigenvector step when comparing
    // implementations) or one computed here from the best compartment-matching
    // principal component. Either way the seeds are derived with the SAME
    // majority-binning rule used for every later round.
    let eigen_track: Track = match &cfg.eigenfile {
        Some(path) => {
            println!("      using supplied eigen track: {}", path.display());
            load_track(path)?
        }
        None if !eigen_init => {
            println!("      eigenvector init disabled: seeding from -a/-b annotations");
            annotation_track(&sizes, &genes, &bbins, ladder[0])
        }
        None => eigen_track_from_hic(cfg, &sizes, &genes, &bbins)?,
    };
    let (a_seed, b_seed) = seeds_from_net(&track_to_net(&eigen_track, ladder[0], None), ladder[0], cfg.seed_split_median, cfg.seed_frac);
    println!(
        "      seeded {} A bins and {} B bins across {} chromosomes",
        a_seed.values().map(|v| v.1.len()).sum::<usize>(),
        b_seed.values().map(|v| v.1.len()).sum::<usize>(),
        a_seed.len()
    );

    // ---- Forward-and-reverse resolution walking ------------------------------
    //
    // Calling compartments at the finest resolution is not a single sweep down
    // the ladder. For each target resolution the walk restarts from the
    // coarsest one, re-deriving every coarser call from the seeds the previous
    // round produced:
    //
    //   target 1 Mb  : PC init          -> 1 Mb
    //   target 500 kb: seeds from above -> 1 Mb -> 500 kb
    //   target 250 kb: seeds from above -> 1 Mb -> 500 kb -> 250 kb
    //
    // and on reaching the finest resolution the whole ladder is recomputed once
    // more with that resolution's output as the initialisation. So a run over
    // n resolutions performs O(n^2) scoring passes by design.
    //
    // The merged output for a target resolution is built from the inner sweep
    // that ended at it — not from the tracks produced on the way down.
    println!("\n[2/4] Forward-and-reverse resolution walking");
    let mut pvals_at_finest: Option<Track> = None;
    let mut merged_tracks: Vec<(u64, Track)> = Vec::new();

    // Round k initialises every one of its passes from the MERGED track that
    // round k-1 produced, re-binned by majority at each resolution. Round 0 is
    // initialised from the eigenvector seeds.
    let finest = *ladder.last().unwrap();
    let mut merged_prev: Option<Track> = None;
    // Coverage masks are per resolution and reused across rounds.
    let mut mask_cache: HashMap<u64, HashMap<String, Vec<bool>>> = HashMap::new();
    // Carried-forward seed state, blended each round when --seed-damping > 0.
    let mut seed_state: HashMap<String, HashMap<u64, f64>> = HashMap::new();

    // Positions the initial annotation supplied, at seed resolution. Used only
    // when --seed-anchor is set; see seeds_from_track.
    let anchor_set: Option<HashMap<String, std::collections::HashSet<u64>>> = if cfg.seed_anchor {
        let sres = if cfg.seed_res > 0 { cfg.seed_res } else { *ladder.last().unwrap() };
        let mut m: HashMap<String, std::collections::HashSet<u64>> = HashMap::new();
        for seeds in [&a_seed, &b_seed] {
            for (chrom, (span, list)) in seeds.iter() {
                let e = m.entry(chrom.clone()).or_default();
                for &p in list {
                    let lo = (p / sres) * sres;
                    let hi = ((p + (*span).max(1) - 1) / sres) * sres;
                    let mut b = lo;
                    while b <= hi.max(lo) {
                        e.insert(b);
                        b += sres;
                    }
                }
            }
        }
        let n: usize = m.values().map(|v| v.len()).sum();
        println!("  seed anchor : {n} positions from the initial annotation");
        Some(m)
    } else {
        None
    };

    for (i, &target) in ladder.iter().enumerate() {
        let t_outer = std::time::Instant::now();

        let mut sweep: Vec<(u64, Track)> = Vec::new();
        let mut coarse_track: Option<Track> = None;

        // FORWARD leg: score this round's target resolution, seeded from the
        // previous round's merge. Its result then seeds the reverse leg below.
        //
        // v1 does `process_resolution $myres` and only afterwards
        // `reprocess_resolutions_with_shifter $myres $i`, which re-runs every
        // resolution 0..i against combined_new{A,B}bins -- and those were just
        // rewritten from the target resolution's own calls. Seeding the sweep
        // from the previous round's merge instead skips that feedback: the
        // coarse passes of the final round never see the finest calls at all.
        let mut round_seed: Option<Track> = None;
        if cfg.reverse_pass && i > 0 {
            let sres0 = if cfg.seed_res > 0 { cfg.seed_res } else { target };
            let (a0, b0) = match &merged_prev {
                Some(m) => seeds_from_net(
                    &track_to_net(m, sres0, anchor_set.as_ref()),
                    sres0, cfg.seed_split_median, cfg.seed_frac),
                None => (a_seed.clone(), b_seed.clone()),
            };
            if !a0.is_empty() && !b0.is_empty() {
                let w0 = window_bins(cfg.window_kb, target);
                let (raw0, _m0, cov0) =
                    score_all_chromosomes(cfg, &sizes, target, w0, &a0, &b0);
                mask_cache.entry(target).or_insert(cov0);
                round_seed = Some(raw0);
            }
        }

        for j in 0..=i {
            let res = ladder[j];
            let window = window_bins(cfg.window_kb, res);

            // Seeds live at `seed_res`, not at the pass's resolution: a coarse
            // pass would otherwise collapse every coarse bin to one label and
            // discard what the finer rounds resolved. `to_bins` expands them
            // back onto whatever grid this pass uses, so a coarse bin holding
            // both A and B seeds contributes to both moment sets -- which is
            // what v1 did, since its seed set stayed the original 500 bp low-GC
            // bins and gene intervals at every resolution.
            let sres = if cfg.seed_res > 0 { cfg.seed_res } else { res };
            let seed_src = round_seed.as_ref().or(merged_prev.as_ref());
            let (a_r, b_r) = match seed_src {
                Some(m) => {
                    let fresh = track_to_net(m, sres, anchor_set.as_ref());
                    if cfg.seed_damping > 0.0 {
                        damp_net(&mut seed_state, fresh, cfg.seed_damping);
                    } else {
                        seed_state = fresh;
                    }
                    seeds_from_net(&seed_state, sres, cfg.seed_split_median, cfg.seed_frac)
                }
                None => (a_seed.clone(), b_seed.clone()),
            };
            if a_r.is_empty() || b_r.is_empty() {
                eprintln!("      NOTE: no usable seeds at {res} bp — skipping this pass");
                continue;
            }

            let (raw, moments, covered) =
                score_all_chromosomes(cfg, &sizes, res, window, &a_r, &b_r);
            mask_cache.entry(res).or_insert(covered);

            let shifted = if j == 0 {
                raw
            } else {
                let masks = if cfg.vc_filter {
                    mask_cache.get(&res)
                } else {
                    None
                };
                apply_shifter(
                    cfg,
                    &sizes,
                    &raw,
                    coarse_track.as_ref(),
                    res,
                    ladder[j - 1],
                    masks,
                )
            };

            if cfg.pcalculation && j == i && i == ladder.len() - 1 {
                pvals_at_finest = Some(compute_pvalues(&moments, res));
            }

            coarse_track = Some(shifted.clone());
            sweep.push((res, shifted));
        }

        if sweep.is_empty() {
            return Err(format!("no resolution produced a track for target {target} bp"));
        }

        let mut round_merged = merge_tracks(&sweep, target);
        if cfg.vc_filter {
            if let Some(m) = mask_cache.get(&target) {
                round_merged = filter_by_coverage(&round_merged, m, target);
            }
        }
        println!(
            "      target {:>9} bp : {} pass(es), {} merged bins  ({:.1}s)",
            target,
            i + 1,
            round_merged.values().map(|v| v.len()).sum::<usize>(),
            t_outer.elapsed().as_secs_f64()
        );

        merged_prev = Some(round_merged);
        merged_tracks = sweep;
    }
    let _ = finest;

    println!("\n[3/4] Merging {} resolutions", merged_tracks.len());
    let finest = *ladder.last().unwrap();
    let mut merged = merge_tracks(&merged_tracks, finest);
    if cfg.vc_filter {
        if let Some(m) = mask_cache.get(&finest) {
            let before: usize = merged.values().map(|v| v.len()).sum();
            merged = filter_by_coverage(&merged, m, finest);
            let after: usize = merged.values().map(|v| v.len()).sum();
            println!(
                "      VC filter: dropped {} uncovered bins ({} remain)",
                before - after,
                after
            );
        }
    }
    Ok(WalkOut { merged, merged_tracks, pvals_at_finest, mask_cache })
}

/// --seed-select both: per chromosome, keep the walk whose merged track has
/// the higher compartment strength ((AA+BB)/2AB, the same check printed at
/// the end of every run). Seeds can settle into a weaker alternative
/// partition on individual chromosomes of shallow maps; across 16 maps the
/// stronger of the two always matched what independent random seeds found.
/// The annotation walk is kept unless the eigenvector walk is stronger by more
/// than SEED_MARGIN, so chromosomes where both seeds reach the same partition
/// do not switch on noise (real disagreements in the 16-map test differed by
/// 0.5-1.2). Chromosomes the check cannot score keep the annotation walk.
const SEED_MARGIN: f64 = 0.05;

fn select_by_strength(
    cfg: &Config,
    sizes: &ChromSizes,
    available: &[u32],
    ann: WalkOut,
    eig: WalkOut,
) -> WalkOut {
    let score = |t: &Track| -> HashMap<String, f64> {
        diagnostic::evaluate(&cfg.hic, &cfg.norm, sizes, t, available)
            .map(|r| r.per_chrom.into_iter().map(|c| (c.chrom, c.strength)).collect())
            .unwrap_or_default()
    };
    let (sa, se) = (score(&ann.merged), score(&eig.merged));
    let mut use_eig: std::collections::HashSet<String> = std::collections::HashSet::new();
    println!("\n  seed selection (compartment strength: annotation / eigenvector)");
    for c in &sizes.names {
        let (a, e) = (sa.get(c).copied(), se.get(c).copied());
        let pick_e = matches!((a, e), (Some(a), Some(e)) if e > a + SEED_MARGIN);
        if pick_e {
            use_eig.insert(c.clone());
        }
        if let (Some(a), Some(e)) = (a, e) {
            let flag = if (a - e).abs() > 0.1 { "  <- seeds disagree" } else { "" };
            println!("        {:<6} {:>6.3} / {:>6.3}   {}{}", c, a, e,
                     if pick_e { "eigenvector" } else { "annotation" }, flag);
        }
    }
    let pick = |a: &Track, e: &Track| -> Track {
        let mut out = Track::new();
        for (c, v) in a {
            if !use_eig.contains(c) {
                out.insert(c.clone(), v.clone());
            }
        }
        for (c, v) in e {
            if use_eig.contains(c) {
                out.insert(c.clone(), v.clone());
            }
        }
        out
    };
    let merged = pick(&ann.merged, &eig.merged);
    let merged_tracks = ann
        .merged_tracks
        .iter()
        .map(|(r, ta)| {
            let te = eig.merged_tracks.iter().find(|(re, _)| re == r).map(|(_, t)| t);
            (*r, te.map_or_else(|| ta.clone(), |te| pick(ta, te)))
        })
        .collect();
    let pvals_at_finest = match (&ann.pvals_at_finest, &eig.pvals_at_finest) {
        (Some(a), Some(e)) => Some(pick(a, e)),
        (a, _) => a.clone(),
    };
    println!("      {} of {} chromosomes take the eigenvector walk", use_eig.len(), sizes.names.len());
    WalkOut { merged, merged_tracks, pvals_at_finest, mask_cache: ann.mask_cache }
}

/// --orient 1: A should be the gene-rich, GC-rich side. Per chromosome, the
/// merged track (100 kb means) is correlated with gene density minus low-GC
/// (B annotation) density, each z-scored; a negative correlation flips the
/// chromosome in every output track. The seed sets which side is called A,
/// so this guards the sign whichever seed won.
fn orient_by_prior(
    merged: &mut Track,
    merged_tracks: &mut Vec<(u64, Track)>,
    sizes: &ChromSizes,
    genes: &BedSet,
    bbins: &BedSet,
) {
    const W: u64 = 100_000;
    let lookup = |b: &BedSet, c: &str| -> Vec<bed::Interval> {
        let alt = c.strip_prefix("chr").map(|x| x.to_string()).unwrap_or_else(|| format!("chr{c}"));
        b.by_chrom.get(c).or_else(|| b.by_chrom.get(&alt)).cloned().unwrap_or_default()
    };
    let z = |v: Vec<f64>| -> Vec<f64> {
        let n = v.len() as f64;
        let m = v.iter().sum::<f64>() / n;
        let sd = (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n).sqrt();
        v.iter().map(|x| if sd > 0.0 { (x - m) / sd } else { 0.0 }).collect()
    };
    let mut flipped = Vec::new();
    for (ci, c) in sizes.names.iter().enumerate() {
        let Some(rows) = merged.get(c) else { continue };
        let n = (sizes.lengths[ci] / W + 1) as usize;
        let dens = |b: &BedSet| -> Vec<f64> {
            let mut d = vec![0.0; n];
            for iv in lookup(b, c) {
                let i = (((iv.start + iv.end) / 2) / W) as usize;
                if i < n {
                    d[i] += 1.0;
                }
            }
            d
        };
        let (g, l) = (z(dens(genes)), z(dens(bbins)));
        let prior: Vec<f64> = g.iter().zip(&l).map(|(a, b)| a - b).collect();
        let (mut sum, mut cnt) = (vec![0.0; n], vec![0.0; n]);
        for &(s, _, v) in rows {
            let i = (s / W) as usize;
            if i < n && v.is_finite() {
                sum[i] += v;
                cnt[i] += 1.0;
            }
        }
        let (mut xs, mut ys) = (Vec::new(), Vec::new());
        for i in 0..n {
            if cnt[i] > 0.0 {
                xs.push(sum[i] / cnt[i]);
                ys.push(prior[i]);
            }
        }
        if let Some(r) = eigen::pearson(&xs, &ys) {
            if r < 0.0 {
                flipped.push(format!("{c} (r = {r:+.2})"));
                for v in merged.get_mut(c).into_iter().flatten() {
                    v.2 = -v.2;
                }
                for (_, t) in merged_tracks.iter_mut() {
                    for v in t.get_mut(c).into_iter().flatten() {
                        v.2 = -v.2;
                    }
                }
            }
        }
    }
    if flipped.is_empty() {
        println!("\n  orientation: every chromosome already has A on the gene-rich, GC-rich side");
    } else {
        println!("\n  orientation: flipped {} chromosome(s) so A is gene-rich and GC-rich: {}",
                 flipped.len(), flipped.join(", "));
    }
}

/// Mean |score| over a track, the bash GIave normalisation constant. Dividing
/// by giave/100 puts the mean magnitude at 100.
fn giave_of(track: &Track) -> f64 {
    let (mut sum, mut n) = (0.0f64, 0usize);
    for v in track.values() {
        for &(_, _, x) in v {
            sum += x.abs();
            n += 1;
        }
    }
    if n == 0 || sum == 0.0 {
        1.0
    } else {
        sum / n as f64
    }
}

/// Bins per sliding window at this resolution: `int(window_kb / (res/1000)) + 1`.
fn window_bins(window_kb: f64, res: u64) -> usize {
    if window_kb <= 0.0 {
        return 1;
    }
    ((window_kb / (res as f64 / 1000.0)) as usize) + 1
}

/// Compute the eigenvector per chromosome and split it into A/B seed bins.
/// Round-0 seeding with no eigenvector, reproducing CRUSH v1.81 -- that version
/// had no --eigenres and seeded purely from the supplied A/B annotations.
///
/// Needed because a correlation-matrix principal component is undefined on maps
/// too sparse to fill it: every Dip-C single cell fails on every chromosome, so
/// the eigenvector path cannot produce single-cell compartments at all.
///
/// The value is (gene bins - low-GC bins) per bin, so the same sign split that
/// every later round uses applies unchanged.
fn annotation_track(
    sizes: &ChromSizes,
    genes: &BedSet,
    bbins: &BedSet,
    res: u64,
) -> Track {
    let mut out: Track = HashMap::new();
    for chrom in &sizes.names {
        let len = match sizes.length_of(chrom) {
            Some(l) => l,
            None => continue,
        };
        let gene_bins = bed::occupied_bins(genes, chrom, res);
        let gc_bins = bed::occupied_bins(bbins, chrom, res);
        let mut v: Vec<(u64, u64, f64)> = Vec::new();
        let mut s = 0u64;
        while s < len {
            let a = *gene_bins.get(&s).unwrap_or(&0) as f64;
            let b = *gc_bins.get(&s).unwrap_or(&0) as f64;
            if a > 0.0 || b > 0.0 {
                v.push((s, (s + res).min(len), a - b));
            }
            s += res;
        }
        if !v.is_empty() {
            out.insert(chrom.clone(), v);
        }
    }
    out
}

fn eigen_track_from_hic(
    cfg: &Config,
    sizes: &ChromSizes,
    genes: &BedSet,
    bbins: &BedSet,
) -> Result<Track, String> {
    let res = cfg.eigen_res;
    let out: Mutex<Track> = Mutex::new(HashMap::new());

    sizes.names.par_iter().for_each(|chrom| {
        let mut f = match HicFile::open(&cfg.hic) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("      NOTE: {chrom}: {e}");
                return;
            }
        };

        // Prefer VC_SQRT-normalised O/E, as the bash pipeline does; fall back to
        // unnormalised O/E for maps that carry no VC_SQRT vector.
        let contacts = match f.straw_cis(ValueKind::Oe, "VC_SQRT", chrom, res as u32) {
            Ok(c) if !c.is_empty() => c,
            _ => match f.straw_cis(ValueKind::Oe, "NONE", chrom, res as u32) {
                Ok(c) if !c.is_empty() => c,
                _ => {
                    if cfg.verbose {
                        eprintln!("      NOTE: no data for {chrom} at {res} bp — skipping");
                    }
                    return;
                }
            },
        };

        let gene_bins = bed::occupied_bins(genes, chrom, res);
        let gc_bins = bed::occupied_bins(bbins, chrom, res);
        let size = sizes.length_of(chrom).unwrap_or(0);

        let input = eigen::EigenInput {
            contacts: &contacts,
            res,
            chrom_size: size,
            gene_bins: &gene_bins,
            gc_bins: &gc_bins,
            n_components: 10,
            max_matrix_bytes: cfg.max_matrix_bytes,
        };

        match eigen::compute_eigen(&input) {
            Ok(r) => {
                if cfg.verbose {
                    println!(
                        "      {chrom}: PC{} (gene r={:.3}, gc r={:.3})",
                        r.chosen_pc, r.gene_corr, r.gc_corr
                    );
                }
                let rows: Vec<(u64, u64, f64)> =
                    r.track.into_iter().map(|(pos, val)| (pos, pos + res, val)).collect();
                out.lock().unwrap().insert(chrom.clone(), rows);
            }
            Err(e) => eprintln!("      NOTE: {chrom}: {e}"),
        }
    });

    let out = out.into_inner().unwrap();
    if out.is_empty() {
        return Err("eigenvector initialisation produced no track on any chromosome".into());
    }
    Ok(out)
}

/// Read a 4-column bedGraph into a Track, skipping any track header lines.
fn load_track(path: &PathBuf) -> Result<Track, String> {
    use std::io::BufRead;
    let f = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut out: Track = HashMap::new();
    for line in std::io::BufReader::new(f).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with("track") || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.trim_end().split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        let (Ok(s), Ok(e), Ok(v)) = (
            f[1].parse::<u64>(),
            f[2].parse::<u64>(),
            f[3].parse::<f64>(),
        ) else {
            continue;
        };
        out.entry(f[0].to_string()).or_default().push((s, e, v));
    }
    if out.is_empty() {
        return Err(format!("{}: no usable bedGraph rows", path.display()));
    }
    Ok(out)
}

type MomentMap = HashMap<String, HashMap<u64, (Moments, Moments)>>;

fn score_all_chromosomes(
    cfg: &Config,
    sizes: &ChromSizes,
    res: u64,
    window: usize,
    a_seed: &Seeds,
    b_seed: &Seeds,
) -> (Track, MomentMap, HashMap<String, Vec<bool>>) {
    let track: Mutex<Track> = Mutex::new(HashMap::new());
    let moms: Mutex<MomentMap> = Mutex::new(HashMap::new());
    let cov: Mutex<HashMap<String, Vec<bool>>> = Mutex::new(HashMap::new());

    sizes.names.par_iter().for_each(|chrom| {
        let Ok(mut f) = HicFile::open(&cfg.hic) else {
            return;
        };
        let contacts = match f.straw_cis(ValueKind::Observed, &cfg.norm, chrom, res as u32) {
            Ok(c) if !c.is_empty() => c,
            _ => return,
        };

        // Seeds carry a span (eigen_res at step 0, the seed resolution later).
        // The bash pipeline bins each seed INTERVAL across every bin it touches
        // (`for (i = int(s/r)*r; i <= int(e/r)*r; i += r)`), so a 100 kb seed
        // marks all 100 finer bins when res is 1 kb. Binning only the interval
        // start would silently drop 99% of the seeds at fine resolutions.
        let to_bins = |v: Option<&(u64, Vec<u64>)>| -> HashMap<u64, u32> {
            let mut m = HashMap::new();
            if let Some((span, list)) = v {
                for &p in list {
                    // Half-open: a seed covering [p, p+span) ends at the last
                    // bin containing p+span-1. Using p+span instead lands on
                    // the NEXT bin boundary and marks one bin too many, which
                    // silently doubled the seed sets whenever span == res.
                    let s = (p / res) * res;
                    let e = ((p + (*span).max(1) - 1) / res) * res;
                    let mut i = s;
                    while i <= e.max(s) {
                        *m.entry(i).or_insert(0) += 1;
                        i += res;
                    }
                }
            }
            m
        };
        let a_bins = to_bins(a_seed.get(chrom));
        let b_bins = to_bins(b_seed.get(chrom));
        if a_bins.is_empty() || b_bins.is_empty() {
            return;
        }

        let input = ChromScoreInput {
            contacts: &contacts,
            chrom_size: sizes.length_of(chrom).unwrap_or(0),
            res,
            min_distance: cfg.min_distance,
            max_distance: cfg.max_distance,
            threshold: cfg.threshold,
            a_bins: &a_bins,
            b_bins: &b_bins,
        };
        let out = score::score_chromosome(&input);
        let mut smoothed = score::sliding_window(&out.scores, res, window);

        // Per-chromosome mean-centring. The bash pipeline applies this at the
        // end of process_oppocheck_statement whenever `myres >= endZ`, and endZ
        // defaults to the finest resolution, so in practice it always runs.
        // Omitting it left crush-rs offset from bash by a constant (measured:
        // mean residual -44.8 at 1 Mb) which flipped the sign of every bin
        // scoring near zero.
        if cfg.end_z.map_or(false, |e| res >= e) && !smoothed.is_empty() {
            let mean: f64 =
                smoothed.iter().map(|r| r.2).sum::<f64>() / smoothed.len() as f64;
            if mean.is_finite() {
                for row in smoothed.iter_mut() {
                    row.2 -= mean;
                }
            }
        }

        // Debug hook: dump the per-bin A/B moments for one chromosome so they
        // can be diffed against bash's A2removedfile / B2removedfile.
        if let Ok(dbg) = std::env::var("CRUSH_DEBUG_CHROM") {
            if dbg == *chrom {
                let dump = |m: &HashMap<u64, Moments>, path: &str| {
                    let mut keys: Vec<_> = m.keys().copied().collect();
                    keys.sort_unstable();
                    let mut w = String::new();
                    for k in keys {
                        let v = m[&k];
                        w.push_str(&format!(
                            "{}\t{:.10e}\t{}\t{:.10e}\t{:.10e}\n",
                            k, v.sum * v.count, v.count as i64, v.mean(), v.sd()
                        ));
                    }
                    let _ = std::fs::write(format!("{path}_{res}.txt"), w);
                };
                dump(&out.a_moments, "rs_A2removed");
                dump(&out.b_moments, "rs_B2removed");
            }
        }

        let mut per_bin: HashMap<u64, (Moments, Moments)> = HashMap::new();
        for (bin, m) in &out.a_moments {
            per_bin.entry(*bin).or_default().0 = *m;
        }
        for (bin, m) in &out.b_moments {
            per_bin.entry(*bin).or_default().1 = *m;
        }

        track.lock().unwrap().insert(chrom.clone(), smoothed);
        moms.lock().unwrap().insert(chrom.clone(), per_bin);
        cov.lock().unwrap().insert(chrom.clone(), out.covered);
    });

    (
        track.into_inner().unwrap(),
        moms.into_inner().unwrap(),
        cov.into_inner().unwrap(),
    )
}

/// Drop bins with no coverage at `res`.
///
/// Filtering inside the shifter is not sufficient on its own: `merge_tracks`
/// expands each resolution onto the finest grid, so a 2.5 Mb bin spanning a
/// centromere repopulates every fine slot beneath it from a pass where that
/// span *did* have coverage. The merged track has to be filtered as well, at
/// the resolution it is expressed in.
fn filter_by_coverage(track: &Track, masks: &HashMap<String, Vec<bool>>, res: u64) -> Track {
    let mut out: Track = HashMap::new();
    for (chrom, rows) in track {
        let Some(mask) = masks.get(chrom) else {
            out.insert(chrom.clone(), rows.clone());
            continue;
        };
        let kept: Vec<(u64, u64, f64)> = rows
            .iter()
            .filter(|(s, _, _)| mask.get((s / res) as usize).copied() == Some(true))
            .copied()
            .collect();
        if !kept.is_empty() {
            out.insert(chrom.clone(), kept);
        }
    }
    out
}

/// Re-centre `fine` with a rolling mean, optionally against the coarser track.
fn apply_shifter(
    cfg: &Config,
    sizes: &ChromSizes,
    fine: &Track,
    coarse: Option<&Track>,
    res: u64,
    coarse_res: u64,
    masks: Option<&HashMap<String, Vec<bool>>>,
) -> Track {
    let w = (3 * (coarse_res / res.max(1))).max(1) as usize;
    let half = w / 2;

    let mut out: Track = HashMap::new();
    for (chrom, rows) in fine {
        let mask = masks.and_then(|m| m.get(chrom));
        let size = sizes.length_of(chrom).unwrap_or(0);
        // Match the bash shifter's allocation, which pads the grid by two
        // coarse bins. Those padding bins are emitted, enter the merge, and are
        // only clipped to the chromosome at the very end -- so omitting them
        // changes the merged bin set and hence the GIave normalisation.
        let n = ((size + 2 * coarse_res) / res) as usize + 1;

        let mut vals = vec![0f64; n];
        for &(s, _, v) in rows {
            let i = (s / res) as usize;
            if i < n {
                vals[i] = v;
            }
        }

        let mut cvals = vec![0f64; n];
        if cfg.shifter_mode == ShifterMode::Coarse {
            if let Some(ct) = coarse.and_then(|c| c.get(chrom)) {
                if cfg.coarse_interp {
                    // Interpolate the coarse score between coarse-bin CENTRES.
                    //
                    // Scattering it as a step function puts a discontinuity at
                    // every coarse bin edge, and because `shifted = vals -
                    // (rolling - cvals)` that jump lands in the output at full
                    // magnitude. It then compounds: each pass's shifted output
                    // becomes the next pass's coarse track, so a step created at
                    // 1 Mb propagates down to 1 kb undiminished. Measured on the
                    // LCL 1 kb track, adjacent-bin sign flips ran 38.8x the
                    // background rate exactly at 1 Mb boundaries (31.2% of pairs)
                    // and 18.1x at 100 kb, while a single-resolution run -- which
                    // has no coarse track at all -- shows no excess (0.97x).
                    //
                    // The regional signal the shifter exists to restore is
                    // low-frequency, so interpolating loses nothing real and
                    // removes the artefact by construction.
                    let mut pts: Vec<(f64, f64)> = ct
                        .iter()
                        .filter(|(s, e, _)| e > s)
                        .map(|&(s, e, v)| (((s + e) as f64) / 2.0 / res as f64, v))
                        .collect();
                    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                    if !pts.is_empty() {
                        let mut k = 0usize;
                        for (i, slot) in cvals.iter_mut().enumerate() {
                            let x = i as f64;
                            while k + 1 < pts.len() && pts[k + 1].0 < x {
                                k += 1;
                            }
                            *slot = if x <= pts[0].0 {
                                pts[0].1
                            } else if x >= pts[pts.len() - 1].0 {
                                pts[pts.len() - 1].1
                            } else {
                                let (x0, y0) = pts[k];
                                let (x1, y1) = pts[(k + 1).min(pts.len() - 1)];
                                if (x1 - x0).abs() < 1e-9 {
                                    y0
                                } else {
                                    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
                                }
                            };
                        }
                    }
                } else {
                    for &(s, e, v) in ct {
                        let a = (s / res) as usize;
                        let b = ((e / res) as usize).min(n);
                        for slot in cvals.iter_mut().take(b).skip(a) {
                            *slot = v;
                        }
                    }
                }
            }
        }

        // O(n) centred rolling mean via prefix sums.
        let mut prefix = vec![0f64; n + 1];
        for i in 0..n {
            prefix[i + 1] = prefix[i] + vals[i];
        }

        let mut rows_out = Vec::with_capacity(n);
        for i in 0..n {
            // Bins the map reports as uncovered carry no information; emitting
            // them would invent a score out of the neighbouring rolling mean.
            if let Some(m) = mask {
                if m.get(i).copied() != Some(true) {
                    continue;
                }
            }
            let lo = i.saturating_sub(half);
            let hi = (i + half).min(n);
            let cnt = (hi - lo).max(1) as f64;
            let rolling = (prefix[hi] - prefix[lo]) / cnt;
            let start = (i as u64) * res;
            rows_out.push((start, start + res, vals[i] - (rolling - cvals[i])));
        }
        out.insert(chrom.clone(), rows_out);
    }
    out
}

/// Derive A/B seeds at resolution `res` from a compartment track.
///
/// Every pass of round k is initialised from the merged track round k-1
/// produced; the coarsest pass of round 0 is initialised from the best
/// compartment-matching principal component.
///
/// Assignment is by **majority**: the values of all track bins falling inside a
/// bin at `res` are summed, and the bin is seeded A if the sum is positive, B
/// if negative. The two seed sets are therefore disjoint at every resolution.
/// Assigning per-bin instead (A if any positive bin, B if any negative) puts
/// nearly every coarse bin into both sets, because fine-scale calls alternate
/// within a coarse bin — which collapses the A-versus-B contrast the score
/// depends on.
///
/// The bash v1 implementation re-intersected the *original* annotation at every
/// step and unioned the result with its own complement, which returns the
/// original set for any input, so the seeds never changed at all (bug B23).
///
/// ## Degeneracy guard
///
/// Splitting on sign alone is a positive feedback loop. If a chromosome drifts
/// so that every bin is positive it yields *no B seeds*, the next round scores
/// it against an all-A seed set, and the drift amplifies. Measured on the LCL
/// map at 100 kb with per-chromosome mean-centring disabled: chr13 reached 100%
/// positive bins, chr9 98%, chrX 91%, chr21 3%, with the score SD collapsing to
/// 26-50 against ~120 on healthy chromosomes.
///
/// Mean-centring masked this by forcing the split back to ~50/50 every round,
/// which is why `endZ` was effectively always-on. That is a fix to the symptom:
/// it subtracts a constant from every reported score to keep the *seeding*
/// stable.
///
/// The guard here targets the loop instead. A chromosome whose seed split falls
/// below `MIN_SEED_FRACTION` in either class is re-split at its median, which
/// cannot degenerate. Chromosomes with a genuine A/B imbalance keep it — the
/// guard only engages where the split would otherwise collapse.
/// Bin a track onto a `res` grid, summing the values that touch each bin.
fn track_to_net(
    track: &Track,
    res: u64,
    anchor: Option<&HashMap<String, std::collections::HashSet<u64>>>,
) -> HashMap<String, HashMap<u64, f64>> {
    let mut out: HashMap<String, HashMap<u64, f64>> = HashMap::new();
    for (chrom, rows) in track {
        let mut net: HashMap<u64, f64> = HashMap::new();
        for &(s, e, v) in rows {
            let lo = (s / res) * res;
            let hi = if e > s { ((e - 1) / res) * res } else { lo };
            let mut bin = lo;
            loop {
                *net.entry(bin).or_insert(0.0) += v;
                if bin >= hi {
                    break;
                }
                bin += res;
            }
        }
        if let Some(all) = anchor {
            let keep = all.get(chrom);
            net.retain(|b, _| keep.map_or(false, |k| k.contains(b)));
        }
        out.insert(chrom.clone(), net);
    }
    out
}

/// Blend this round's binned track into the carried-forward seed state.
///
/// `state <- (1 - d) * fresh + d * state`, per bin. The walk still follows the
/// data, but a difference introduced in one round can only move the seeds by
/// (1-d) of its size, so it compounds geometrically more slowly across rounds
/// instead of being handed on at full weight.
fn damp_net(
    state: &mut HashMap<String, HashMap<u64, f64>>,
    fresh: HashMap<String, HashMap<u64, f64>>,
    d: f64,
) {
    for (chrom, net) in fresh {
        let e = state.entry(chrom).or_default();
        for (bin, v) in net {
            let prev = e.get(&bin).copied();
            let blended = match prev {
                Some(p) => (1.0 - d) * v + d * p,
                None => v,
            };
            e.insert(bin, blended);
        }
    }
}

fn seeds_from_net(
    per_chrom_in: &HashMap<String, HashMap<u64, f64>>,
    res: u64,
    median_split: bool,
    seed_frac: f64,
) -> (Seeds, Seeds) {
    // Below this fraction in either class, re-split that chromosome at its own
    // median: a last-resort floor against a chromosome producing no seeds of
    // one class at all.
    const MIN_SEED_FRACTION: f64 = 0.05;

    let mut a: Seeds = HashMap::new();
    let mut b: Seeds = HashMap::new();
    let per_chrom: Vec<(String, HashMap<u64, f64>)> =
        per_chrom_in.iter().map(|(c, n)| (c.clone(), n.clone())).collect();

    // Split each chromosome at its own median when `median_split`, otherwise at
    // zero (the historical behaviour).
    //
    // A genome-wide threshold was tried and is worse: it reintroduces
    // saturation (1 chromosome, min SD 31.8 against 80.9) and strengthens the
    // gene-density anti-correlation (-0.65 against -0.57). The per-chromosome
    // median is what actually removes the runaway.
    for (chrom, net) in per_chrom {
        let threshold = if median_split && !net.is_empty() {
            let mut vals: Vec<f64> = net.values().copied().collect();
            vals.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
            vals[vals.len() / 2]
        } else {
            0.0
        };
        let mut av: Vec<u64> = Vec::new();
        let mut bv: Vec<u64> = Vec::new();
        for (&bin, &sum) in &net {
            if sum > threshold {
                av.push(bin);
            } else if sum < threshold {
                bv.push(bin);
            }
        }

        // Floor: a chromosome that still collapses is re-split at its own
        // median, which cannot degenerate.
        let total = (av.len() + bv.len()) as f64;
        if total > 0.0 {
            let frac = av.len().min(bv.len()) as f64 / total;
            if frac < MIN_SEED_FRACTION {
                let mut vals: Vec<(u64, f64)> = net.iter().map(|(&k, &v)| (k, v)).collect();
                vals.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal));
                let mid = vals.len() / 2;
                bv = vals[..mid].iter().map(|p| p.0).collect();
                av = vals[mid..].iter().map(|p| p.0).collect();
                eprintln!(
                    "      NOTE: {chrom} seed split degenerate at {res} bp \
                     ({:.1}% minority) — re-split at its own median",
                    frac * 100.0
                );
            }
        }

        // Keep only the most confident seeds: the bins furthest from the
        // split point. Bins sitting on the threshold are close to a coin flip,
        // and because round k+1 seeds from round k's track, that indecision is
        // fed forward and amplified. Measured on Wi38, replicate-to-replicate
        // disagreement barely falls from 1 kb to 100 kb (28.7 -> 22.1) despite
        // 100x more reads per bin, so the dominant noise is this feedback and
        // not counting statistics.
        if seed_frac < 1.0 {
            let keep = |v: &mut Vec<u64>, desc: bool| {
                let mut with: Vec<(u64, f64)> = v
                    .iter()
                    .map(|&b| (b, (net.get(&b).copied().unwrap_or(0.0) - threshold).abs()))
                    .collect();
                with.sort_by(|x, y| y.1.partial_cmp(&x.1).unwrap_or(std::cmp::Ordering::Equal));
                let n = ((with.len() as f64) * seed_frac).ceil().max(1.0) as usize;
                let _ = desc;
                *v = with.into_iter().take(n).map(|p| p.0).collect();
            };
            keep(&mut av, true);
            keep(&mut bv, false);
        }

        av.sort_unstable();
        bv.sort_unstable();
        if !av.is_empty() {
            a.insert(chrom.clone(), (res, av));
        }
        if !bv.is_empty() {
            b.insert(chrom.clone(), (res, bv));
        }
    }
    (a, b)
}

/// Convert a bedGraph to bigWig with `bedGraphToBigWig`, returning the bigWig
/// path on success.
///
/// bigWig is what genome browsers actually want, and a 1 kb genome-wide track
/// is ~90 MB as text against ~27 MB as bigWig with random access. The
/// conversion is shelled out rather than implemented here: writing the format
/// correctly is a lot of surface area for something every genomics box already
/// has from UCSC's toolset.
///
/// Failure is never fatal — the bedGraph is already on disk and is the real
/// output. A missing converter downgrades the run, it does not break it.
fn bedgraph_to_bigwig(
    bedgraph: &PathBuf,
    sizes: &ChromSizes,
    keep_bedgraph: bool,
    chr_prefix: bool,
) -> Option<PathBuf> {
    use std::process::Command;

    // --chr-prefix rewrites bare names ("1") to UCSC style ("chr1") in the
    // bigWig only. Off by default: this dataset already carries two
    // assembly-naming traps, and a tool that silently renames chromosomes
    // would add a third.
    let pfx = |n: &str| -> String {
        if chr_prefix && !n.starts_with("chr") {
            format!("chr{n}")
        } else {
            n.to_string()
        }
    };

    if Command::new("bedGraphToBigWig").arg("-version").output().is_err() {
        eprintln!(
            "      NOTE: bedGraphToBigWig not on PATH — leaving bedGraph only.\n\
             \x20           Install UCSC's bedGraphToBigWig, or pass --bigwig 0 to silence this."
        );
        return None;
    }

    let dir = bedgraph.parent()?.to_path_buf();
    let sizes_path = dir.join(".crush_chrom.sizes");
    let sorted = dir.join(".crush_sorted.bedgraph");
    {
        let f = File::create(&sizes_path).ok()?;
        let mut w = BufWriter::new(f);
        for name in &sizes.names {
            if let Some(len) = sizes.length_of(name) {
                writeln!(w, "{}\t{len}", pfx(name)).ok()?;
            }
        }
        w.flush().ok()?;
    }

    let ok = Command::new("sh")
        .arg("-c")
        .arg(if chr_prefix {
            format!(
                "grep -v '^track' {bg} | awk 'BEGIN{{OFS=\"\\t\"}} \
                 {{print ($1 ~ /^chr/ ? $1 : \"chr\"$1), $2, $3, $4}}' \
                 | LC_ALL=C sort -k1,1 -k2,2n > {out}",
                bg = bedgraph.display(),
                out = sorted.display()
            )
        } else {
            format!(
                "grep -v '^track' {bg} | LC_ALL=C sort -k1,1 -k2,2n > {out}",
                bg = bedgraph.display(),
                out = sorted.display()
            )
        })
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !ok {
        let _ = std::fs::remove_file(&sizes_path);
        return None;
    }

    let bw = bedgraph.with_extension("bw");
    let status = Command::new("bedGraphToBigWig")
        .arg(&sorted)
        .arg(&sizes_path)
        .arg(&bw)
        .status();
    let _ = std::fs::remove_file(&sorted);
    let _ = std::fs::remove_file(&sizes_path);

    match status {
        Ok(s) if s.success() => {
            if !keep_bedgraph {
                let _ = std::fs::remove_file(bedgraph);
            }
            Some(bw)
        }
        _ => {
            eprintln!("      NOTE: bedGraphToBigWig failed — leaving bedGraph only.");
            None
        }
    }
}

fn compute_pvalues(moments: &MomentMap, res: u64) -> Track {
    let mut out: Track = HashMap::new();
    for (chrom, bins) in moments {
        let mut rows: Vec<(u64, u64, f64)> = bins
            .iter()
            .filter_map(|(&bin, (a, b))| stats::welch_p(a, b).map(|p| (bin, bin + res, p)))
            .collect();
        rows.sort_unstable_by_key(|r| r.0);
        if !rows.is_empty() {
            out.insert(chrom.clone(), rows);
        }
    }
    out
}

/// Expand every resolution onto the finest grid and average per bin.
///
/// Values are divided by the interval width before expansion, so coarse
/// resolutions contribute proportionally less — the resolution weighting the
/// bash merge applies.
fn merge_tracks(tracks: &[(u64, Track)], finest: u64) -> Track {
    let mut acc: HashMap<String, HashMap<u64, (f64, u32)>> = HashMap::new();

    for (_res, track) in tracks {
        for (chrom, rows) in track {
            let e = acc.entry(chrom.clone()).or_default();
            for &(s, end, v) in rows {
                if end <= s {
                    continue;
                }
                let width = end - s;
                let contrib = v / width as f64;
                let mut b = (s / finest) * finest;
                while b < end {
                    let slot = e.entry(b).or_insert((0.0, 0));
                    slot.0 += contrib;
                    slot.1 += 1;
                    b += finest;
                }
            }
        }
    }

    let mut out: Track = HashMap::new();
    for (chrom, bins) in acc {
        let mut rows: Vec<(u64, u64, f64)> = bins
            .into_iter()
            .map(|(b, (sum, n))| (b, b + finest, sum / n as f64))
            .collect();
        rows.sort_unstable_by_key(|r| r.0);
        out.insert(chrom, rows);
    }
    out
}

fn write_track(
    path: &PathBuf,
    track: &Track,
    sizes: &ChromSizes,
    scale: Option<f64>,
    clamp: Option<f64>,
) -> Result<(), String> {
    let f = File::create(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut w = BufWriter::with_capacity(1 << 20, f);

    writeln!(
        w,
        "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 \
         viewLimits=-100:100 autoScale off"
    )
    .map_err(|e| e.to_string())?;

    // Emit in chromosome-sizes order so output is stable and diffable.
    for chrom in &sizes.names {
        let Some(rows) = track.get(chrom) else {
            continue;
        };
        let limit = sizes.length_of(chrom).unwrap_or(u64::MAX);
        for &(s, e, v) in rows {
            if e > limit {
                continue;
            }
            let val = match scale {
                Some(k) if k != 0.0 => v / k,
                _ => v,
            };
            if !val.is_finite() {
                continue;
            }
            // Clamp rather than drop: an extreme bin is still a real call, and
            // removing it would leave a hole in the track. Only its magnitude
            // is untrustworthy, and only the magnitude is what breaks autoscale.
            let val = match clamp {
                Some(c) if c > 0.0 => val.clamp(-c, c),
                _ => val,
            };
            writeln!(w, "{chrom}\t{s}\t{e}\t{val:.6}").map_err(|e| e.to_string())?;
        }
    }
    w.flush().map_err(|e| e.to_string())?;
    Ok(())
}
