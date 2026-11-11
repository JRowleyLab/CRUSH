//! Native reader for Juicer `.hic` files (versions 6-9).
//!
//! This replaces the `hicstraw` Python round-trip in the bash pipeline. The
//! format is not formally specified anywhere; this implementation follows the
//! reference `straw` C++ reader:
//!   - header: magic, version, footer position, genome id, attributes,
//!     chromosome table, base-pair resolutions, fragment resolutions
//!   - footer: master index (byte offsets of each chromosome-pair matrix),
//!     expected-value vectors, normalised expected-value vectors, and the
//!     index of per-chromosome normalisation vectors
//!   - matrix: a list of resolution levels ("zoom data"), each carrying a
//!     block grid; blocks are individually zlib-deflated
//!
//! v9 differs from v8 in three ways that matter here: block coordinates may be
//! stored along the diagonal rather than in a rectangular grid, counts may be
//! stored as short/int as well as float, and the normalisation index carries
//! 64-bit sizes.

use flate2::read::ZlibDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    NotHic,
    UnsupportedVersion(i32),
    UnknownChrom(String),
    NoResolution(u32),
    NoNorm(String),
    Corrupt(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error reading .hic file: {e}"),
            Error::NotHic => write!(f, "not a .hic file (bad magic bytes)"),
            Error::UnsupportedVersion(v) => {
                write!(f, ".hic version {v} is not supported (need 6-9)")
            }
            Error::UnknownChrom(c) => write!(f, "chromosome '{c}' not present in the .hic file"),
            Error::NoResolution(r) => write!(f, "resolution {r} bp not stored in the .hic file"),
            Error::NoNorm(n) => write!(f, "normalisation '{n}' not available for this record"),
            Error::Corrupt(m) => write!(f, "corrupt .hic structure: {m}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Error::Io(e)
    }
}

// ---------------------------------------------------------------------------
// Little-endian primitive reads
// ---------------------------------------------------------------------------

trait ReadLe: Read {
    fn u8v(&mut self) -> io::Result<u8> {
        let mut b = [0u8; 1];
        self.read_exact(&mut b)?;
        Ok(b[0])
    }
    fn i16v(&mut self) -> io::Result<i16> {
        let mut b = [0u8; 2];
        self.read_exact(&mut b)?;
        Ok(i16::from_le_bytes(b))
    }
    fn i32v(&mut self) -> io::Result<i32> {
        let mut b = [0u8; 4];
        self.read_exact(&mut b)?;
        Ok(i32::from_le_bytes(b))
    }
    fn i64v(&mut self) -> io::Result<i64> {
        let mut b = [0u8; 8];
        self.read_exact(&mut b)?;
        Ok(i64::from_le_bytes(b))
    }
    fn f32v(&mut self) -> io::Result<f32> {
        let mut b = [0u8; 4];
        self.read_exact(&mut b)?;
        Ok(f32::from_le_bytes(b))
    }
    fn f64v(&mut self) -> io::Result<f64> {
        let mut b = [0u8; 8];
        self.read_exact(&mut b)?;
        Ok(f64::from_le_bytes(b))
    }
    /// NUL-terminated string.
    fn cstr(&mut self) -> io::Result<String> {
        let mut out = Vec::new();
        loop {
            let c = self.u8v()?;
            if c == 0 {
                break;
            }
            out.push(c);
        }
        Ok(String::from_utf8_lossy(&out).into_owned())
    }
}

impl<R: Read> ReadLe for R {}

// ---------------------------------------------------------------------------
// Public data types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Chromosome {
    pub name: String,
    pub length: i64,
    pub index: usize,
}

/// One contact record: bin start coordinates in bp, plus the value.
#[derive(Debug, Clone, Copy)]
pub struct Contact {
    pub x: u32,
    pub y: u32,
    pub value: f32,
}

/// Which value to report for each contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValueKind {
    /// Raw (optionally normalised) counts.
    Observed,
    /// Counts divided by the genome-wide expected value at that distance.
    Oe,
}

#[derive(Debug, Clone)]
struct BlockIndexEntry {
    position: i64,
    size: i32,
}

#[derive(Debug, Clone)]
struct ZoomData {
    bin_size: i32,
    sum_counts: f32,
    block_bin_count: i32,
    block_column_count: i32,
    blocks: HashMap<i32, BlockIndexEntry>,
}

#[derive(Debug, Clone)]
struct ExpectedValues {
    bin_size: i32,
    values: Vec<f64>,
    /// chromosome index -> scaling factor
    chr_scale: HashMap<i32, f64>,
}

#[derive(Debug, Clone)]
struct NormVectorIndex {
    position: i64,
    n_bytes: i64,
}

pub struct HicFile {
    reader: BufReader<File>,
    pub version: i32,
    pub genome_id: String,
    pub chromosomes: Vec<Chromosome>,
    chrom_by_name: HashMap<String, usize>,
    pub bp_resolutions: Vec<u32>,
    /// (chr1_index, chr2_index) -> byte offset of that matrix record
    master_index: HashMap<String, i64>,
    /// unit -> list of expected-value vectors, one per resolution
    expected: HashMap<String, Vec<ExpectedValues>>,
    /// (norm, unit) -> list of expected-value vectors
    norm_expected: HashMap<(String, String), Vec<ExpectedValues>>,
    /// (norm, chr_index, unit, bin_size) -> where the vector lives
    norm_index: HashMap<(String, i32, String, i32), NormVectorIndex>,
}

impl HicFile {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::with_capacity(1 << 20, file);

        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic[..3] != b"HIC" || magic[3] != 0 {
            return Err(Error::NotHic);
        }

        let version = reader.i32v()?;
        if !(6..=9).contains(&version) {
            return Err(Error::UnsupportedVersion(version));
        }

        let footer_position = reader.i64v()?;
        let genome_id = reader.cstr()?;

        // v9 added a normalised-expected-vector index pointer here.
        let _nvi_position = if version >= 9 {
            let p = reader.i64v()?;
            let _len = reader.i64v()?;
            Some(p)
        } else {
            None
        };

        // Attribute dictionary (ignored, but must be consumed to stay in sync).
        let n_attributes = reader.i32v()?;
        for _ in 0..n_attributes {
            let _k = reader.cstr()?;
            let _v = reader.cstr()?;
        }

        let n_chroms = reader.i32v()?;
        let mut chromosomes = Vec::with_capacity(n_chroms as usize);
        let mut chrom_by_name = HashMap::new();
        for i in 0..n_chroms {
            let name = reader.cstr()?;
            // v9 widened chromosome lengths to 64 bits.
            let length = if version >= 9 {
                reader.i64v()?
            } else {
                reader.i32v()? as i64
            };
            chrom_by_name.insert(name.to_lowercase(), i as usize);
            chromosomes.push(Chromosome {
                name,
                length,
                index: i as usize,
            });
        }

        let n_bp_res = reader.i32v()?;
        let mut bp_resolutions = Vec::with_capacity(n_bp_res as usize);
        for _ in 0..n_bp_res {
            bp_resolutions.push(reader.i32v()? as u32);
        }
        // Fragment resolutions are read but unused: CRUSH is BP-only.
        let n_frag_res = reader.i32v()?;
        for _ in 0..n_frag_res {
            let _ = reader.i32v()?;
        }

        let mut hic = HicFile {
            reader,
            version,
            genome_id,
            chromosomes,
            chrom_by_name,
            bp_resolutions,
            master_index: HashMap::new(),
            expected: HashMap::new(),
            norm_expected: HashMap::new(),
            norm_index: HashMap::new(),
        };
        hic.read_footer(footer_position)?;
        Ok(hic)
    }

    fn read_footer(&mut self, position: i64) -> Result<()> {
        self.reader.seek(SeekFrom::Start(position as u64))?;

        // Total footer size: i32 in v8 and earlier, i64 from v9.
        if self.version >= 9 {
            let _ = self.reader.i64v()?;
        } else {
            let _ = self.reader.i32v()?;
        }

        let n_entries = self.reader.i32v()?;
        for _ in 0..n_entries {
            let key = self.reader.cstr()?;
            let pos = self.reader.i64v()?;
            let _size = self.reader.i32v()?;
            self.master_index.insert(key, pos);
        }

        // Expected value vectors (no normalisation).
        let n_ev = self.reader.i32v()?;
        for _ in 0..n_ev {
            let ev = self.read_expected_vector()?;
            self.expected.entry(ev.0).or_default().push(ev.1);
        }

        // Normalised expected value vectors. Older files may simply end here.
        let n_nev = match self.reader.i32v() {
            Ok(v) => v,
            Err(_) => return Ok(()),
        };
        for _ in 0..n_nev {
            let norm = self.reader.cstr()?;
            let ev = self.read_expected_vector()?;
            self.norm_expected
                .entry((norm, ev.0))
                .or_default()
                .push(ev.1);
        }

        // Index of per-chromosome normalisation vectors.
        let n_norm = match self.reader.i32v() {
            Ok(v) => v,
            Err(_) => return Ok(()),
        };
        for _ in 0..n_norm {
            let norm = self.reader.cstr()?;
            let chr_idx = self.reader.i32v()?;
            let unit = self.reader.cstr()?;
            let bin_size = self.reader.i32v()?;
            let pos = self.reader.i64v()?;
            let n_bytes = if self.version >= 9 {
                self.reader.i64v()?
            } else {
                self.reader.i32v()? as i64
            };
            self.norm_index
                .insert((norm, chr_idx, unit, bin_size), NormVectorIndex {
                    position: pos,
                    n_bytes,
                });
        }

        Ok(())
    }

    fn read_expected_vector(&mut self) -> Result<(String, ExpectedValues)> {
        let unit = self.reader.cstr()?;
        let bin_size = self.reader.i32v()?;

        // v9 stores the vector as f32 with an i64 length; v8 as f64 with i32.
        let n_values = if self.version >= 9 {
            self.reader.i64v()?
        } else {
            self.reader.i32v()? as i64
        };
        let mut values = Vec::with_capacity(n_values as usize);
        for _ in 0..n_values {
            values.push(if self.version >= 9 {
                self.reader.f32v()? as f64
            } else {
                self.reader.f64v()?
            });
        }

        let n_scale = if self.version >= 9 {
            self.reader.i32v()?
        } else {
            self.reader.i32v()?
        };
        let mut chr_scale = HashMap::new();
        for _ in 0..n_scale {
            let chr_idx = self.reader.i32v()?;
            let factor = if self.version >= 9 {
                self.reader.f32v()? as f64
            } else {
                self.reader.f64v()?
            };
            chr_scale.insert(chr_idx, factor);
        }

        Ok((
            unit,
            ExpectedValues {
                bin_size,
                values,
                chr_scale,
            },
        ))
    }

    pub fn chromosome(&self, name: &str) -> Option<&Chromosome> {
        // .hic files are inconsistent about the chr prefix; match either way.
        let lower = name.to_lowercase();
        if let Some(&i) = self.chrom_by_name.get(&lower) {
            return Some(&self.chromosomes[i]);
        }
        let alt = if let Some(stripped) = lower.strip_prefix("chr") {
            stripped.to_string()
        } else {
            format!("chr{lower}")
        };
        self.chrom_by_name.get(&alt).map(|&i| &self.chromosomes[i])
    }

    /// Read all intrachromosomal contacts for `chrom` at `resolution`.
    ///
    /// This is the Rust equivalent of
    /// `hicstraw.straw(kind, norm, file, chrom, chrom, "BP", res)`.
    pub fn straw_cis(
        &mut self,
        kind: ValueKind,
        norm: &str,
        chrom: &str,
        resolution: u32,
    ) -> Result<Vec<Contact>> {
        let chr = self
            .chromosome(chrom)
            .ok_or_else(|| Error::UnknownChrom(chrom.to_string()))?
            .clone();

        if !self.bp_resolutions.contains(&resolution) {
            return Err(Error::NoResolution(resolution));
        }

        let key = format!("{}_{}", chr.index, chr.index);
        let matrix_pos = *self
            .master_index
            .get(&key)
            .ok_or_else(|| Error::Corrupt(format!("no matrix record for {}", chr.name)))?;

        let zoom = self.read_matrix_zoom(matrix_pos, resolution)?;

        // Normalisation vector, if requested.
        let norm_vec = if norm.eq_ignore_ascii_case("NONE") {
            None
        } else {
            Some(self.read_norm_vector(norm, chr.index as i32, resolution)?)
        };

        // Expected-value vector, if we are reporting O/E.
        let expected = if kind == ValueKind::Oe {
            Some(self.expected_vector(norm, resolution, chr.index as i32)?)
        } else {
            None
        };

        let mut positions: Vec<(i64, i32)> = zoom
            .blocks
            .values()
            .map(|b| (b.position, b.size))
            .collect();
        // Sequential offsets keep the BufReader from thrashing.
        positions.sort_unstable();

        let mut out = Vec::new();
        for (pos, size) in positions {
            let records = self.read_block(pos, size)?;
            for (bx, by, v) in records {
                let mut value = v;
                if let Some(nv) = &norm_vec {
                    let (a, b) = (bx as usize, by as usize);
                    if a >= nv.len() || b >= nv.len() {
                        continue;
                    }
                    let (na, nb) = (nv[a], nv[b]);
                    if !(na.is_finite() && nb.is_finite()) || na == 0.0 || nb == 0.0 {
                        continue;
                    }
                    value /= (na * nb) as f32;
                }
                if let Some(ev) = &expected {
                    let d = (bx - by).unsigned_abs() as usize;
                    let e = if d < ev.len() { ev[d] } else { 0.0 };
                    if e == 0.0 || !e.is_finite() {
                        continue;
                    }
                    value /= e as f32;
                }
                if !value.is_finite() {
                    continue;
                }
                out.push(Contact {
                    x: (bx as u32).saturating_mul(resolution),
                    y: (by as u32).saturating_mul(resolution),
                    value,
                });
            }
        }

        Ok(out)
    }

    fn read_matrix_zoom(&mut self, position: i64, resolution: u32) -> Result<ZoomData> {
        self.reader.seek(SeekFrom::Start(position as u64))?;
        let _chr1 = self.reader.i32v()?;
        let _chr2 = self.reader.i32v()?;
        let n_resolutions = self.reader.i32v()?;

        for _ in 0..n_resolutions {
            let unit = self.reader.cstr()?;
            let _zoom_index = self.reader.i32v()?;
            let _sum_counts = self.reader.f32v()?;
            let _occupied = self.reader.f32v()?;
            let _std_dev = self.reader.f32v()?;
            let _percent95 = self.reader.f32v()?;
            let bin_size = self.reader.i32v()?;
            let block_bin_count = self.reader.i32v()?;
            let block_column_count = self.reader.i32v()?;
            let block_count = self.reader.i32v()?;

            let wanted = unit == "BP" && bin_size == resolution as i32;

            let mut blocks = HashMap::new();
            for _ in 0..block_count {
                let block_number = self.reader.i32v()?;
                let block_position = self.reader.i64v()?;
                let block_size = self.reader.i32v()?;
                if wanted {
                    blocks.insert(
                        block_number,
                        BlockIndexEntry {
                            position: block_position,
                            size: block_size,
                        },
                    );
                }
            }

            if wanted {
                return Ok(ZoomData {
                    bin_size,
                    sum_counts: _sum_counts,
                    block_bin_count,
                    block_column_count,
                    blocks,
                });
            }
        }

        Err(Error::NoResolution(resolution))
    }

    /// Decompress and parse a single block into (binX, binY, value) triples.
    fn read_block(&mut self, position: i64, size: i32) -> Result<Vec<(i32, i32, f32)>> {
        if size <= 0 {
            return Ok(Vec::new());
        }
        self.reader.seek(SeekFrom::Start(position as u64))?;
        let mut compressed = vec![0u8; size as usize];
        self.reader.read_exact(&mut compressed)?;

        let mut raw = Vec::new();
        ZlibDecoder::new(&compressed[..]).read_to_end(&mut raw)?;
        let mut cur = io::Cursor::new(raw);

        let n_records = cur.i32v()?;
        let mut out = Vec::with_capacity(n_records.max(0) as usize);

        if self.version < 7 {
            // v6: a flat list of (binX, binY, counts) as i32/i32/f32.
            for _ in 0..n_records {
                let x = cur.i32v()?;
                let y = cur.i32v()?;
                let v = cur.f32v()?;
                out.push((x, y, v));
            }
            return Ok(out);
        }

        let bin_x_offset = cur.i32v()?;
        let bin_y_offset = cur.i32v()?;

        // v9 allows short as well as float counts; v7/v8 are float-only.
        let use_short = if self.version >= 9 {
            cur.u8v()? == 0
        } else {
            cur.u8v()? == 0
        };
        // v9 additionally allows short bin coordinates.
        let (use_short_bin_x, use_short_bin_y) = if self.version >= 9 {
            (cur.u8v()? == 0, cur.u8v()? == 0)
        } else {
            (true, true)
        };

        let block_type = cur.u8v()?;

        let mut read_count = |c: &mut io::Cursor<Vec<u8>>| -> io::Result<f32> {
            if use_short {
                let v = c.i16v()?;
                // -32768 is the sentinel for "no value" in v9 short blocks.
                Ok(if v == i16::MIN { f32::NAN } else { v as f32 })
            } else {
                c.f32v()
            }
        };

        match block_type {
            1 => {
                // Row-major list of runs.
                let row_count = if use_short_bin_y {
                    cur.i16v()? as i32
                } else {
                    cur.i32v()?
                };
                for _ in 0..row_count {
                    let y_off = if use_short_bin_y {
                        cur.i16v()? as i32
                    } else {
                        cur.i32v()?
                    };
                    let bin_y = y_off + bin_y_offset;
                    let col_count = if use_short_bin_x {
                        cur.i16v()? as i32
                    } else {
                        cur.i32v()?
                    };
                    for _ in 0..col_count {
                        let x_off = if use_short_bin_x {
                            cur.i16v()? as i32
                        } else {
                            cur.i32v()?
                        };
                        let bin_x = x_off + bin_x_offset;
                        let v = read_count(&mut cur)?;
                        if v.is_finite() {
                            out.push((bin_x, bin_y, v));
                        }
                    }
                }
            }
            2 => {
                // Dense representation over a square patch.
                let n_pts = cur.i32v()?;
                let w = cur.i16v()? as i32;
                for i in 0..n_pts {
                    let row = i / w;
                    let col = i - row * w;
                    let bin_x = bin_x_offset + col;
                    let bin_y = bin_y_offset + row;
                    let v = read_count(&mut cur)?;
                    if v.is_finite() {
                        out.push((bin_x, bin_y, v));
                    }
                }
            }
            _ => {
                return Err(Error::Corrupt(format!(
                    "unknown block type {block_type} at offset {position}"
                )))
            }
        }

        Ok(out)
    }

    fn read_norm_vector(&mut self, norm: &str, chr_index: i32, resolution: u32) -> Result<Vec<f64>> {
        let key = (
            norm.to_string(),
            chr_index,
            "BP".to_string(),
            resolution as i32,
        );
        let entry = self
            .norm_index
            .get(&key)
            .cloned()
            .ok_or_else(|| Error::NoNorm(norm.to_string()))?;

        self.reader.seek(SeekFrom::Start(entry.position as u64))?;
        let n = if self.version >= 9 {
            self.reader.i64v()?
        } else {
            self.reader.i32v()? as i64
        };
        let mut v = Vec::with_capacity(n as usize);
        for _ in 0..n {
            v.push(if self.version >= 9 {
                self.reader.f32v()? as f64
            } else {
                self.reader.f64v()?
            });
        }
        Ok(v)
    }

    /// Expected-value vector for this (norm, resolution, chromosome), already
    /// multiplied by the per-chromosome scaling factor.
    fn expected_vector(&self, norm: &str, resolution: u32, chr_index: i32) -> Result<Vec<f64>> {
        let list = if norm.eq_ignore_ascii_case("NONE") {
            self.expected.get("BP")
        } else {
            self.norm_expected
                .get(&(norm.to_string(), "BP".to_string()))
        };

        let list = list.ok_or_else(|| {
            Error::Corrupt("no expected-value vectors stored in this .hic file".into())
        })?;

        let ev = list
            .iter()
            .find(|e| e.bin_size == resolution as i32)
            .ok_or(Error::NoResolution(resolution))?;

        let scale = ev.chr_scale.get(&chr_index).copied().unwrap_or(1.0);
        if scale == 0.0 || !scale.is_finite() {
            return Ok(ev.values.clone());
        }
        Ok(ev.values.iter().map(|v| v / scale).collect())
    }
}
