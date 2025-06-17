use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::hash::Hash;
use std::io;
use std::io::ErrorKind::InvalidData;
use std::io::prelude::*;
use std::num::ParseIntError;
use std::str::FromStr;

fn newkey<T: FromStr>(id: &str) -> Result<T, T::Err> {
    id.parse::<T>()
}

/// 32 bits are enough for most sequence lengths...
type Iv = Interval<u32, bool>;
///
/// T: model ID type (Rfam, Dfam, etc. use numerics, allowing us to save space)
type Matches<T> = HashMap<String, HashMap<T, Vec<Iv>>>;

/// How many hits a model has
type OwnHits<T> = HashMap<T, u64>;

/// Skim a FASTA file to extract sequence name/start-end
fn skim_fasta<F: BufRead, T: FromStr + Eq + Clone + Send + Sync + Hash>(
    m: &mut Matches<T>,
    h: &mut OwnHits<T>,
    mut f: F,
    model_name: &str,
    use_id_version: bool,
) -> io::Result<()> {
    let mut buffer = String::new();
    let mut count: u64 = 0;
    let model: T = newkey(model_name)
        .map_err(|_: <T as FromStr>::Err| io::Error::new(InvalidData, model_name))?;
    while f.read_line(&mut buffer)? > 0 {
        if buffer.starts_with('>') {
            // IDIDID.ver/start-end garbage garbage
            let endpos = buffer.find(' ').unwrap_or(buffer.len());
            buffer.truncate(endpos);
            let slash = buffer
                .find('/')
                .ok_or(io::Error::new(InvalidData, "No '/' found in FASTA ID line"))?;
            let dash = buffer
                .rfind('-')
                .ok_or(io::Error::new(InvalidData, "No '-' found in FASTA ID line"))?;
            let dot = buffer.find('.').unwrap_or(slash);
            let id = if use_id_version {
                &buffer[1..dot] // Skip '>'
            } else {
                &buffer[1..slash] // Skip '>'
            };
            let mut start: u32 = buffer[slash + 1..dash]
                .parse::<u32>()
                .map_err(|e| io::Error::new(InvalidData, e.to_string()))?;
            let mut stop: u32 = buffer[dash + 1..]
                .parse::<u32>()
                .map_err(|e| io::Error::new(InvalidData, e.to_string()))?;
            let rev = start > stop;
            if rev {
                std::mem::swap(&mut start, &mut stop);
            }
            let entry = m
                .entry(id.to_string())
                .or_insert_with(HashMap::new)
                .entry(model.clone())
                .or_insert_with(Vec::new);
            entry.push(Iv {
                start,
                stop,
                val: rev,
            });
            count += 1;
        }
        buffer.clear();
    }
    h.entry(model)
        .and_modify(|e| *e += count as u64)
        .or_insert(count as u64);
    Ok(())
}

/*
#seqname	modelacc	modelname	bitsc	eval	hmmstart	hmmend	hmmlen	strand	alistart    aliend	    envstart	envend	    seqlen
chr10	    DF000000001	MIR	        100.0	2.1e-25	6	        255	    262	    -	    112667059	112666809	112667062	112666804	133797422
*/
/// Skim a Dfam TSV file to extract sequence name/start-end
fn skim_dfam<F: BufRead, T: FromStr + Eq + Clone + Send + Sync + Hash>(
    m: &mut Matches<T>,
    h: &mut OwnHits<T>,
    mut f: F,
    seqname_prefix: &str,
    trim_acc: bool,
) -> io::Result<()> {
    let mut buffer = String::new();
    while f.read_line(&mut buffer)? > 0 {
        if !buffer.starts_with('#') && !buffer.trim().is_empty() {
            let parts: Vec<&str> = buffer.split('\t').collect();
            if parts.len() < 14 {
                continue; // Not enough columns, skip this line
            }
            let id = seqname_prefix.to_owned() + parts[0].trim();
            let mut model = parts[1].trim();
            if trim_acc {
                if let Some(pos) = model.find(|c: char| '0' <= c && c <= '9') {
                    model = &model[pos..];
                }
            }
            let key: T = newkey(model)
                .map_err(|e: <T as FromStr>::Err| io::Error::new(InvalidData, model))?;
            let rev = parts[7].trim() == "-";
            let start: u32 = parts[9]
                .trim()
                .parse()
                .map_err(|e: ParseIntError| io::Error::new(InvalidData, e.to_string()))?;
            let stop: u32 = parts[10]
                .trim()
                .parse()
                .map_err(|e: ParseIntError| io::Error::new(InvalidData, e.to_string()))?;
            let entry = m
                .entry(id)
                .or_insert_with(HashMap::new)
                .entry(key.clone())
                .or_insert_with(Vec::new);
            entry.push(Iv {
                start,
                stop,
                val: rev,
            });
            h.entry(key).and_modify(|e| *e += 1).or_insert(1);
        }
        buffer.clear();
    }
    Ok(())
}

/// Interpro TSV
// (!)Accession	Source Database	Name	Tax ID	Tax Name	Length	Entry-Accession	(!)Matches
// A0A010Q6J1	unreviewed	Piwi domain-containing protein	1445577	Colletotrichum fioriniae PJ7	973	PF02171	603..907
fn skim_interpro<F: BufRead, T: FromStr + Eq + Clone + Send + Sync + Hash>(
    m: &mut Matches<T>,
    h: &mut OwnHits<T>,
    mut f: F,
    model_name: &str,
) -> io::Result<()> {
    let mut buffer = String::new();
    let mut count: u64 = 0;
    let key: T = newkey(model_name)
        .map_err(|_: <T as FromStr>::Err| io::Error::new(InvalidData, model_name))?;
    while f.read_line(&mut buffer)? > 0 {
        if !buffer.starts_with('#') && !buffer.trim().is_empty() {
            let parts: Vec<&str> = buffer.split('\t').collect();
            if parts.len() < 9 {
                continue; // Not enough columns, skip this line
            }
            let id = parts[0].trim();
            let start_stop: Vec<&str> = parts[8].split("..").collect();
            if start_stop.len() != 2 {
                continue; // Invalid range format
            }
            let start: u32 = start_stop[0]
                .parse()
                .map_err(|e: ParseIntError| io::Error::new(InvalidData, e.to_string()))?;
            let stop: u32 = start_stop[1]
                .parse()
                .map_err(|e: ParseIntError| io::Error::new(InvalidData, e.to_string()))?;
            let rev = false; // InterPro entries are not strand-specific
            let entry = m
                .entry(id.to_string())
                .or_insert_with(HashMap::new)
                .entry(key.clone())
                .or_insert_with(Vec::new);
            entry.push(Iv {
                start,
                stop,
                val: rev,
            });
            count += 1;
        }
        buffer.clear();
    }
    h.entry(key)
        .and_modify(|e| *e += count as u64)
        .or_insert(count as u64);
    Ok(())
}

fn skim_filename<T: FromStr + Eq + Clone + Send + Sync + Hash>(
    m: &mut Matches<T>,
    h: &mut OwnHits<T>,
    filename: &str,
    use_id_version: bool,
    trim_acc: bool,
) -> io::Result<()> {
    let file = std::fs::File::open(filename)?;
    let reader = io::BufReader::new(file);
    let (mut basename, ext) = filename.rsplit_once('.').unwrap_or((filename, ""));
    // check file extension
    if basename.starts_with("RF") {
        if trim_acc {
            if let Some(pos) = basename.find(|c: char| '0' <= c && c <= '9') {
                basename = &basename[pos..];
            }
        }
        skim_fasta(m, h, reader, basename, use_id_version)
    } else if basename.starts_with("DF") {
        let (genome_name, _) = basename.rsplit_once('.').unwrap_or((basename, ""));
        let (_, genome_name) = genome_name.rsplit_once('.').unwrap_or(("", genome_name));
        skim_dfam(m, h, reader, &(genome_name.to_owned() + "."), trim_acc)
    } else if basename.starts_with("proteins-matching") {
        if trim_acc {
            if let Some(pos) = basename.find(|c: char| '0' <= c && c <= '9') {
                basename = &basename[pos..];
            }
        }
        skim_interpro(m, h, reader, basename)
    } else {
        Err(io::Error::new(
            InvalidData,
            "I only directly understand filenames starting with DF or RF",
        ))
    }
}

/// Prepare intervals by sorting and merging overlaps
fn prep_intervals<T: FromStr + Eq + Clone + Send + Sync + Hash>(matches: &mut Matches<T>) {
    for (_, ivs) in matches.iter_mut() {
        for (_, intervals) in ivs.iter_mut() {
            let mut lapper = Lapper::new(intervals.clone());
            lapper.merge_overlaps();
            *intervals = lapper.intervals;
        }
    }
}

/// Count of observed overlapping regions (M^M)
type Commons<T> = HashMap<(T, T), u64>;
/// Obs and total count of any matched regions (M^M + 1)
type ComTot<T> = (Commons<T>, u64);

fn pairkey<T: Clone + PartialOrd>(a: &T, b: &T) -> (T, T) {
    if a < b {
        (a.clone(), b.clone())
    } else {
        (b.clone(), a.clone())
    }
}

/// Count overlaps between models
fn count_overlaps<T: FromStr + Eq + Clone + Send + Sync + Hash + PartialOrd>(
    matches: &Matches<T>,
) -> ComTot<T> {
    let mut commons: Commons<T> = HashMap::new();
    let mut total: u64 = 0;
    for (_, ms) in matches.iter() {
        for (m1, ivs1) in ms.iter() {
            let lapper = Lapper::new(ivs1.clone());
            total += ivs1.len() as u64;
            for (m2, ivs2) in ms.iter() {
                let overlaps = if m1 == m2 {
                    ivs1.len() as u64
                } else {
                    ivs2.into_iter()
                        .map(|iv| lapper.count(iv.start, iv.stop) as u64)
                        .sum()
                };
                *commons.entry(pairkey(m1, m2)).or_insert(0) += overlaps as u64;
            }
        }
    }
    (commons, total)
}

/// Count overlaps without considering intervals
fn count_overlaps_no_interval<T: FromStr + Eq + Clone + Send + Sync + Hash + PartialOrd>(
    matches: &Matches<T>,
) -> ComTot<T> {
    let mut commons: Commons<T> = HashMap::new();
    let mut total: u64 = 0;
    for (_, ms) in matches.iter() {
        for (m1, ivs1) in ms.iter() {
            total += ivs1.len() as u64;
            for (m2, _) in ms.iter() {
                *commons.entry(pairkey(m1, m2)).or_insert(0) += 1;
            }
        }
    }
    (commons, total)
}

/// Count of observed overlapping regions (M^M)
type Scores2D<T> = HashMap<(T, T), f64>;

/// Sraw(A,B) = Obs / (1 + (NA * NB) / Tot)
fn raw_score<T: FromStr + Eq + Clone + Send + Sync + Hash + PartialOrd>(
    commons: &ComTot<T>,
    hits: &OwnHits<T>,
) -> Scores2D<T> {
    let mut scores = HashMap::new();
    for ((m1, m2), count) in commons.0.iter() {
        scores.insert(
            (m1.clone(), m2.clone()),
            *count as f64 / (1.0 + (hits[m1] as f64 * hits[m2] as f64) / commons.1 as f64),
        );
    }

    scores
}

/// S(A) = Sum(Sraw(A, x) forall x in M) / (100 + Tot)
type Scores1D<T> = HashMap<T, f64>;

fn model_score<T: FromStr + Eq + Clone + Send + Sync + Hash>(
    scores: &Scores2D<T>,
    tot: u64,
) -> Scores1D<T> {
    let mut model_scores: Scores1D<T> = HashMap::new();
    for ((m1, m2), score) in scores.iter() {
        *model_scores.entry(m1.clone()).or_insert(0.0) += score;
        if m1 == m2 {
            continue;
        }
        *model_scores.entry(m2.clone()).or_insert(0.0) += score;
    }
    for score in model_scores.values_mut() {
        *score /= 100.0 + tot as f64;
    }
    model_scores
}

/// Snorm(A,B) = Sraw(A,B) / max(S(A), S(B))
fn norm_score<T: FromStr + Eq + Clone + Send + Sync + Hash>(
    raw_scores: &Scores2D<T>,
    model_scores: &Scores1D<T>,
) -> Scores2D<T> {
    let mut norm_scores = HashMap::new();
    for ((m1, m2), score) in raw_scores.iter() {
        let max_score = model_scores[m1].max(model_scores[m2]);
        if max_score > 0.0 {
            norm_scores.insert((m1.clone(), m2.clone()), score / max_score);
        } else {
            norm_scores.insert((m1.clone(), m2.clone()), 0.0);
        }
    }
    norm_scores
}

/// Raw and Normalized scores for model pairs
type Scores2DN<T> = HashMap<(T, T), (f64, f64)>;

fn scores<T: FromStr + Eq + Clone + Send + Sync + Hash + PartialOrd>(
    commons: &ComTot<T>,
    hits: &OwnHits<T>,
) -> Scores2DN<T> {
    let raw = raw_score(commons, hits);
    let model_scores = model_score(&raw, commons.1);
    let norm = norm_score(&raw, &model_scores);

    // Combine raw and normalized scores
    let mut combined_scores = HashMap::new();
    for ((m1, m2), score) in norm.iter() {
        combined_scores.insert(
            (m1.clone(), m2.clone()),
            (raw[&(m1.clone(), m2.clone())], *score),
        );
    }

    combined_scores
}
