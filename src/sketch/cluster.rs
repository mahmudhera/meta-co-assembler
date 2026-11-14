//! Turns a sparse similarity graph into capped, deterministic groups.
//!
//! Samples are indices into inputs sorted by sample ID. Groups are returned
//! ordered by their first sample, members in sample order; every sample
//! appears in exactly one group.

use std::collections::BTreeMap;

use super::Edge;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Policy {
    /// Average-linkage merging along every link while the caps allow.
    Merge,
    /// Average-linkage merging along links of at least this similarity.
    Threshold(f64),
    /// Seed with the sample sharing the most hashes, add its neighbours.
    Seeded,
}

/// Upper bounds per group; `None` is unlimited.
#[derive(Clone, Copy, Debug, Default)]
pub struct Caps {
    pub samples: Option<usize>,
    pub bases: Option<u64>,
}

impl Caps {
    fn admit(&self, samples: usize, bases: u64) -> bool {
        self.samples.is_none_or(|cap| samples <= cap) && self.bases.is_none_or(|cap| bases <= cap)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Cluster {
    pub members: Vec<usize>,
    /// The seed of a seeded group.
    pub seed: Option<usize>,
}

pub fn cluster(
    bases: &[u64],
    shared: &[u64],
    edges: &[Edge],
    policy: Policy,
    caps: Caps,
) -> Vec<Cluster> {
    let mut clusters = match policy {
        Policy::Merge => average_linkage(bases, edges.iter(), caps),
        Policy::Threshold(minimum) => average_linkage(
            bases,
            edges.iter().filter(|e| e.similarity >= minimum),
            caps,
        ),
        Policy::Seeded => seeded(bases, shared, edges, caps),
    };
    for found in &mut clusters {
        found.members.sort_unstable();
    }
    clusters.sort_by_key(|found| found.members[0]);
    clusters
}

/// Repeatedly merges the pair of linked groups with the highest mean link
/// similarity (absent links count as 0) whose union fits the caps; ties go
/// to the pair whose first samples come first. A group is identified by its
/// first sample.
fn average_linkage<'a>(
    bases: &[u64],
    edges: impl Iterator<Item = &'a Edge>,
    caps: Caps,
) -> Vec<Cluster> {
    let mut members: Vec<Vec<usize>> = (0..bases.len()).map(|i| vec![i]).collect();
    let mut group_bases = bases.to_vec();
    let mut links: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    for edge in edges {
        *links.entry((edge.a, edge.b)).or_default() += edge.similarity;
    }
    loop {
        let mut best: Option<((usize, usize), f64)> = None;
        for (&(g, h), &sum) in &links {
            let size = members[g].len() + members[h].len();
            if !caps.admit(size, group_bases[g] + group_bases[h]) {
                continue;
            }
            let mean = sum / (members[g].len() * members[h].len()) as f64;
            if best.is_none_or(|(_, top)| mean > top) {
                best = Some(((g, h), mean));
            }
        }
        let Some(((g, h), _)) = best else {
            break;
        };
        // g < h, so g stays the first sample of the merged group.
        let moved = std::mem::take(&mut members[h]);
        members[g].extend(moved);
        group_bases[g] += group_bases[h];
        links.remove(&(g, h));
        let touching: Vec<(usize, usize)> = links
            .keys()
            .filter(|&&(a, b)| a == h || b == h)
            .copied()
            .collect();
        for key in touching {
            let sum = links.remove(&key).expect("key was just listed");
            let other = if key.0 == h { key.1 } else { key.0 };
            *links.entry((g.min(other), g.max(other))).or_default() += sum;
        }
    }
    members
        .into_iter()
        .filter(|m| !m.is_empty())
        .map(|members| Cluster {
            members,
            seed: None,
        })
        .collect()
}

/// Seeds in order of shared hashes (most first, then sample order); each
/// seed takes its ungrouped neighbours by similarity while the caps allow.
fn seeded(bases: &[u64], shared: &[u64], edges: &[Edge], caps: Caps) -> Vec<Cluster> {
    let samples = bases.len();
    let mut neighbours: Vec<Vec<(f64, usize)>> = vec![Vec::new(); samples];
    for edge in edges {
        neighbours[edge.a].push((edge.similarity, edge.b));
        neighbours[edge.b].push((edge.similarity, edge.a));
    }
    for list in &mut neighbours {
        list.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));
    }
    let mut order: Vec<usize> = (0..samples).collect();
    order.sort_by(|&a, &b| shared[b].cmp(&shared[a]).then(a.cmp(&b)));
    let mut grouped = vec![false; samples];
    let mut clusters = Vec::new();
    for seed in order {
        if grouped[seed] {
            continue;
        }
        grouped[seed] = true;
        let mut members = vec![seed];
        let mut total = bases[seed];
        for &(_, candidate) in &neighbours[seed] {
            if grouped[candidate] {
                continue;
            }
            if caps.samples.is_some_and(|cap| members.len() >= cap) {
                break;
            }
            if caps.bases.is_some_and(|cap| total + bases[candidate] > cap) {
                continue;
            }
            grouped[candidate] = true;
            members.push(candidate);
            total += bases[candidate];
        }
        clusters.push(Cluster {
            members,
            seed: Some(seed),
        });
    }
    clusters
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edge(a: usize, b: usize, similarity: f64) -> Edge {
        Edge {
            a,
            b,
            shared: 1,
            similarity,
        }
    }

    /// Two tight triangles {0,2,4} and {1,3,5}, weakly linked by 4-5, and an
    /// unlinked sample 6.
    fn graph() -> Vec<Edge> {
        vec![
            edge(0, 2, 0.8),
            edge(0, 4, 0.7),
            edge(2, 4, 0.9),
            edge(1, 3, 0.6),
            edge(1, 5, 0.6),
            edge(3, 5, 0.5),
            edge(4, 5, 0.05),
        ]
    }

    fn groups(clusters: &[Cluster]) -> Vec<Vec<usize>> {
        clusters.iter().map(|c| c.members.clone()).collect()
    }

    fn caps(samples: Option<usize>, bases: Option<u64>) -> Caps {
        Caps { samples, bases }
    }

    #[test]
    fn merging_groups_similar_samples_and_pools_linked_samples_without_a_cap() {
        let bases = [10; 7];
        let shared = [5; 7];
        let capped = cluster(
            &bases,
            &shared,
            &graph(),
            Policy::Merge,
            caps(Some(3), None),
        );
        assert_eq!(groups(&capped), [vec![0, 2, 4], vec![1, 3, 5], vec![6]]);
        let open = cluster(&bases, &shared, &graph(), Policy::Merge, Caps::default());
        assert_eq!(
            groups(&open),
            [vec![0, 1, 2, 3, 4, 5], vec![6]],
            "an unlinked sample stays alone"
        );
        let by_bases = cluster(
            &bases,
            &shared,
            &graph(),
            Policy::Merge,
            caps(None, Some(20)),
        );
        assert_eq!(
            groups(&by_bases),
            [vec![0], vec![1, 3], vec![2, 4], vec![5], vec![6]],
            "the most similar pairs fill the base cap first"
        );
        assert!(capped.iter().all(|c| c.seed.is_none()));
    }

    #[test]
    fn a_threshold_keeps_weakly_linked_samples_apart_even_when_they_fit() {
        let bases = [10; 7];
        let shared = [5; 7];
        let clusters = cluster(
            &bases,
            &shared,
            &graph(),
            Policy::Threshold(0.5),
            Caps::default(),
        );
        assert_eq!(groups(&clusters), [vec![0, 2, 4], vec![1, 3, 5], vec![6]]);
        let strict = cluster(
            &bases,
            &shared,
            &graph(),
            Policy::Threshold(0.75),
            Caps::default(),
        );
        assert_eq!(
            groups(&strict),
            [vec![0, 2, 4], vec![1], vec![3], vec![5], vec![6]],
            "0 joins 2-4 through the 0.8 link"
        );
    }

    #[test]
    fn seeds_with_most_shared_hashes_take_their_nearest_neighbours() {
        let bases = [10; 7];
        let shared = [1, 9, 1, 1, 1, 1, 1];
        let clusters = cluster(
            &bases,
            &shared,
            &graph(),
            Policy::Seeded,
            caps(Some(2), None),
        );
        assert_eq!(
            clusters
                .iter()
                .map(|c| (c.seed.unwrap(), c.members.clone()))
                .collect::<Vec<_>>(),
            [
                (0, vec![0, 2]),
                (1, vec![1, 3]),
                (4, vec![4, 5]),
                (6, vec![6])
            ],
            "1 seeds first and takes 3 (tie with 5 goes to the earlier sample)"
        );
        let no_cap = cluster(&bases, &shared, &graph(), Policy::Seeded, Caps::default());
        assert_eq!(groups(&no_cap), [vec![0, 2, 4], vec![1, 3, 5], vec![6]]);
        let by_bases = cluster(
            &bases,
            &shared,
            &graph(),
            Policy::Seeded,
            caps(None, Some(20)),
        );
        assert!(by_bases.iter().all(|c| c.members.len() <= 2));
    }

    #[test]
    fn clustering_is_deterministic_and_covers_every_sample_once() {
        let bases = [10, 20, 30, 40, 50, 60, 70];
        let shared = [5; 7];
        for policy in [Policy::Merge, Policy::Threshold(0.1), Policy::Seeded] {
            let first = cluster(&bases, &shared, &graph(), policy, caps(Some(2), Some(100)));
            let mut reversed = graph();
            reversed.reverse();
            assert_eq!(
                first,
                cluster(&bases, &shared, &reversed, policy, caps(Some(2), Some(100)))
            );
            let mut all: Vec<usize> = first.iter().flat_map(|c| c.members.clone()).collect();
            all.sort_unstable();
            assert_eq!(all, (0..7).collect::<Vec<_>>());
            assert!(
                first.iter().all(|c| c.members.len() <= 2
                    && c.members.iter().map(|&m| bases[m]).sum::<u64>() <= 100)
            );
        }
    }
}
