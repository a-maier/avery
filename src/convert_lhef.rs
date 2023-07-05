use ahash::AHashMap;
use itertools::izip;
use lhef::{HEPRUP, HEPEUP, status::{INCOMING, OUTGOING, INTERMEDIATE_SPACELIKE, INTERMEDIATE_RESONANCE, INTERMEDIATE_DOC, INCOMING_BEAM}};
use particle_id::ParticleID;
use petgraph::{prelude::DiGraph, visit::NodeIndexable};

use crate::event::{Event, SampleInfo, Beam, CrossSection, WeightInfo, Scales, Particle, Status::{*, self}};

impl From<HEPRUP> for SampleInfo {
    fn from(source: HEPRUP) -> Self {
        let beam = izip!(source.IDBMUP, source.EBMUP)
            .map(|(id, e)| Beam {
                id: Some(ParticleID::new(id)),
                energy: Some(e),
            }).collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let cross_sections = source.XSECUP
            .into_iter()
            .zip(source.XERRUP.into_iter())
            .map(|(mean, err)| CrossSection{ mean, err: Some(err) })
            .collect();
        Self {
            generators: Default::default(), // TODO: update with LHEF version 3
            beam,
            pdf: source.PDFSUP.map(Some),
            cross_sections,
            process_ids: source.LPRUP,
            weight_type: Some(source.IDWTUP),
            info: source.info,
            attr: source.attr,
            ..Default::default()
        }
    }
}

impl From<(HEPRUP, HEPEUP)> for Event {
    fn from(source: (HEPRUP, HEPEUP)) -> Self {
        let (header, event) = source;
        // TODO: update with LHEF version 3
        let particles = izip!(
            event.IDUP,
            event.ISTUP,
            event.PUP,
            event.ICOLUP,
            event.SPINUP,
            event.VTIMUP,
        );
        let particles = particles.map(|(id, status, p, col, spin, lifetime)| {
            let m = p[4];
            let p = [p[3], p[0], p[1], p[2]];
            Particle {
                id: Some(ParticleID::new(id)),
                p: Some(p),
                m: Some(m),
                status: Some(from_i32(status)),
                spin: Some(spin),
                col: Some(col),
                lifetime: Some(lifetime),
                ..Default::default()
            }
        }).collect();
        let mut vertices: AHashMap<(i32, Option<i32>), Vec<_>> = AHashMap::new();
        for (n, parents) in event.MOTHUP.into_iter().enumerate() {
            if parents[0] == 0 {
                continue;
            }
            if parents[1] == 0 || parents[1] == parents[0] {
                vertices.entry((parents[0] - 1, None)).or_default().push(n);
            } else {
                vertices.entry((parents[0] - 1, Some(parents[1] - 1))).or_default().push(n);
            }
        }
        let vertices = Vec::from_iter(vertices);
        let mut end_vertices = AHashMap::new();
        for (end, ((in0, in1), _out)) in vertices.iter().enumerate() {
            let old_end = end_vertices.insert(*in0 as usize, end);
            assert!(old_end.is_none());
            if let Some(in1) = in1 {
                let old_end = end_vertices.insert(*in1 as usize, end);
                assert!(old_end.is_none());
            }
        }
        let mut g = DiGraph::new();
        for _ in 0..vertices.len() {
            g.add_node(Default::default());
        }
        // particles with parents (usually intermediate or outgoing)
        for (start, ((_, _), out)) in vertices.iter().enumerate() {
            for out in out {
                let end = if let Some(end) = end_vertices.remove(out) {
                    end
                } else {
                    g.add_node(Default::default());
                    g.node_count() - 1
                };
                let start = g.from_index(start);
                let end = g.from_index(end);
                g.add_edge(start, end, *out);
            }
        }
        // orphan particles (usually incoming)
        for (particle, end) in end_vertices {
            g.add_node(Default::default());
            let start = g.node_count() - 1;
            let start = g.from_index(start);
            let end = g.from_index(end);
            g.add_edge(start, end, particle);
        }

        Self {
            sample_info: header.into(),
            weights: vec![WeightInfo{
                weight: Some(event.XWGTUP),
                ..Default::default()
            }],
            scales: Scales{
                mu_r: Some(event.SCALUP),
                ..Default::default()
            },
            alpha_s: Some(event.AQCDUP),
            alpha: Some(event.AQEDUP),
            process_id: Some(event.IDRUP),
            particles,
            info: event.info,
            attr: event.attr,
            topology: g,
            ..Default::default()
        }
    }
}

fn from_i32(status: i32) -> Status {
    match status {
        INCOMING => Incoming,
        OUTGOING => Outgoing,
        INTERMEDIATE_SPACELIKE => IntermediateSpacelike,
        INTERMEDIATE_RESONANCE => IntermediateResonance,
        INTERMEDIATE_DOC => IntermediateDoc,
        INCOMING_BEAM => IncomingBeam,
        s => Unknown(s)
    }
}

impl From<HEPEUP> for Event {
    fn from(source: HEPEUP) -> Self {
        Self::from((HEPRUP{
            IDBMUP: Default::default(),
            EBMUP: Default::default(),
            PDFGUP: Default::default(),
            PDFSUP: Default::default(),
            IDWTUP: Default::default(),
            NPRUP: Default::default(),
            XSECUP: Default::default(),
            XERRUP: Default::default(),
            XMAXUP: Default::default(),
            LPRUP: Default::default(),
            info: Default::default(),
            attr: Default::default(),
        }, source))
    }
}

impl From<HEPRUP> for Event {
    fn from(source: HEPRUP) -> Self {
        Self::from((source, HEPEUP{
            NUP: Default::default(),
            IDRUP: Default::default(),
            XWGTUP: Default::default(),
            SCALUP: Default::default(),
            AQEDUP: Default::default(),
            AQCDUP: Default::default(),
            IDUP: Default::default(),
            ISTUP: Default::default(),
            MOTHUP: Default::default(),
            ICOLUP: Default::default(),
            PUP: Default::default(),
            VTIMUP: Default::default(),
            SPINUP: Default::default(),
            info: Default::default(),
            attr: Default::default(),
        }))
    }
}

impl From<Event> for (HEPRUP, HEPEUP) {
    fn from(mut source: Event) -> Self {
        let info = std::mem::take(&mut source.sample_info);
        (info.into(), source.into())
    }
}

impl From<Event> for HEPEUP {
    fn from(source: Event) -> Self {
        let mut parents = vec![[0, 0]; source.topology.edge_count()];
        let g = &source.topology;
        for vx in g.node_indices() {
            use petgraph::Direction::{Incoming, Outgoing};
            let incoming = Vec::from_iter(
                g.edges_directed(vx, Incoming).map(|e| *e.weight())
            );
            if incoming.is_empty() {
                continue;
            }
            assert!(incoming.len() <= 2);
            for out in g.edges_directed(vx, Outgoing) {
                let out = *out.weight();
                parents[out][0] = 1 + incoming[0]as i32;
                if let Some(parent) = incoming.last() {
                    parents[out][1] = 1 + *parent as i32;
                }
            }
        }
        Self {
            NUP: source.particles.len() as i32,
            IDRUP: source.process_id.unwrap_or_default(),
            XWGTUP: source.weights.first()
                .map(|w| w.weight.unwrap_or_default())
                .unwrap_or_default(),
            SCALUP: source.scales.mu_r.unwrap_or_default(),
            AQEDUP: source.alpha_s.unwrap_or_default(),
            AQCDUP: source.alpha.unwrap_or_default(),
            IDUP: source.particles.iter().map(|p| p.id.unwrap_or(ParticleID::new(0)).id()).collect(),
            ISTUP: source.particles.iter().map(|p| to_i32(p.status.unwrap_or(Unknown(0)))).collect(),
            MOTHUP: parents,
            ICOLUP: source.particles.iter().map(|p| p.col.unwrap_or_default()).collect(),
            PUP: source.particles.iter().map(|p| {
                let m = p.m.unwrap_or_default();
                let p = p.p.unwrap_or_default();
                [p[1], p[2], p[3], p[0], m]
            }).collect(),
            VTIMUP: source.particles.iter().map(|p| p.lifetime.unwrap_or_default()).collect(),
            SPINUP: source.particles.iter().map(|p| p.spin.unwrap_or_default()).collect(),
            info: source.info,
            attr: source.attr,
        }
    }
}

fn to_i32(status: Status) -> i32 {
    match status {
        Incoming => INCOMING,
        Outgoing => OUTGOING,
        IntermediateSpacelike => INTERMEDIATE_SPACELIKE,
        IntermediateResonance => INTERMEDIATE_RESONANCE,
        IntermediateDoc => INTERMEDIATE_DOC,
        IncomingBeam => INCOMING_BEAM,
        Unknown(s) => s,
    }
}

impl From<SampleInfo> for HEPRUP {
    fn from(source: SampleInfo) -> Self {
        let nsub = std::cmp::max(source.process_ids.len(), 1);
        Self {
            IDBMUP: source.beam.map(|b| b.id.map(|id| id.id()).unwrap_or_default()),
            EBMUP: source.beam.map(|b| b.energy.unwrap_or_default()),
            PDFGUP: Default::default(),
            PDFSUP: source.pdf.map(|p| p.unwrap_or_default()),
            IDWTUP: source.weight_type.unwrap_or_default(),
            NPRUP: nsub as i32,
            XSECUP: source.cross_sections.iter().map(|xs| xs.mean).collect(),
            XERRUP: source.cross_sections.iter().map(|xs| xs.err.unwrap_or_default()).collect(),
            XMAXUP: vec![0.; nsub], // TODO
            LPRUP: source.process_ids,
            info: source.info,
            attr: source.attr,
        }
    }
}
