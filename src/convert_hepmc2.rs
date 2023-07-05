use ahash::AHashMap;
use hepmc2::event::{FourVector, Vertex, PdfInfo};
use itertools::izip;
use particle_id::ParticleID;
use petgraph::{visit::NodeIndexable, prelude::DiGraph};

use crate::event::{Event, WeightInfo, Scales, SampleInfo, Particle, Status, CrossSection, Beam};

const HEPMC_OUTGOING: i32 = 1;
const HEPMC_DECAYED: i32 = 2;
const HEPMC_DOC: i32 = 3;
const HEPMC_INCOMING: i32 = 4;

impl From<hepmc2::Event> for Event {
    fn from(mut source: hepmc2::Event) -> Self {
        let efact = if source.energy_unit == hepmc2::event::EnergyUnit::MEV {
            1e-3
        } else {
            1.
        };
        let mut beam = [Beam::default(), Beam::default()];
        let root_vx = source.vertices.iter().position(
            |vx| vx.particles_in.is_empty() && vx.particles_out.len() <= 2
        );
        if let Some(root_vx) = root_vx {
            let root_vx = source.vertices.swap_remove(root_vx);
            for (beam, beam_p) in izip!(&mut beam, root_vx.particles_out.into_iter()) {
                beam.energy = Some(beam_p.p[0]);
                beam.id = Some(ParticleID::new(beam_p.id));
            }
        }
        let sample_info = SampleInfo {
            pdf: source.pdf_info.pdf_id.map(Some),
            cross_sections: vec![source.xs.into()],
            beam,
            ..Default::default()
        };
        let weights = izip!(
            source.weights,
            source.weight_names
        ).map(|(weight, name)| {
            WeightInfo {
                weight: Some(weight),
                name: Some(name),
                ..Default::default()
            }
        }).collect();
        let mut barcode_to_idx = AHashMap::new();
        for (idx, vx) in source.vertices.iter().enumerate() {
            let seen = barcode_to_idx.insert(vx.barcode, idx);
            assert!(seen.is_none());
        }
        let mut particles = Vec::new();
        let mut topology: DiGraph<crate::event::Vertex, _> = DiGraph::new();
        for _ in 0..source.vertices.len() {
            topology.add_node(Default::default());
        }
        for vx in source.vertices {
            let idx = *barcode_to_idx.get(&vx.barcode).unwrap();
            {
                let idx = topology.from_index(idx);
                let node = topology.node_weight_mut(idx).unwrap();
                node.status = Some(vx.status);
                node.x = Some(vx.x);
                node.y = Some(vx.y);
                node.z = Some(vx.z);
                node.t = Some(vx.t);
                node.weights = vx.weights;
            }
            // TODO: here we assume that everything is sane
            //       make it more robust?
            let mut vx_particles = Vec::new();
            // only consider incident particles with incoming status,
            // all other incident particles will be treated when
            // looking at the vertex where they originate
            let incoming = vx.particles_in.into_iter()
                .filter(|p| p.status == HEPMC_INCOMING);
            for incoming in incoming {
                topology.add_node(Default::default());
                let start = topology.from_index(topology.node_count() - 1);
                let end = topology.from_index(idx);
                topology.add_edge(start, end, particles.len() + vx_particles.len());
                vx_particles.push(incoming);
            }
            for out in vx.particles_out {
                let end_vtx = barcode_to_idx.get(&out.end_vtx).copied();
                let end_vtx = end_vtx.unwrap_or_else(|| {
                    topology.add_node(Default::default());
                    topology.node_count() - 1
                });
                let start = topology.from_index(idx);
                let end = topology.from_index(end_vtx);
                topology.add_edge(start, end, particles.len() + vx_particles.len());
                vx_particles.push(out);
            }

            for particle in vx_particles {
                let conv = Particle {
                    id: Some(ParticleID::new(particle.id)),
                    p: Some(particle.p.0.map(|p| efact * p)),
                    m: Some(efact * particle.m),
                    status: Some(from_i32(particle.status)),
                    flows: particle.flows,
                    ..Default::default()

                };
                particles.push(conv)
            }
        }

        Self {
            id: Some(source.number),
            sample_info,
            weights,
            scales: Scales {
                mu_r: Some(efact * source.scale),
                mu_f: Some(efact * source.pdf_info.scale),
                mu_ps: Default::default(),
            },
            alpha_s: Some(source.alpha_qcd),
            alpha: Some(source.alpha_qed),
            process_id: Some(source.signal_process_id),
            particles,
            topology,
            mpi: Some(source.mpi),
            random_states: source.random_states,
            heavy_ion_info: source.heavy_ion_info,
            ..Default::default()
        }
    }
}

impl From<Event> for hepmc2::Event {
    fn from(source: Event) -> Self {
        let mut vertices = Vec::new();
        let mut beam_particles = Vec::with_capacity(2);
        for beam in source.sample_info.beam {
            let e = beam.energy.unwrap_or_default();
            let p = hepmc2::event::Particle {
                id: beam.id.map(|id| id.id()).unwrap_or_default(),
                p: FourVector([e, 0., 0., 0.]),
                status: HEPMC_INCOMING,
                ..Default::default()
            };
            beam_particles.push(p);
        }
        let root_vx = Vertex {
            barcode: 0,
            particles_out: beam_particles,
            ..Default::default()
        };
        vertices.push(root_vx);

        let nparticles = source.particles.len();
        let mut particles = Vec::with_capacity(nparticles);
        for particle in source.particles {
            let mut p = hepmc2::event::Particle {
                id: particle.id.map(|p| p.id()).unwrap_or_default(),
                p: hepmc2::event::FourVector(particle.p.unwrap_or_default()),
                m: particle.m.unwrap_or_default(),
                theta: particle.theta.unwrap_or_default(),
                phi: particle.phi.unwrap_or_default(),
                flows: particle.flows,
                ..Default::default()
            };
            p.status = particle.status.map(to_i32).unwrap_or_default();
            particles.push(p);
        }
        let g = &source.topology;
        // ensure there is always at least one vertex
        if g.node_count() == 0 {
            const END_VTX: i32 = -1;
            let (incoming, outgoing) = particles.into_iter()
                .partition(|p| p.status == HEPMC_INCOMING);
            let mut vx = Vertex {
                barcode: END_VTX,
                particles_in: incoming,
                particles_out: outgoing,
                ..Default::default()
            };
            for incoming in &mut vx.particles_in {
                incoming.end_vtx = END_VTX;
            }
            vertices.push(vx);
        } else {
            for (n, vx) in g.node_indices().enumerate() {
                use petgraph::Direction::{Incoming, Outgoing};
                // TODO: is this barcode correct?
                let barcode = - (n as i32) - 1;
                let incoming = g.edges_directed(vx, Incoming);
                let outgoing = g.edges_directed(vx, Outgoing);
                // HepMC does not like vertices with no incoming
                // particles and does not include the end vertex for
                // final-state particles
                if incoming.clone().count() == 0
                    || outgoing.count() == 0 {
                        continue;
                    }
                for n in incoming.map(|e| *e.weight()) {
                    particles[n].end_vtx = barcode;
                }
            }

            for (n, vx) in g.node_indices().enumerate() {
                let barcode = - (n as i32) - 1;
                use petgraph::Direction::{Incoming, Outgoing};

                let incoming = Vec::from_iter(
                    g.edges_directed(vx, Incoming)
                        .map(|e| particles[*e.weight()].clone())
                );
                if incoming.is_empty() {
                    continue;
                }
                let outgoing = Vec::from_iter(
                    g.edges_directed(vx, Outgoing)
                        .map(|e| particles[*e.weight()].clone())
                );
                if outgoing.is_empty() {
                    continue;
                }
                assert!(incoming.iter().all(|p| p.end_vtx == barcode));
                let vx = g.node_weight(vx).unwrap();

                let vx = Vertex {
                    barcode,
                    particles_in: incoming,
                    particles_out: outgoing,
                    status: vx.status.unwrap_or_default(),
                    x: vx.x.unwrap_or_default(),
                    y: vx.y.unwrap_or_default(),
                    z: vx.z.unwrap_or_default(),
                    t: vx.t.unwrap_or_default(),
                    weights: vx.weights.clone(),

                };
                vertices.push(vx);
            }
        }
        let mut xs = hepmc2::event::CrossSection::default();
        for channel_xs in source.sample_info.cross_sections {
            xs.cross_section += channel_xs.mean;
            if let Some(err) = channel_xs.err {
                xs.cross_section_error += err * err;
            }
        }
        xs.cross_section_error = xs.cross_section_error.sqrt();
        let pdf_info = PdfInfo {
            parton_id: source.sample_info.beam
                .map(|b| b.id.map(|id| id.id()).unwrap_or_default()),
            x: Default::default(), // TODO: compute
            scale: source.scales.mu_f.unwrap_or_default(),
            xf: Default::default(), // TODO
            pdf_id: source.sample_info.pdf.map(|p| p.unwrap_or_default()),
        };
        Self {
            number: source.id.unwrap_or_default(),
            mpi: source.mpi.unwrap_or_default(),
            scale: source.scales.mu_r.unwrap_or_default(),
            alpha_qcd: source.alpha_s.unwrap_or_default(),
            alpha_qed: source.alpha.unwrap_or_default(),
            signal_process_id: source.process_id.unwrap_or_default(),
            signal_process_vertex: Default::default(),
            random_states: source.random_states,
            weights: source.weights.iter().map(|w| w.weight.unwrap_or_default()).collect(),
            weight_names: source.weights.into_iter().map(|w| w.name.unwrap_or_default()).collect(),
            vertices,
            xs,
            pdf_info,
            energy_unit: hepmc2::event::EnergyUnit::GEV,
            length_unit: hepmc2::event::LengthUnit::MM,
            heavy_ion_info: source.heavy_ion_info,
        }
    }
}

fn from_i32(status: i32) -> Status {
    use Status::*;
    match status {
        HEPMC_INCOMING => Incoming,
        HEPMC_OUTGOING => Outgoing,
        HEPMC_DECAYED => IntermediateResonance,
        HEPMC_DOC => IntermediateDoc,
        s => Unknown(s),
    }
}

fn to_i32(status: Status) -> i32 {
    use Status::*;
    match status {
        Incoming | IncomingBeam => HEPMC_INCOMING,
        IntermediateResonance
            | IntermediateSpacelike => HEPMC_DECAYED,
        IntermediateDoc => HEPMC_DOC,
        Outgoing =>  HEPMC_OUTGOING,
        Unknown(s) => s,
    }
}

impl From<hepmc2::event::CrossSection> for CrossSection {
    fn from(source: hepmc2::event::CrossSection) -> Self {
        Self {
            mean: source.cross_section,
            err: Some(source.cross_section_error),
        }
    }
}

impl From<CrossSection> for hepmc2::event::CrossSection {
    fn from(source: CrossSection) -> Self {
        Self {
            cross_section: source.mean,
            cross_section_error: source.err.unwrap_or_default(),
        }
    }
}
