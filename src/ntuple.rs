// TODO: support for x1p, x2p, alphas_power
use itertools::izip;
use ntuple::event::Part;
use particle_id::ParticleID;

use crate::{
    event::{
        Beam, Particle, SampleInfo, Scales,
        Status::{self, Outgoing},
        WeightInfo,
    },
    util::{extract_inc_info, IncomingInfo},
    Event,
};

impl From<ntuple::Event> for Event {
    fn from(source: ntuple::Event) -> Self {
        let ntuple::Event {
            id,
            nparticle,
            px,
            py,
            pz,
            energy,
            alphas,
            pdg_code,
            weight,
            weight2,
            me_weight,
            me_weight2,
            x1,
            x2,
            id1,
            id2,
            fac_scale,
            ren_scale,
            user_weights,
            part,
            ..
        } = source;
        let nout = nparticle as usize;
        let mut particles = Vec::with_capacity(nout + 2);
        let e_tot = energy.iter().sum::<f32>() as f64;
        let pz_tot = pz.iter().sum::<f32>() as f64;
        let e_in = [(e_tot - pz_tot) / 2., (e_tot + pz_tot) / 2.];
        particles.push(Particle {
            id: Some(ParticleID::new(id1)),
            p: Some([0., 0., -e_in[0], e_in[0]]),
            status: Some(Status::Incoming),
            ..Default::default()
        });
        particles.push(Particle {
            id: Some(ParticleID::new(id2)),
            p: Some([0., 0., e_in[1], e_in[1]]),
            status: Some(Status::Incoming),
            ..Default::default()
        });
        for (pdg_code, e, px, py, pz) in izip!(pdg_code, energy, px, py, pz) {
            particles.push(Particle {
                id: Some(ParticleID::new(pdg_code)),
                p: Some([e as f64, px as f64, py as f64, pz as f64]),
                status: Some(Status::Outgoing),
                ..Default::default()
            });
        }

        let beam = [
            Beam {
                id: None,
                energy: Some(e_in[0] / x1),
            },
            Beam {
                id: None,
                energy: Some(e_in[1] / x2),
            },
        ];

        let sample_info = SampleInfo {
            beam,
            process_ids: vec![0, 1, 2, 3],
            ..Default::default()
        };

        let scales = Scales {
            mu_r: Some(ren_scale),
            mu_f: Some(fac_scale),
            mu_ps: None,
        };

        let mut weights = Vec::with_capacity(user_weights.len() + 4);
        weights.push(WeightInfo {
            weight: Some(weight),
            ..Default::default()
        });
        weights.push(WeightInfo {
            weight: Some(weight2),
            name: Some("2".to_owned()),
            ..Default::default()
        });
        weights.push(WeightInfo {
            weight: Some(me_weight),
            name: Some("ME".to_owned()),
            ..Default::default()
        });
        weights.push(WeightInfo {
            weight: Some(me_weight2),
            name: Some("ME2".to_owned()),
            ..Default::default()
        });
        weights.extend(user_weights.into_iter().map(|w| WeightInfo {
            weight: Some(w),
            ..Default::default()
        }));

        Self {
            id: Some(id),
            sample_info,
            weights,
            scales,
            alpha_s: Some(alphas),
            process_id: Some(part_to_id(part)),
            particles,
            ..Default::default()
        }
    }
}

impl From<Event> for ntuple::Event {
    fn from(source: Event) -> Self {
        let IncomingInfo { parton_id, x } = extract_inc_info(&source);

        let outgoing = source
            .particles
            .iter()
            .filter(|p| p.status == Some(Outgoing));
        let nparticle = outgoing.clone().count();
        let mut pdg_code = Vec::with_capacity(nparticle);
        let mut energy = Vec::with_capacity(nparticle);
        let mut px = Vec::with_capacity(nparticle);
        let mut py = Vec::with_capacity(nparticle);
        let mut pz = Vec::with_capacity(nparticle);

        for out in outgoing {
            pdg_code.push(out.id.map(|id| id.id()).unwrap_or_default());
            let p = out.p.unwrap_or_default();
            energy.push(p[0] as f32);
            px.push(p[1] as f32);
            py.push(p[2] as f32);
            pz.push(p[3] as f32);
        }

        let mut weights = source.weights;
        let mut me_weight2 = 0.;
        let mut me_weight = 0.;
        let mut weight2 = 0.;

        weights.retain(|w| match w.name.as_deref() {
            Some("ME") => {
                me_weight = w.weight.unwrap_or_default();
                false
            }
            Some("ME2") => {
                me_weight2 = w.weight.unwrap_or_default();
                false
            }
            Some("2") => {
                weight2 = w.weight.unwrap_or_default();
                false
            }
            _ => true,
        });

        let mut weights = weights.into_iter();
        let weight = weights
            .next()
            .map(|w| w.weight.unwrap_or_default())
            .unwrap_or_default();
        let user_weights = weights.filter_map(|w| w.weight).collect();

        Self {
            id: source.id.unwrap_or_default(),
            nparticle: nparticle as i32,
            px,
            py,
            pz,
            energy,
            alphas: source.alpha_s.unwrap_or_default(),
            pdg_code,
            weight,
            weight2,
            me_weight,
            me_weight2,
            x1: x[0],
            x2: x[1],
            x1p: 0., // TODO
            x2p: 0., // TODO
            id1: parton_id[0],
            id2: parton_id[1],
            fac_scale: source.scales.mu_f.unwrap_or_default(),
            ren_scale: source.scales.mu_r.unwrap_or_default(),
            user_weights,
            part: id_to_part(source.process_id),
            alphas_power: 0, // TODO
        }
    }
}

fn id_to_part(process_id: Option<i32>) -> Part {
    match process_id.unwrap_or_default() {
        1 => Part::I,
        2 => Part::R,
        3 => Part::V,
        _ => Part::B,
    }
}

fn part_to_id(part: Part) -> i32 {
    match part {
        Part::B => 0,
        Part::I => 1,
        Part::R => 2,
        Part::V => 3,
    }
}
