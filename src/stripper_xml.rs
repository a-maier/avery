use stripper_xml::{Id, Momentum};

use crate::{
    event::{Particle, Reweight, Scales, Status, WeightInfo},
    Event,
};

impl From<stripper_xml::SubEvent> for Event {
    fn from(ev: stripper_xml::SubEvent) -> Self {
        let weight = WeightInfo {
            weight: Some(ev.weight),
            ..Default::default()
        };
        Self {
            weights: vec![weight],
            scales: Scales {
                mu_r: Some(ev.mu_r),
                mu_f: Some(ev.mu_f),
                mu_ps: None,
            },
            particles: ev.particles.into_iter().map(|p| p.into()).collect(),
            reweights: ev.reweight.into_iter().map(|rw| rw.into()).collect(),
            ..Default::default()
        }
    }
}

impl From<Event> for stripper_xml::SubEvent {
    fn from(ev: Event) -> Self {
        let particles = ev
            .particles
            .into_iter()
            .filter_map(|p| match p.status {
                Some(Status::Incoming | Status::Outgoing) => Some(p.into()),
                _ => None,
            })
            .collect();
        Self {
            weight: ev
                .weights
                .first()
                .map(|w| w.weight.unwrap_or_default())
                .unwrap_or_default(),
            mu_r: ev.scales.mu_r.unwrap_or_default(),
            mu_f: ev.scales.mu_r.unwrap_or_default(),
            particles,
            reweight: ev.reweights.into_iter().map(|rw| rw.into()).collect(),
        }
    }
}

impl From<Particle> for stripper_xml::Particle {
    fn from(particle: Particle) -> Self {
        Self {
            id: Id {
                status: particle.status.unwrap().into(),
                pdg_id: particle.id.unwrap(),
            },
            momentum: Momentum(particle.p.unwrap_or_default()),
        }
    }
}

impl From<stripper_xml::Particle> for Particle {
    fn from(particle: stripper_xml::Particle) -> Self {
        Self {
            id: Some(particle.id.pdg_id),
            p: Some(particle.momentum.0),
            status: Some(particle.id.status.into()),
            ..Default::default()
        }
    }
}

impl From<Reweight> for stripper_xml::Reweight {
    fn from(rw: Reweight) -> Self {
        Self {
            channel: rw.channel,
            reweights: stripper_xml::Reweights {
                x1: rw.x1,
                x2: rw.x2,
                log_coeff: rw.log_coeff,
            },
        }
    }
}

impl From<stripper_xml::Reweight> for Reweight {
    fn from(rw: stripper_xml::Reweight) -> Self {
        Self {
            channel: rw.channel,
            x1: rw.reweights.x1,
            x2: rw.reweights.x2,
            log_coeff: rw.reweights.log_coeff,
        }
    }
}

impl From<Status> for stripper_xml::Status {
    fn from(status: Status) -> Self {
        match status {
            Status::Incoming => stripper_xml::Status::Incoming,
            Status::Outgoing => stripper_xml::Status::Outgoing,
            _ => panic!("Can only convert incoming and outgoing particles to STRIPPER-XML format")
        }
    }
}

impl From<stripper_xml::Status> for Status {
    fn from(status: stripper_xml::Status) -> Self {
        match status {
            stripper_xml::Status::Incoming => Status::Incoming,
            stripper_xml::Status::Outgoing => Status::Outgoing,
        }
    }
}
