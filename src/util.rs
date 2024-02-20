use crate::{event::Status, Event};

pub(crate) struct IncomingInfo {
    pub parton_id: [i32; 2],
    pub x: [f64; 2],
}

pub(crate) fn extract_inc_info(ev: &Event) -> IncomingInfo {
    let incoming = ev
        .particles
        .iter()
        .filter(|p| p.status == Some(Status::Incoming));
    let mut parton_id = [0, 0];
    let mut x = [0., 0.];
    for particle in incoming {
        let Some(p) = particle.p else { continue };
        let idx = if p[3] < 0. { 0 } else { 1 };
        parton_id[idx] = particle.id.map(|id| id.id()).unwrap_or_default();
        let beam = &ev.sample_info.beam[idx];
        if let Some(e) = beam.energy {
            x[idx] = p[0] / e;
        }
    }
    IncomingInfo { parton_id, x }
}
