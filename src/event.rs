use std::collections::{HashMap, BTreeMap};

use particle_id::ParticleID;
use petgraph::prelude::DiGraph;

/// Scattering event
#[derive(Clone, Debug, Default)]
pub struct Event {
    /// Event id
    pub id: Option<i32>,
    /// Global information about the current event sample
    pub sample_info: SampleInfo,
    /// Event weights
    pub weights: Vec<WeightInfo>,
    /// Scale settings
    pub scales: Scales,
    /// Value of the QCD coupling α_s
    pub alpha_s: Option<f64>,
    /// Value of the QED coupling α
    pub alpha: Option<f64>,
    /// ID of the process this event belongs to
    pub process_id: Option<i32>,
    /// Particles involved in the scattering
    pub particles: Vec<Particle>,
    /// Optional event information
    pub info: String,
    /// Optional additional structured information
    pub attr: HashMap<String, String>,
    /// Multiparton interaction
    pub mpi: Option<i32>,
    pub random_states: Vec<i32>,
    #[cfg(feature = "hepmc2")]
    pub heavy_ion_info: Option<hepmc2::event::HeavyIonInfo>,
    /// Event topology
    ///
    /// Edge weights correspond to the index in the `particles` vector.
    ///
    /// TODO: HepMC vertex information
    pub topology: DiGraph<Vertex, usize>
}

/// Global information about an event sample
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SampleInfo {
    /// Generators used to produce the sample
    pub generators: Vec<String>,
    /// Number of events
    pub nevents: Option<usize>,
    /// Collider beam information
    pub beam: [Beam; 2],
    /// Parton Distribution Functions used
    pub pdf: [Option<i32>; 2],
    /// Cross sections for the various subprocesses
    pub cross_sections: Vec<CrossSection>,
    /// Subprocess IDs
    pub process_ids: Vec<i32>,
    /// How to interprete the event weights
    pub weight_type: Option<i32>,
    /// Optional run information
    pub info: String,
    /// Optional additional structured information
    pub attr: HashMap<String, String>,
}

/// Event weight information
#[derive(Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct WeightInfo {
    /// The weight itself
    pub weight: Option<f64>,
    /// Weight name
    pub name: Option<String>,
    /// Factor multiplying the central renormalisation scale
    pub mu_r_factor: Option<f64>,
    /// Factor multiplying the central factorisation scale
    pub mu_f_factor: Option<f64>,
    /// PDF id
    pub pdf: Option<i32>,
    /// PDF id for second beam, if different from first beam
    pub pdf2: Option<i32>,

}

// Scales associated with event
#[derive(Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct Scales {
    /// The renormalisation scale
    pub mu_r: Option<f64>,
    /// The factorisation scale
    pub mu_f: Option<f64>,
    /// The suggested parton shower starting scale
    pub mu_ps: Option<f64>,
}

/// A particle
#[derive(Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct Particle {
    /// Particle type
    pub id:  Option<ParticleID>,
    /// Four-momentum
    pub p: Option<[f64; 4]>,
    /// Mass
    pub m: Option<f64>,
    /// Status
    pub status: Option<Status>,
    /// Spin angle
    pub spin: Option<f64>,
    /// Colour flow (LHEF)
    pub col: Option<[i32; 2]>,
    /// Colour flow (HepMC)
    pub flows: BTreeMap<i32, i32>,
    /// Lifetime
    pub lifetime: Option<f64>,
    /// θ polarisation angle
    pub theta: Option<f64>,
    /// φ polarisation angle
    pub phi: Option<f64>,
}

/// Beam information
#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct Beam {
    /// Particle type
    pub id: Option<ParticleID>,
    /// Energy in GeV
    pub energy: Option<f64>,
}

#[derive(Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct CrossSection {
    /// Mean value for the cross section
    pub mean: f64,
    /// Cross section errors
    pub err: Option<f64>,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum Status {
    /// Incoming particle
    Incoming,
    /// Outgoing particle
    Outgoing,
    /// Intermediate space-like propagator defining an x and Q^2 which should be preserved
    IntermediateSpacelike,
    /// Intermediate resonance, mass should be preserved
    IntermediateResonance,
    /// Intermediate resonance, for documentation only
    IntermediateDoc,
    /// Incoming beam particles at time t = −∞
    IncomingBeam,
    /// Unknown
    Unknown(i32),
}

#[derive(Clone, Debug, Default, PartialEq, PartialOrd)]
pub struct Vertex {
    pub status: Option<i32>,
    pub x: Option<f64>,
    pub y: Option<f64>,
    pub z: Option<f64>,
    pub t: Option<f64>,
    pub weights: Vec<f64>,
}
