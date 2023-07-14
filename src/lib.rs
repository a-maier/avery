pub mod event;
#[cfg(feature = "lhef")]
pub mod lhef;
#[cfg(feature = "hepmc2")]
pub mod hepmc2;

pub use crate::event::Event;

pub fn convert<T, U>(event: T) -> U
where
    T: Into<Event>,
    U: From<Event>,
{
    U::from(event.into())
}
