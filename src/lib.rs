pub mod event;
#[cfg(feature = "lhef")]
pub mod convert_lhef;
#[cfg(feature = "hepmc2")]
pub mod convert_hepmc2;

use crate::event::Event;

pub fn convert<T, U>(event: T) -> U
where
    T: Into<Event>,
    U: From<Event>,
{
    U::from(event.into())
}
