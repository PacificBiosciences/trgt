use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TrgtPreset {
    WGS,
    Targeted,
}

impl FromStr for TrgtPreset {
    type Err = &'static str;
    fn from_str(preset: &str) -> Result<Self, Self::Err> {
        match preset {
            "wgs" => Ok(TrgtPreset::WGS),
            "targeted" => Ok(TrgtPreset::Targeted),
            _ => Err("Invalid preset. Options are: wgs, targeted"),
        }
    }
}
