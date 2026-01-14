# xsteamrs

An IAPWS-IF97 water/steam thermophysical properties library for Rust. It provides two usage styles:

- **IF97 low-level functions**: region-based `*_p_t` / `*_rho_t` functions using IF97 conventional SI units (`p: MPa`, `T: K`, `h/u: kJ/kg`, `s: kJ/(kg·K)`, `v: m³/kg`).
- **Public string-dispatch API**: `props` / `props_si`, dispatching by function name strings (e.g. `"h_pT"`, `"T_ph"`).

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
xsteamrs = "0.1"
```

## Quick start

### XSteam-style units: `props`

- Pressure `p`: bar
- Temperature `T`: °C

```rust
use xsteamrs::props;

fn main() -> Result<(), xsteamrs::If97Error> {
    let h = props("h_pT", 1.0, 100.0)?;
    println!("h_pT(p=1 bar, T=100°C) = {h} kJ/kg");

    let t_c = props("T_ph", 1.0, 500.0)?;
    println!("T_ph(p=1 bar, h=500 kJ/kg) = {t_c} °C");

    Ok(())
}
```

### SI pressure unit: `props_si`

- Pressure `p`: Pa
- Temperature `T`: K

```rust
use xsteamrs::props_si;

fn main() -> Result<(), xsteamrs::If97Error> {
    let p_pa = 101_325.0;
    let t_k = props_si("Tsat_p", p_pa, 0.0)?;
    println!("Tsat_p(p=101325 Pa) = {t_k} K");

    let t_k = 373.15;
    let p_sat_pa = props_si("psat_T", t_k, 0.0)?;
    println!("psat_T(T=373.15 K) = {p_sat_pa} Pa");

    Ok(())
}
```

## Low-level IF97 example

If you want to work directly with IF97 conventional SI units (`p: MPa`, `T: K`), call the low-level functions:

```rust
use xsteamrs::{h1_p_t, psat_mpa_from_t_k};

fn main() -> Result<(), xsteamrs::If97Error> {
    let p_sat_mpa = psat_mpa_from_t_k(373.15)?;
    println!("psat(T=373.15 K) = {p_sat_mpa} MPa");

    let h = h1_p_t(0.101_325, 373.15)?;
    println!("h1_p_t(p=0.101325 MPa, T=373.15 K) = {h} kJ/kg");

    Ok(())
}
```

## Build and test

```bash
cargo fmt
cargo test
cargo clippy --all-targets --all-features -- -D warnings
```

Generate and open local documentation:

```bash
cargo doc --no-deps --open
```

## Documentation

After publishing to crates.io, docs are built and hosted by docs.rs:

- https://docs.rs/xsteamrs
