use xsteamrs::{If97Result, h1_p_t, props, psat_mpa_from_t_k};

fn main() -> If97Result<()> {
    let h = props("h_pT", 1.0, 100.0)?;
    println!("h_pT(p=1 bar, T=100°C) = {h} kJ/kg");

    let t_c = props("T_ph", 1.0, 100.0)?;
    println!("T_ph(p=1 bar, h=100 kJ/kg) = {t_c} °C");

    let p_sat_mpa = psat_mpa_from_t_k(373.15)?;
    println!("psat(T=373.15 K) = {p_sat_mpa} MPa");

    let h_r1 = h1_p_t(0.101_325, 373.15)?;
    println!("h1_p_t(p=0.101325 MPa, T=373.15 K) = {h_r1} kJ/kg");

    Ok(())
}
