//! IAPWS-IF97 水/水蒸气物性计算（Rust）
//!
//! 这份库提供两种调用方式：
//!
//! - **面向 IF97 的底层函数**：以 Region 为单位的 `*_p_t` / `*_rho_t` 等函数，输入输出均为 IF97 常用 SI 单位（`p: MPa`、`T: K`、`h/u: kJ/kg`、`s: kJ/(kg·K)`、`v: m³/kg`）。
//! - **对外接口**：[`props`] / [`props_si`]，通过字符串函数名（如 `"h_pT"`、`"T_ph"`）做分发。
//!   - [`props`]：单位与 XSteam.m 默认一致（`p: bar`、`T: °C`，其余大多为 IF97 SI 单位）
//!   - [`props_si`]：输入为国际单位制（`p: Pa`、`T: K`，其余同 IF97 常用 SI 单位）；且温度输出为 `K`、压力输出为 `Pa`
//!
//! # 快速上手：用 XSteam 风格接口
//!
//! - 函数名大小写不敏感，会做 `trim()` + `to_ascii_lowercase()`
//! - 不同函数名对应的输入单位不同（下面表格给出最常见的几类）
//!
//! ```no_run
//! use xsteamrs::props;
//!
//! fn rel_close(a: f64, b: f64, rel: f64) -> bool {
//!     if b == 0.0 { (a - b).abs() <= rel } else { ((a - b) / b).abs() <= rel }
//! }
//!
//! // 例1：给定压力( bar )与温度( °C )，求比焓 h (kJ/kg)
//! let h = props("h_pT", 1.0, 100.0).unwrap();
//! assert!(rel_close(h, 2_675.0, 1e-3));
//!
//! // 例2：给定压力( bar )与比焓( kJ/kg )，求温度( °C )
//! let t_c = props("T_ph", 1.0, 100.0).unwrap();
//! assert!(t_c.is_finite());
//! ```
//!
//! ## XSteam 风格接口的单位约定（常用子集）
//!
//! - `*_p*`：`p` 用 **bar**
//! - `*_T*` / `*_t*`：`T` 用 **°C**（注意：底层 IF97 函数使用 K）
//! - `h/u`：**kJ/kg**
//! - `s`：**kJ/(kg·K)**
//! - `v`：**m³/kg**
//! - `rho`：**kg/m³**
//! - `Cp/Cv`：**kJ/(kg·K)**
//! - `w`：**m/s**
//! - `my`（动力黏度）：按 XSteam.m 默认输出（通常为 `Pa·s`）
//! - `tc`（导热系数）：按 XSteam.m 默认输出（通常为 `W/(m·K)`）
//! - `st`（表面张力）：按 XSteam.m 默认输出（通常为 `N/m`）
//!
//! 支持的函数名集合与 XSteam.m 的公开接口一致，可参考
//! [XSteam.m#L3826-L3829](file:///d:/5.code/xsteamrs/XSteam.m#L3826-L3829) 的 `fun` 列表。
//!
//! # 使用底层 IF97 函数（更类型/单位明确）
//!
//! 底层函数通常以 `p_mpa` 与 `t_k` 命名参数，避免混淆单位：
//! - 如果只知道 `(p, T)`：先用 [`region_from_p_t`] 判区，再调用 `h1_p_t` / `h2_p_t` / `h3_p_t` / `h5_p_t` 等
//! - 如果只知道 `(p, h)`：先用 [`region_from_p_h`] 判区，再调用 `t1_p_h` / `t2_p_h` / `t3_p_h` / `t5_p_h` 等
//! - 如果只知道 `(p, s)`：先用 [`region_from_p_s`] 判区，再调用 `t1_p_s` / `t2_p_s` / `t3_p_s` / `t5_p_s` 等
//!
//! ```no_run
//! use xsteamrs::{h1_p_t, psat_mpa_from_t_k};
//!
//! // 饱和压力：输入 T(K)，输出 p(MPa)
//! let p_sat = psat_mpa_from_t_k(373.15).unwrap();
//! assert!(p_sat > 0.0);
//!
//! // Region1：输入 p(MPa), T(K)，输出 h(kJ/kg)
//! let h = h1_p_t(0.101_325, 373.15).unwrap();
//! assert!(h.is_finite());
//! ```
//!
//! # 错误处理
//!
//! 推荐使用 [`props`] / [`props_si`] 返回的 [`If97Result`] 显式处理错误。
//!
#[derive(Debug, Clone, Copy, PartialEq)]
/// IF97 计算过程中可能出现的错误。
pub enum If97Error {
    /// 输入为 `NaN` 或 `±Inf`。
    NonFiniteInput,
    /// 温度（K）超出 IF97 支持范围。
    TemperatureOutOfRange { t_k: f64 },
    /// 压力（MPa）超出 IF97 支持范围。
    PressureOutOfRange { p_mpa: f64 },
    /// 计算中间量无效（如超出公式适用域、迭代不收敛、溢出等）。
    InvalidIntermediateValue,
}

/// IF97 计算结果类型。
pub type If97Result<T> = Result<T, If97Error>;

fn sqrt_nonneg(x: f64) -> If97Result<f64> {
    if x.is_nan() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    if x < 0.0 {
        if x > -1e-12 {
            return Ok(0.0);
        }
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(x.sqrt())
}

/// 给定温度 `T`（K），计算饱和压力 `psat`（MPa）。
///
/// 适用范围：`273.15..=647.096 K`。
pub fn psat_mpa_from_t_k(t_k: f64) -> If97Result<f64> {
    if !t_k.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if !(273.15..=647.096).contains(&t_k) {
        return Err(If97Error::TemperatureOutOfRange { t_k });
    }

    let teta = t_k - 0.238_555_575_678_49 / (t_k - 650.175_348_447_98);
    let a = teta * teta + 1_167.052_145_276_7 * teta - 724_213.167_032_06;
    let b = -17.073_846_940_092 * teta * teta + 12_020.824_702_47 * teta - 3_232_555.032_233_3;
    let c = 14.915_108_613_53 * teta * teta - 4_823.265_736_159_1 * teta + 405_113.405_420_57;
    let disc = b * b - 4.0 * a * c;
    let sqrt_disc = sqrt_nonneg(disc)?;
    let denom = -b + sqrt_disc;
    if denom == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let p = (2.0 * c / denom).powi(4);
    if !p.is_finite() || p < 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(p)
}

/// 给定饱和压力 `p`（MPa），计算饱和温度 `Tsat`（K）。
///
/// 适用范围：约 `0.000611..=22.06395 MPa`。
pub fn tsat_k_from_p_mpa(p_mpa: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if !(0.000_611..=22.063_95).contains(&p_mpa) {
        return Err(If97Error::PressureOutOfRange { p_mpa });
    }

    let beta = p_mpa.powf(0.25);
    let e = beta * beta - 17.073_846_940_092 * beta + 14.915_108_613_53;
    let f = 1_167.052_145_276_7 * beta * beta + 12_020.824_702_47 * beta - 4_823.265_736_159_1;
    let g = -724_213.167_032_06 * beta * beta - 3_232_555.032_233_3 * beta + 405_113.405_420_57;
    let disc = f * f - 4.0 * e * g;
    let sqrt_disc = sqrt_nonneg(disc)?;
    let denom = -f - sqrt_disc;
    if denom == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let d = 2.0 * g / denom;
    let k = 650.175_348_447_98;
    let disc2 = (k + d) * (k + d) - 4.0 * (-0.238_555_575_678_49 + k * d);
    let sqrt_disc2 = sqrt_nonneg(disc2)?;
    let t = (k + d - sqrt_disc2) / 2.0;
    if !t.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(t)
}

/// XSteam 风格：给定温度 `T`（°C），计算饱和压力 `psat`（bar）。
pub fn psat_bar_from_t_c(t_c: f64) -> If97Result<f64> {
    psat_mpa_from_t_k(t_c + 273.15).map(|p_mpa| p_mpa * 10.0)
}

/// XSteam 风格：给定饱和压力 `p`（bar），计算饱和温度 `Tsat`（°C）。
pub fn tsat_c_from_p_bar(p_bar: f64) -> If97Result<f64> {
    tsat_k_from_p_mpa(p_bar / 10.0).map(|t_k| t_k - 273.15)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Region {
    R1,
    R2,
    R3,
    R4,
    R5,
}

pub fn b23_p_mpa_from_t_k(t_k: f64) -> If97Result<f64> {
    if !t_k.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let p = 348.051_856_289_69 - 1.167_185_987_997_5 * t_k + 1.019_297_003_932_6e-3 * t_k * t_k;
    if !p.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(p)
}

pub fn b23_t_k_from_p_mpa(p_mpa: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let disc = (p_mpa - 13.918_839_778_87) / 1.019_297_003_932_6e-3;
    let root = sqrt_nonneg(disc)?;
    let t = 572.544_598_627_46 + root;
    if !t.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(t)
}

pub fn region_from_p_t(p_mpa: f64, t_k: f64) -> If97Result<Region> {
    if !p_mpa.is_finite() || !t_k.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if t_k > 1073.15 && p_mpa < 10.0 && t_k < 2273.15 && p_mpa > 0.000_611 {
        return Ok(Region::R5);
    }
    if !(273.15 < t_k && t_k <= 1073.15 && 0.000_611 < p_mpa && p_mpa <= 100.0) {
        return Err(If97Error::InvalidIntermediateValue);
    }

    if t_k > 623.15 {
        if p_mpa > b23_p_mpa_from_t_k(t_k)? {
            if t_k < 647.096 {
                let ps = psat_mpa_from_t_k(t_k)?;
                if (p_mpa - ps).abs() < 1e-5 {
                    return Ok(Region::R4);
                }
            }
            Ok(Region::R3)
        } else {
            Ok(Region::R2)
        }
    } else {
        let ps = psat_mpa_from_t_k(t_k)?;
        if (p_mpa - ps).abs() < 1e-5 {
            Ok(Region::R4)
        } else if p_mpa > ps {
            Ok(Region::R1)
        } else {
            Ok(Region::R2)
        }
    }
}

const R_KJ_KG_K: f64 = 0.461_526;

fn to_si_p_mpa_from_bar(p_bar: f64) -> f64 {
    p_bar / 10.0
}

fn from_si_p_bar_from_mpa(p_mpa: f64) -> f64 {
    p_mpa * 10.0
}

fn to_si_p_pa_from_bar(p_bar: f64) -> f64 {
    p_bar * 100_000.0
}

fn from_si_p_bar_from_pa(p_pa: f64) -> f64 {
    p_pa / 100_000.0
}

fn to_si_t_k_from_c(t_c: f64) -> f64 {
    t_c + 273.15
}

fn from_si_t_c_from_k(t_k: f64) -> f64 {
    t_k - 273.15
}

fn p3sat_mpa_from_h(h: f64) -> If97Result<f64> {
    if !h.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let ii: [i32; 14] = [0, 1, 1, 1, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36];
    let ji: [i32; 14] = [0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3, 18, 8, 24];
    let ni: [f64; 14] = [
        0.600_073_641_753_024,
        -9.362_036_548_498_57,
        24.659_079_859_414_7,
        -107.014_222_858_224,
        -91_582_131_580_576.8,
        -8_623.320_117_006_62,
        -23.583_734_474_003_2,
        2.523_049_693_841_28e17,
        -3.897_187_719_977_19e18,
        -3.337_757_136_452_96e22,
        35_649_946_963.632_8,
        -1.485_475_447_206_41e26,
        3.306_115_148_387_98e18,
        8.136_412_944_678_29e37,
    ];
    let hs = h / 2600.0;
    let mut ps = 0.0;
    for i in 0..14 {
        ps += ni[i] * (hs - 1.02).powi(ii[i]) * (hs - 0.608).powi(ji[i]);
    }
    let p = ps * 22.0;
    if !p.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(p)
}

fn h4l_p_mpa(p_mpa: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if !(0.000_611_657..=22.063_95).contains(&p_mpa) {
        return Err(If97Error::PressureOutOfRange { p_mpa });
    }
    let ts = tsat_k_from_p_mpa(p_mpa)?;
    if p_mpa < 16.529 {
        return h1_p_t(p_mpa, ts);
    }
    let mut low = 1_670.858_218;
    let mut high = 2_087.235_001_648_64;
    let mut hs = (low + high) / 2.0;
    for _ in 0..80 {
        hs = (low + high) / 2.0;
        let ps = p3sat_mpa_from_h(hs)?;
        if (p_mpa - ps).abs() <= 1e-6 {
            return Ok(hs);
        }
        if ps > p_mpa {
            high = hs;
        } else {
            low = hs;
        }
    }
    Ok(hs)
}

fn h4v_p_mpa(p_mpa: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if !(0.000_611_657..=22.063_95).contains(&p_mpa) {
        return Err(If97Error::PressureOutOfRange { p_mpa });
    }
    let ts = tsat_k_from_p_mpa(p_mpa)?;
    if p_mpa < 16.529 {
        return h2_p_t(p_mpa, ts);
    }
    let mut low = 2_087.235_001_648_64;
    let mut high = 2_563.592_004 + 5.0;
    let mut hs = (low + high) / 2.0;
    for _ in 0..120 {
        hs = (low + high) / 2.0;
        let ps = p3sat_mpa_from_h(hs)?;
        if (p_mpa - ps).abs() <= 1e-6 {
            return Ok(hs);
        }
        if ps < p_mpa {
            high = hs;
        } else {
            low = hs;
        }
    }
    Ok(hs)
}

fn x4_ph(p_mpa: f64, h: f64) -> If97Result<f64> {
    let hv = h4v_p_mpa(p_mpa)?;
    let hl = h4l_p_mpa(p_mpa)?;
    if h > hv {
        Ok(1.0)
    } else if h < hl {
        Ok(0.0)
    } else {
        Ok((h - hl) / (hv - hl))
    }
}

fn x4_ps(p_mpa: f64, s: f64) -> If97Result<f64> {
    let ts = tsat_k_from_p_mpa(p_mpa)?;
    let (sv, sl) = if p_mpa < 16.529 {
        (s2_p_t(p_mpa, ts)?, s1_p_t(p_mpa, ts)?)
    } else {
        let rho_v = 1.0 / v3_p_h(p_mpa, h4v_p_mpa(p_mpa)?)?;
        let rho_l = 1.0 / v3_p_h(p_mpa, h4l_p_mpa(p_mpa)?)?;
        (s3_rho_t(rho_v, ts)?, s3_rho_t(rho_l, ts)?)
    };
    if s < sl {
        Ok(0.0)
    } else if s > sv {
        Ok(1.0)
    } else {
        Ok((s - sl) / (sv - sl))
    }
}

fn v4l_p_mpa(p_mpa: f64) -> If97Result<f64> {
    let ts = tsat_k_from_p_mpa(p_mpa)?;
    if p_mpa < 16.529 {
        v1_p_t(p_mpa, ts)
    } else {
        v3_p_h(p_mpa, h4l_p_mpa(p_mpa)?)
    }
}

fn v4v_p_mpa(p_mpa: f64) -> If97Result<f64> {
    let ts = tsat_k_from_p_mpa(p_mpa)?;
    if p_mpa < 16.529 {
        v2_p_t(p_mpa, ts)
    } else {
        v3_p_h(p_mpa, h4v_p_mpa(p_mpa)?)
    }
}

fn p4_h_s(h: f64, s: f64) -> If97Result<f64> {
    if !h.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let mut low = 0.000_611_657;
    let mut high = 22.063_95;
    let mut p = (low + high) / 2.0;
    for _ in 0..120 {
        p = (low + high) / 2.0;
        let x = x4_ph(p, h)?;
        let ts = tsat_k_from_p_mpa(p)?;
        let (sv, sl) = if p < 16.529 {
            (s2_p_t(p, ts)?, s1_p_t(p, ts)?)
        } else {
            let rho_v = 1.0 / v3_p_h(p, h4v_p_mpa(p)?)?;
            let rho_l = 1.0 / v3_p_h(p, h4l_p_mpa(p)?)?;
            (s3_rho_t(rho_v, ts)?, s3_rho_t(rho_l, ts)?)
        };
        let s_mix = x * sv + (1.0 - x) * sl;
        if (s - s_mix).abs() <= 1e-6 && (high - low) <= 1e-7 {
            return Ok(p);
        }
        if s_mix < s {
            high = p;
        } else {
            low = p;
        }
    }
    Ok(p)
}

fn t4_h_s(h: f64, s: f64) -> If97Result<f64> {
    let p = p4_h_s(h, s)?;
    tsat_k_from_p_mpa(p)
}

fn region_from_p_h(p_mpa: f64, h: f64) -> If97Result<Region> {
    if !p_mpa.is_finite() || !h.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if !(0.000_611_657..=100.0).contains(&p_mpa) {
        return Err(If97Error::PressureOutOfRange { p_mpa });
    }
    if p_mpa <= 22.063_95 {
        let hl = h4l_p_mpa(p_mpa).ok();
        let hv = h4v_p_mpa(p_mpa).ok();
        if let (Some(hl), Some(hv)) = (hl, hv)
            && h >= hl
            && h <= hv
        {
            return Ok(Region::R4);
        }
    }
    if let Ok(t) = t1_p_h(p_mpa, h)
        && region_from_p_t(p_mpa, t).ok() == Some(Region::R1)
    {
        return Ok(Region::R1);
    }
    if let Ok(t) = t2_p_h(p_mpa, h)
        && region_from_p_t(p_mpa, t).ok() == Some(Region::R2)
    {
        return Ok(Region::R2);
    }
    if let Ok(t) = t3_p_h(p_mpa, h)
        && region_from_p_t(p_mpa, t).ok() == Some(Region::R3)
    {
        return Ok(Region::R3);
    }
    if p_mpa <= 10.0
        && let Ok(t) = t5_p_h(p_mpa, h)
        && region_from_p_t(p_mpa, t).ok() == Some(Region::R5)
    {
        return Ok(Region::R5);
    }
    Err(If97Error::InvalidIntermediateValue)
}

fn region_from_p_s(p_mpa: f64, s: f64) -> If97Result<Region> {
    if !p_mpa.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if !(0.000_611_657..=100.0).contains(&p_mpa) {
        return Err(If97Error::PressureOutOfRange { p_mpa });
    }
    if p_mpa <= 22.063_95 {
        let ts = tsat_k_from_p_mpa(p_mpa).ok();
        if let Ok(ts) = ts.ok_or(If97Error::InvalidIntermediateValue) {
            let (sv, sl) = if p_mpa < 16.529 {
                (s2_p_t(p_mpa, ts)?, s1_p_t(p_mpa, ts)?)
            } else {
                let rho_v = 1.0 / v3_p_h(p_mpa, h4v_p_mpa(p_mpa)?)?;
                let rho_l = 1.0 / v3_p_h(p_mpa, h4l_p_mpa(p_mpa)?)?;
                (s3_rho_t(rho_v, ts)?, s3_rho_t(rho_l, ts)?)
            };
            if s >= sl && s <= sv {
                return Ok(Region::R4);
            }
        }
    }
    if let Ok(t) = t1_p_s(p_mpa, s)
        && region_from_p_t(p_mpa, t).ok() == Some(Region::R1)
    {
        return Ok(Region::R1);
    }
    if let Ok(t) = t2_p_s(p_mpa, s)
        && region_from_p_t(p_mpa, t).ok() == Some(Region::R2)
    {
        return Ok(Region::R2);
    }
    if let Ok(t) = t3_p_s(p_mpa, s)
        && region_from_p_t(p_mpa, t).ok() == Some(Region::R3)
    {
        return Ok(Region::R3);
    }
    if p_mpa <= 10.0
        && let Ok(t) = t5_p_s(p_mpa, s)
        && region_from_p_t(p_mpa, t).ok() == Some(Region::R5)
    {
        return Ok(Region::R5);
    }
    Err(If97Error::InvalidIntermediateValue)
}

fn p_from_h_s(h: f64, s: f64) -> If97Result<f64> {
    if let Ok(p) = p1_h_s(h, s)
        && let Ok(t) = t1_p_h(p, h)
        && region_from_p_t(p, t).ok() == Some(Region::R1)
    {
        return Ok(p);
    }
    if let Ok(p) = p2_h_s(h, s)
        && let Ok(t) = t2_p_h(p, h)
        && region_from_p_t(p, t).ok() == Some(Region::R2)
    {
        return Ok(p);
    }
    if let Ok(p) = p3_h_s(h, s)
        && let Ok(t) = t3_p_h(p, h)
        && region_from_p_t(p, t).ok() == Some(Region::R3)
    {
        return Ok(p);
    }
    let p = p4_h_s(h, s)?;
    Ok(p)
}

const R5_J0: [i32; 6] = [0, 1, -3, -2, -1, 2];
const R5_N0: [f64; 6] = [
    -13.179_983_674_201,
    6.854_084_163_443_4,
    -0.024_805_148_933_466,
    0.369_015_349_803_33,
    -3.116_131_821_392_5,
    -0.329_616_265_389_17,
];
const R5_IR: [i32; 5] = [1, 1, 1, 2, 3];
const R5_JR: [i32; 5] = [0, 1, 3, 9, 3];
const R5_NR: [f64; 5] = [
    -1.256_318_358_959_2e-4,
    2.177_467_871_457_1e-3,
    -0.004_594_282_089_991,
    -3.972_482_835_956_9e-6,
    1.291_922_828_978_4e-7,
];

fn r5_common(p_mpa: f64, t_k: f64) -> If97Result<(f64, f64)> {
    if !p_mpa.is_finite() || !t_k.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if p_mpa <= 0.0 || t_k <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let pi = p_mpa;
    let tau = 1000.0 / t_k;
    Ok((pi, tau))
}

pub fn v5_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r5_common(p_mpa, t_k)?;
    let gamma0_pi = 1.0 / pi;
    let mut gammar_pi = 0.0;
    for i in 0..5 {
        gammar_pi += R5_NR[i] * (R5_IR[i] as f64) * pi.powi(R5_IR[i] - 1) * tau.powi(R5_JR[i]);
    }
    let v = R_KJ_KG_K * t_k / p_mpa * pi * (gamma0_pi + gammar_pi) / 1000.0;
    if !v.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(v)
}

pub fn h5_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r5_common(p_mpa, t_k)?;
    let mut gamma0_tau = 0.0;
    for i in 0..6 {
        gamma0_tau += R5_N0[i] * (R5_J0[i] as f64) * tau.powi(R5_J0[i] - 1);
    }
    let mut gammar_tau = 0.0;
    for i in 0..5 {
        gammar_tau += R5_NR[i] * pi.powi(R5_IR[i]) * (R5_JR[i] as f64) * tau.powi(R5_JR[i] - 1);
    }
    let h = R_KJ_KG_K * t_k * tau * (gamma0_tau + gammar_tau);
    if !h.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(h)
}

pub fn u5_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r5_common(p_mpa, t_k)?;
    let gamma0_pi = 1.0 / pi;
    let mut gamma0_tau = 0.0;
    for i in 0..6 {
        gamma0_tau += R5_N0[i] * (R5_J0[i] as f64) * tau.powi(R5_J0[i] - 1);
    }
    let mut gammar_pi = 0.0;
    let mut gammar_tau = 0.0;
    for i in 0..5 {
        gammar_pi += R5_NR[i] * (R5_IR[i] as f64) * pi.powi(R5_IR[i] - 1) * tau.powi(R5_JR[i]);
        gammar_tau += R5_NR[i] * pi.powi(R5_IR[i]) * (R5_JR[i] as f64) * tau.powi(R5_JR[i] - 1);
    }
    let u = R_KJ_KG_K * t_k * (tau * (gamma0_tau + gammar_tau) - pi * (gamma0_pi + gammar_pi));
    if !u.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(u)
}

pub fn s5_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r5_common(p_mpa, t_k)?;
    let mut gamma0 = pi.ln();
    let mut gamma0_tau = 0.0;
    for i in 0..6 {
        gamma0 += R5_N0[i] * tau.powi(R5_J0[i]);
        gamma0_tau += R5_N0[i] * (R5_J0[i] as f64) * tau.powi(R5_J0[i] - 1);
    }
    let mut gammar = 0.0;
    let mut gammar_tau = 0.0;
    for i in 0..5 {
        gammar += R5_NR[i] * pi.powi(R5_IR[i]) * tau.powi(R5_JR[i]);
        gammar_tau += R5_NR[i] * pi.powi(R5_IR[i]) * (R5_JR[i] as f64) * tau.powi(R5_JR[i] - 1);
    }
    let s = R_KJ_KG_K * (tau * (gamma0_tau + gammar_tau) - (gamma0 + gammar));
    if !s.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(s)
}

pub fn cp5_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r5_common(p_mpa, t_k)?;
    let mut gamma0_tautau = 0.0;
    for i in 0..6 {
        gamma0_tautau +=
            R5_N0[i] * (R5_J0[i] as f64) * ((R5_J0[i] - 1) as f64) * tau.powi(R5_J0[i] - 2);
    }
    let mut gammar_tautau = 0.0;
    for i in 0..5 {
        gammar_tautau += R5_NR[i]
            * pi.powi(R5_IR[i])
            * (R5_JR[i] as f64)
            * ((R5_JR[i] - 1) as f64)
            * tau.powi(R5_JR[i] - 2);
    }
    let cp = -R_KJ_KG_K * tau * tau * (gamma0_tautau + gammar_tautau);
    if !cp.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(cp)
}

pub fn cv5_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r5_common(p_mpa, t_k)?;
    let mut gamma0_tautau = 0.0;
    for i in 0..6 {
        gamma0_tautau +=
            R5_N0[i] * (R5_J0[i] as f64) * ((R5_J0[i] - 1) as f64) * tau.powi(R5_J0[i] - 2);
    }
    let mut gammar_pi = 0.0;
    let mut gammar_pitau = 0.0;
    let mut gammar_pipi = 0.0;
    let mut gammar_tautau = 0.0;
    for i in 0..5 {
        gammar_pi += R5_NR[i] * (R5_IR[i] as f64) * pi.powi(R5_IR[i] - 1) * tau.powi(R5_JR[i]);
        gammar_pipi += R5_NR[i]
            * (R5_IR[i] as f64)
            * ((R5_IR[i] - 1) as f64)
            * pi.powi(R5_IR[i] - 2)
            * tau.powi(R5_JR[i]);
        gammar_pitau += R5_NR[i]
            * (R5_IR[i] as f64)
            * pi.powi(R5_IR[i] - 1)
            * (R5_JR[i] as f64)
            * tau.powi(R5_JR[i] - 1);
        gammar_tautau += R5_NR[i]
            * pi.powi(R5_IR[i])
            * (R5_JR[i] as f64)
            * ((R5_JR[i] - 1) as f64)
            * tau.powi(R5_JR[i] - 2);
    }
    let denom = 1.0 - pi * pi * gammar_pipi;
    if denom == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let cv = R_KJ_KG_K
        * (-(tau * tau * (gamma0_tautau + gammar_tautau))
            - (1.0 + pi * gammar_pi - tau * pi * gammar_pitau).powi(2) / denom);
    if !cv.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(cv)
}

pub fn w5_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r5_common(p_mpa, t_k)?;
    let mut gamma0_tautau = 0.0;
    for i in 0..6 {
        gamma0_tautau +=
            R5_N0[i] * (R5_J0[i] as f64) * ((R5_J0[i] - 1) as f64) * tau.powi(R5_J0[i] - 2);
    }
    let mut gammar_pi = 0.0;
    let mut gammar_pitau = 0.0;
    let mut gammar_pipi = 0.0;
    let mut gammar_tautau = 0.0;
    for i in 0..5 {
        gammar_pi += R5_NR[i] * (R5_IR[i] as f64) * pi.powi(R5_IR[i] - 1) * tau.powi(R5_JR[i]);
        gammar_pipi += R5_NR[i]
            * (R5_IR[i] as f64)
            * ((R5_IR[i] - 1) as f64)
            * pi.powi(R5_IR[i] - 2)
            * tau.powi(R5_JR[i]);
        gammar_pitau += R5_NR[i]
            * (R5_IR[i] as f64)
            * pi.powi(R5_IR[i] - 1)
            * (R5_JR[i] as f64)
            * tau.powi(R5_JR[i] - 1);
        gammar_tautau += R5_NR[i]
            * pi.powi(R5_IR[i])
            * (R5_JR[i] as f64)
            * ((R5_JR[i] - 1) as f64)
            * tau.powi(R5_JR[i] - 2);
    }
    let term1 = 1.0 + 2.0 * pi * gammar_pi + pi * pi * gammar_pi * gammar_pi;
    let denom = (1.0 - pi * pi * gammar_pipi)
        + (1.0 + pi * gammar_pi - tau * pi * gammar_pitau).powi(2)
            / (tau * tau * (gamma0_tautau + gammar_tautau));
    if denom == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let w2 = 1000.0 * R_KJ_KG_K * t_k * term1 / denom;
    let w = sqrt_nonneg(w2)?;
    if !w.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(w)
}

pub fn t5_p_h(p_mpa: f64, h: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !h.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let mut low = 1073.15;
    let mut high = 2273.15;
    let mut t = (low + high) / 2.0;
    for _ in 0..80 {
        t = (low + high) / 2.0;
        let hs = h5_p_t(p_mpa, t)?;
        if (hs - h).abs() <= 1e-5 {
            return Ok(t);
        }
        if hs > h {
            high = t;
        } else {
            low = t;
        }
    }
    Ok(t)
}

pub fn t5_p_s(p_mpa: f64, s: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let mut low = 1073.15;
    let mut high = 2273.15;
    let mut t = (low + high) / 2.0;
    for _ in 0..100 {
        t = (low + high) / 2.0;
        let ss = s5_p_t(p_mpa, t)?;
        if (ss - s).abs() <= 1e-5 {
            return Ok(t);
        }
        if ss > s {
            high = t;
        } else {
            low = t;
        }
    }
    Ok(t)
}

fn my_allregions_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let rho = match region_from_p_t(p_mpa, t_k)? {
        Region::R1 => 1.0 / v1_p_t(p_mpa, t_k)?,
        Region::R2 => 1.0 / v2_p_t(p_mpa, t_k)?,
        Region::R3 => rho3_p_t(p_mpa, t_k)?,
        Region::R5 => 1.0 / v5_p_t(p_mpa, t_k)?,
        Region::R4 => return Err(If97Error::InvalidIntermediateValue),
    };
    let rhos = rho / 317.763;
    let ts = t_k / 647.226;
    if t_k > 1173.15
        || (t_k > 873.15 && p_mpa > 300.0)
        || (t_k > 423.15 && p_mpa > 350.0)
        || p_mpa > 500.0
    {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let my0 =
        ts.sqrt() / (1.0 + 0.978_197 / ts + 0.579_829 / (ts * ts) - 0.202_354 / (ts * ts * ts));
    let h0: [f64; 6] = [
        0.513_204_7,
        0.320_565_6,
        0.0,
        0.0,
        -0.778_256_7,
        0.188_544_7,
    ];
    let h1: [f64; 6] = [0.215_177_8, 0.731_788_3, 1.241_044, 1.476_783, 0.0, 0.0];
    let h2: [f64; 6] = [-0.281_810_7, -1.070_786, -1.263_184, 0.0, 0.0, 0.0];
    let h3: [f64; 6] = [0.177_806_4, 0.460_504, 0.234_037_9, -0.492_417_9, 0.0, 0.0];
    let h4: [f64; 6] = [-0.041_766_1, 0.0, 0.0, 0.160_043_5, 0.0, 0.0];
    let h5: [f64; 6] = [0.0, -0.015_783_86, 0.0, 0.0, 0.0, 0.0];
    let h6: [f64; 6] = [0.0, 0.0, 0.0, -0.003_629_481, 0.0, 0.0];
    let inv_ts_minus1 = 1.0 / ts - 1.0;
    let dr = rhos - 1.0;
    let mut sum = 0.0;
    for i in 0..6 {
        let a = inv_ts_minus1.powi(i as i32);
        sum += h0[i] * a;
        sum += h1[i] * a * dr;
        sum += h2[i] * a * dr.powi(2);
        sum += h3[i] * a * dr.powi(3);
        sum += h4[i] * a * dr.powi(4);
        sum += h5[i] * a * dr.powi(5);
        sum += h6[i] * a * dr.powi(6);
    }
    let my1 = (rhos * sum).exp();
    let mys = my0 * my1;
    let my = mys * 0.000_055_071;
    if !my.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(my)
}

fn my_allregions_p_h(p_mpa: f64, h: f64) -> If97Result<f64> {
    let region = region_from_p_h(p_mpa, h)?;
    if region == Region::R4 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let (t_k, rho) = match region {
        Region::R1 => {
            let t = t1_p_h(p_mpa, h)?;
            (t, 1.0 / v1_p_t(p_mpa, t)?)
        }
        Region::R2 => {
            let t = t2_p_h(p_mpa, h)?;
            (t, 1.0 / v2_p_t(p_mpa, t)?)
        }
        Region::R3 => {
            let t = t3_p_h(p_mpa, h)?;
            (t, 1.0 / v3_p_h(p_mpa, h)?)
        }
        Region::R5 => {
            let t = t5_p_h(p_mpa, h)?;
            (t, 1.0 / v5_p_t(p_mpa, t)?)
        }
        Region::R4 => unreachable!(),
    };
    let rhos = rho / 317.763;
    let ts = t_k / 647.226;
    if t_k > 1173.15
        || (t_k > 873.15 && p_mpa > 300.0)
        || (t_k > 423.15 && p_mpa > 350.0)
        || p_mpa > 500.0
    {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let my0 =
        ts.sqrt() / (1.0 + 0.978_197 / ts + 0.579_829 / (ts * ts) - 0.202_354 / (ts * ts * ts));
    let h0: [f64; 6] = [
        0.513_204_7,
        0.320_565_6,
        0.0,
        0.0,
        -0.778_256_7,
        0.188_544_7,
    ];
    let h1: [f64; 6] = [0.215_177_8, 0.731_788_3, 1.241_044, 1.476_783, 0.0, 0.0];
    let h2: [f64; 6] = [-0.281_810_7, -1.070_786, -1.263_184, 0.0, 0.0, 0.0];
    let h3: [f64; 6] = [0.177_806_4, 0.460_504, 0.234_037_9, -0.492_417_9, 0.0, 0.0];
    let h4: [f64; 6] = [-0.041_766_1, 0.0, 0.0, 0.160_043_5, 0.0, 0.0];
    let h5: [f64; 6] = [0.0, -0.015_783_86, 0.0, 0.0, 0.0, 0.0];
    let h6: [f64; 6] = [0.0, 0.0, 0.0, -0.003_629_481, 0.0, 0.0];
    let inv_ts_minus1 = 1.0 / ts - 1.0;
    let dr = rhos - 1.0;
    let mut sum = 0.0;
    for i in 0..6 {
        let a = inv_ts_minus1.powi(i as i32);
        sum += h0[i] * a;
        sum += h1[i] * a * dr;
        sum += h2[i] * a * dr.powi(2);
        sum += h3[i] * a * dr.powi(3);
        sum += h4[i] * a * dr.powi(4);
        sum += h5[i] * a * dr.powi(5);
        sum += h6[i] * a * dr.powi(6);
    }
    let my1 = (rhos * sum).exp();
    let mys = my0 * my1;
    let my = mys * 0.000_055_071;
    if !my.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(my)
}

fn tc_ptrho(p_mpa: f64, t_k: f64, rho: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !t_k.is_finite() || !rho.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if t_k < 273.15 {
        return Err(If97Error::TemperatureOutOfRange { t_k });
    }
    if t_k < 773.15 {
        if p_mpa > 100.0 {
            return Err(If97Error::PressureOutOfRange { p_mpa });
        }
    } else if t_k <= 923.15 {
        if p_mpa > 70.0 {
            return Err(If97Error::PressureOutOfRange { p_mpa });
        }
    } else if p_mpa > 40.0 {
        return Err(If97Error::PressureOutOfRange { p_mpa });
    }
    let t = t_k / 647.26;
    let rho_r = rho / 317.7;
    let tc0 =
        t.sqrt() * (0.010_281_1 + 0.029_962_1 * t + 0.015_614_6 * t * t - 0.004_224_64 * t * t * t);
    let tc1 =
        -0.397_07 + 0.400_302 * rho_r + 1.06 * (-0.171_587 * (rho_r + 2.392_19).powi(2)).exp();
    let dt = (t - 1.0).abs() + 0.003_089_76;
    let q = 2.0 + 0.082_299_4 / dt.powf(3.0 / 5.0);
    let s = if t >= 1.0 {
        1.0 / dt
    } else {
        10.093_2 / dt.powf(3.0 / 5.0)
    };
    let tc2 = (0.070_130_9 / t.powi(10) + 0.011_852)
        * rho_r.powf(9.0 / 5.0)
        * (0.642_857 * (1.0 - rho_r.powf(14.0 / 5.0))).exp()
        + 0.001_699_37 * s * rho_r.powf(q) * ((q / (1.0 + q)) * (1.0 - rho_r.powf(1.0 + q))).exp()
        - 1.02 * (-4.117_17 * t.powf(3.0 / 2.0) - 6.179_37 / rho_r.powi(5)).exp();
    let tc = tc0 + tc1 + tc2;
    if !tc.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(tc)
}

fn surface_tension_t(t_k: f64) -> If97Result<f64> {
    if !t_k.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let tc = 647.096;
    if t_k < 0.01 || t_k > tc {
        return Err(If97Error::TemperatureOutOfRange { t_k });
    }
    let b = 0.2358;
    let bb = -0.625;
    let my = 1.256;
    let tau = 1.0 - t_k / tc;
    let st = b * tau.powf(my) * (1.0 + bb * tau);
    if !st.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(st)
}

/// 对外接口（XSteam 单位：bar / °C）。
///
/// - `fun`：函数名（大小写不敏感，内部会 `trim()` + `to_ascii_lowercase()`）
/// - `in1`/`in2`：与 `fun` 对应的两个输入参数（单位遵循 XSteam.m 约定）
/// - 返回值：成功时为对应物性值；失败时返回 [`If97Error`]
///
/// 常见单位约定：
/// - `p`：bar
/// - `T`：°C
/// - `h/u`：kJ/kg
/// - `s`：kJ/(kg·K)
///
/// 支持的函数名集合可参考
/// [XSteam.m#L3826-L3829](file:///d:/5.code/xsteamrs/XSteam.m#L3826-L3829)。
pub fn props(fun: &str, in1: f64, in2: f64) -> If97Result<f64> {
    let f = fun.trim().to_ascii_lowercase();
    match f.as_str() {
        "tsat_p" => Ok(from_si_t_c_from_k(tsat_k_from_p_mpa(
            to_si_p_mpa_from_bar(in1),
        )?)),
        "t_ph" => {
            let p = to_si_p_mpa_from_bar(in1);
            let h = in2;
            match region_from_p_h(p, h)? {
                Region::R1 => Ok(from_si_t_c_from_k(t1_p_h(p, h)?)),
                Region::R2 => Ok(from_si_t_c_from_k(t2_p_h(p, h)?)),
                Region::R3 => Ok(from_si_t_c_from_k(t3_p_h(p, h)?)),
                Region::R4 => Ok(from_si_t_c_from_k(tsat_k_from_p_mpa(p)?)),
                Region::R5 => Ok(from_si_t_c_from_k(t5_p_h(p, h)?)),
            }
        }
        "t_ps" => {
            let p = to_si_p_mpa_from_bar(in1);
            let s = in2;
            match region_from_p_s(p, s)? {
                Region::R1 => Ok(from_si_t_c_from_k(t1_p_s(p, s)?)),
                Region::R2 => Ok(from_si_t_c_from_k(t2_p_s(p, s)?)),
                Region::R3 => Ok(from_si_t_c_from_k(t3_p_s(p, s)?)),
                Region::R4 => Ok(from_si_t_c_from_k(tsat_k_from_p_mpa(p)?)),
                Region::R5 => Ok(from_si_t_c_from_k(t5_p_s(p, s)?)),
            }
        }
        "t_hs" => {
            let h = in1;
            let s = in2;
            let p = p_from_h_s(h, s)?;
            let r = region_from_p_h(p, h)?;
            match r {
                Region::R1 => Ok(from_si_t_c_from_k(t1_p_h(p, h)?)),
                Region::R2 => Ok(from_si_t_c_from_k(t2_p_h(p, h)?)),
                Region::R3 => Ok(from_si_t_c_from_k(t3_p_h(p, h)?)),
                Region::R4 => Ok(from_si_t_c_from_k(t4_h_s(h, s)?)),
                Region::R5 => Ok(from_si_t_c_from_k(t5_p_h(p, h)?)),
            }
        }
        "psat_t" => Ok(from_si_p_bar_from_mpa(psat_mpa_from_t_k(
            to_si_t_k_from_c(in1),
        )?)),
        "p_hs" => Ok(from_si_p_bar_from_mpa(p_from_h_s(in1, in2)?)),
        "hv_p" => Ok(h4v_p_mpa(to_si_p_mpa_from_bar(in1))?),
        "hl_p" => Ok(h4l_p_mpa(to_si_p_mpa_from_bar(in1))?),
        "hv_t" => {
            let p = psat_mpa_from_t_k(to_si_t_k_from_c(in1))?;
            Ok(h4v_p_mpa(p)?)
        }
        "hl_t" => {
            let p = psat_mpa_from_t_k(to_si_t_k_from_c(in1))?;
            Ok(h4l_p_mpa(p)?)
        }
        "h_pt" => {
            let p = to_si_p_mpa_from_bar(in1);
            let t = to_si_t_k_from_c(in2);
            match region_from_p_t(p, t)? {
                Region::R1 => h1_p_t(p, t),
                Region::R2 => h2_p_t(p, t),
                Region::R3 => h3_p_t(p, t),
                Region::R5 => h5_p_t(p, t),
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
            }
        }
        "h_ps" => {
            let p = to_si_p_mpa_from_bar(in1);
            let s = in2;
            match region_from_p_s(p, s)? {
                Region::R1 => h1_p_t(p, t1_p_s(p, s)?),
                Region::R2 => h2_p_t(p, t2_p_s(p, s)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_s(p, s)?;
                    let t = t3_p_s(p, s)?;
                    h3_rho_t(rho, t)
                }
                Region::R4 => {
                    let x = x4_ps(p, s)?;
                    let hl = h4l_p_mpa(p)?;
                    let hv = h4v_p_mpa(p)?;
                    Ok(x * hv + (1.0 - x) * hl)
                }
                Region::R5 => h5_p_t(p, t5_p_s(p, s)?),
            }
        }
        "h_px" => {
            let p = to_si_p_mpa_from_bar(in1);
            let x = in2;
            if !(0.0..=1.0).contains(&x) || p >= 22.064 {
                return Err(If97Error::InvalidIntermediateValue);
            }
            let hl = h4l_p_mpa(p)?;
            let hv = h4v_p_mpa(p)?;
            Ok(hl + x * (hv - hl))
        }
        "h_prho" => {
            let p = to_si_p_mpa_from_bar(in1);
            let rho = in2;
            let v = 1.0 / rho;
            if p <= 22.063_95 {
                let v_l = v4l_p_mpa(p)?;
                let v_v = v4v_p_mpa(p)?;
                if v >= v_l && v <= v_v {
                    let hl = h4l_p_mpa(p)?;
                    let hv = h4v_p_mpa(p)?;
                    let x = (v - v_l) / (v_v - v_l);
                    return Ok((1.0 - x) * hl + x * hv);
                }
            }
            Err(If97Error::InvalidIntermediateValue)
        }
        "h_tx" => {
            let t = to_si_t_k_from_c(in1);
            let x = in2;
            if !(0.0..=1.0).contains(&x) {
                return Err(If97Error::InvalidIntermediateValue);
            }
            let p = psat_mpa_from_t_k(t)?;
            let hl = h4l_p_mpa(p)?;
            let hv = h4v_p_mpa(p)?;
            Ok(hl + x * (hv - hl))
        }
        "vv_p" => Ok(v4v_p_mpa(to_si_p_mpa_from_bar(in1))?),
        "vl_p" => Ok(v4l_p_mpa(to_si_p_mpa_from_bar(in1))?),
        "vv_t" => Ok(v4v_p_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?)?),
        "vl_t" => Ok(v4l_p_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?)?),
        "v_pt" => {
            let p = to_si_p_mpa_from_bar(in1);
            let t = to_si_t_k_from_c(in2);
            match region_from_p_t(p, t)? {
                Region::R1 => v1_p_t(p, t),
                Region::R2 => v2_p_t(p, t),
                Region::R3 => {
                    let rho = rho3_p_t(p, t)?;
                    Ok(1.0 / rho)
                }
                Region::R5 => v5_p_t(p, t),
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
            }
        }
        "v_ph" => {
            let p = to_si_p_mpa_from_bar(in1);
            let h = in2;
            match region_from_p_h(p, h)? {
                Region::R1 => v1_p_t(p, t1_p_h(p, h)?),
                Region::R2 => v2_p_t(p, t2_p_h(p, h)?),
                Region::R3 => v3_p_h(p, h),
                Region::R4 => {
                    let x = x4_ph(p, h)?;
                    let v_l = v4l_p_mpa(p)?;
                    let v_v = v4v_p_mpa(p)?;
                    Ok((1.0 - x) * v_l + x * v_v)
                }
                Region::R5 => v5_p_t(p, t5_p_h(p, h)?),
            }
        }
        "v_ps" => {
            let p = to_si_p_mpa_from_bar(in1);
            let s = in2;
            match region_from_p_s(p, s)? {
                Region::R1 => v1_p_t(p, t1_p_s(p, s)?),
                Region::R2 => v2_p_t(p, t2_p_s(p, s)?),
                Region::R3 => v3_p_s(p, s),
                Region::R4 => {
                    let x = x4_ps(p, s)?;
                    let v_l = v4l_p_mpa(p)?;
                    let v_v = v4v_p_mpa(p)?;
                    Ok((1.0 - x) * v_l + x * v_v)
                }
                Region::R5 => v5_p_t(p, t5_p_s(p, s)?),
            }
        }
        "rhov_p" => Ok(1.0 / v4v_p_mpa(to_si_p_mpa_from_bar(in1))?),
        "rhol_p" => Ok(1.0 / v4l_p_mpa(to_si_p_mpa_from_bar(in1))?),
        "rhov_t" => Ok(1.0 / v4v_p_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?)?),
        "rhol_t" => Ok(1.0 / v4l_p_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?)?),
        "rho_pt" => Ok(1.0 / props("v_pT", in1, in2)?),
        "rho_ph" => Ok(1.0 / props("v_ph", in1, in2)?),
        "rho_ps" => Ok(1.0 / props("v_ps", in1, in2)?),
        "sv_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                s2_p_t(p, ts)
            } else {
                let rho_v = 1.0 / v3_p_h(p, h4v_p_mpa(p)?)?;
                s3_rho_t(rho_v, ts)
            }
        }
        "sl_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                s1_p_t(p, ts)
            } else {
                let rho_l = 1.0 / v3_p_h(p, h4l_p_mpa(p)?)?;
                s3_rho_t(rho_l, ts)
            }
        }
        "sv_t" => props(
            "sV_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "sl_t" => props(
            "sL_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "s_pt" => {
            let p = to_si_p_mpa_from_bar(in1);
            let t = to_si_t_k_from_c(in2);
            match region_from_p_t(p, t)? {
                Region::R1 => s1_p_t(p, t),
                Region::R2 => s2_p_t(p, t),
                Region::R3 => {
                    let rho = rho3_p_t(p, t)?;
                    s3_rho_t(rho, t)
                }
                Region::R5 => s5_p_t(p, t),
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
            }
        }
        "s_ph" => {
            let p = to_si_p_mpa_from_bar(in1);
            let h = in2;
            match region_from_p_h(p, h)? {
                Region::R1 => s1_p_t(p, t1_p_h(p, h)?),
                Region::R2 => s2_p_t(p, t2_p_h(p, h)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_h(p, h)?;
                    let t = t3_p_h(p, h)?;
                    s3_rho_t(rho, t)
                }
                Region::R4 => {
                    let x = x4_ph(p, h)?;
                    let ts = tsat_k_from_p_mpa(p)?;
                    let (sv, sl) = if p < 16.529 {
                        (s2_p_t(p, ts)?, s1_p_t(p, ts)?)
                    } else {
                        let rho_v = 1.0 / v3_p_h(p, h4v_p_mpa(p)?)?;
                        let rho_l = 1.0 / v3_p_h(p, h4l_p_mpa(p)?)?;
                        (s3_rho_t(rho_v, ts)?, s3_rho_t(rho_l, ts)?)
                    };
                    Ok(x * sv + (1.0 - x) * sl)
                }
                Region::R5 => s5_p_t(p, t5_p_h(p, h)?),
            }
        }
        "uv_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                u2_p_t(p, ts)
            } else {
                let rho_v = 1.0 / v3_p_h(p, h4v_p_mpa(p)?)?;
                u3_rho_t(rho_v, ts)
            }
        }
        "ul_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                u1_p_t(p, ts)
            } else {
                let rho_l = 1.0 / v3_p_h(p, h4l_p_mpa(p)?)?;
                u3_rho_t(rho_l, ts)
            }
        }
        "uv_t" => props(
            "uV_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "ul_t" => props(
            "uL_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "u_pt" => {
            let p = to_si_p_mpa_from_bar(in1);
            let t = to_si_t_k_from_c(in2);
            match region_from_p_t(p, t)? {
                Region::R1 => u1_p_t(p, t),
                Region::R2 => u2_p_t(p, t),
                Region::R3 => {
                    let rho = rho3_p_t(p, t)?;
                    u3_rho_t(rho, t)
                }
                Region::R5 => u5_p_t(p, t),
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
            }
        }
        "u_ph" => {
            let p = to_si_p_mpa_from_bar(in1);
            let h = in2;
            match region_from_p_h(p, h)? {
                Region::R1 => u1_p_t(p, t1_p_h(p, h)?),
                Region::R2 => u2_p_t(p, t2_p_h(p, h)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_h(p, h)?;
                    let t = t3_p_h(p, h)?;
                    u3_rho_t(rho, t)
                }
                Region::R4 => {
                    let x = x4_ph(p, h)?;
                    let ts = tsat_k_from_p_mpa(p)?;
                    let (uv, ul) = if p < 16.529 {
                        (u2_p_t(p, ts)?, u1_p_t(p, ts)?)
                    } else {
                        let rho_v = 1.0 / v3_p_h(p, h4v_p_mpa(p)?)?;
                        let rho_l = 1.0 / v3_p_h(p, h4l_p_mpa(p)?)?;
                        (u3_rho_t(rho_v, ts)?, u3_rho_t(rho_l, ts)?)
                    };
                    Ok(x * uv + (1.0 - x) * ul)
                }
                Region::R5 => u5_p_t(p, t5_p_h(p, h)?),
            }
        }
        "u_ps" => {
            let p = to_si_p_mpa_from_bar(in1);
            let s = in2;
            match region_from_p_s(p, s)? {
                Region::R1 => u1_p_t(p, t1_p_s(p, s)?),
                Region::R2 => u2_p_t(p, t2_p_s(p, s)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_s(p, s)?;
                    let t = t3_p_s(p, s)?;
                    u3_rho_t(rho, t)
                }
                Region::R4 => {
                    let x = x4_ps(p, s)?;
                    let ts = tsat_k_from_p_mpa(p)?;
                    let (uv, ul) = if p < 16.529 {
                        (u2_p_t(p, ts)?, u1_p_t(p, ts)?)
                    } else {
                        let rho_v = 1.0 / v3_p_h(p, h4v_p_mpa(p)?)?;
                        let rho_l = 1.0 / v3_p_h(p, h4l_p_mpa(p)?)?;
                        (u3_rho_t(rho_v, ts)?, u3_rho_t(rho_l, ts)?)
                    };
                    Ok(x * uv + (1.0 - x) * ul)
                }
                Region::R5 => u5_p_t(p, t5_p_s(p, s)?),
            }
        }
        "cpv_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                cp2_p_t(p, ts)
            } else {
                let rho_v = 1.0 / v3_p_h(p, h4v_p_mpa(p)?)?;
                cp3_rho_t(rho_v, ts)
            }
        }
        "cpl_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                cp1_p_t(p, ts)
            } else {
                let rho_l = 1.0 / v3_p_h(p, h4l_p_mpa(p)?)?;
                cp3_rho_t(rho_l, ts)
            }
        }
        "cpv_t" => props(
            "CpV_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "cpl_t" => props(
            "CpL_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "cp_pt" => {
            let p = to_si_p_mpa_from_bar(in1);
            let t = to_si_t_k_from_c(in2);
            match region_from_p_t(p, t)? {
                Region::R1 => cp1_p_t(p, t),
                Region::R2 => cp2_p_t(p, t),
                Region::R3 => cp3_p_t(p, t),
                Region::R5 => cp5_p_t(p, t),
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
            }
        }
        "cp_ph" => {
            let p = to_si_p_mpa_from_bar(in1);
            let h = in2;
            match region_from_p_h(p, h)? {
                Region::R1 => cp1_p_t(p, t1_p_h(p, h)?),
                Region::R2 => cp2_p_t(p, t2_p_h(p, h)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_h(p, h)?;
                    let t = t3_p_h(p, h)?;
                    cp3_rho_t(rho, t)
                }
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
                Region::R5 => cp5_p_t(p, t5_p_h(p, h)?),
            }
        }
        "cp_ps" => {
            let p = to_si_p_mpa_from_bar(in1);
            let s = in2;
            match region_from_p_s(p, s)? {
                Region::R1 => cp1_p_t(p, t1_p_s(p, s)?),
                Region::R2 => cp2_p_t(p, t2_p_s(p, s)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_s(p, s)?;
                    let t = t3_p_s(p, s)?;
                    cp3_rho_t(rho, t)
                }
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
                Region::R5 => cp5_p_t(p, t5_p_s(p, s)?),
            }
        }
        "cvv_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                cv2_p_t(p, ts)
            } else {
                let rho_v = 1.0 / v3_p_h(p, h4v_p_mpa(p)?)?;
                cv3_rho_t(rho_v, ts)
            }
        }
        "cvl_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                cv1_p_t(p, ts)
            } else {
                let rho_l = 1.0 / v3_p_h(p, h4l_p_mpa(p)?)?;
                cv3_rho_t(rho_l, ts)
            }
        }
        "cvv_t" => props(
            "CvV_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "cvl_t" => props(
            "CvL_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "cv_pt" => {
            let p = to_si_p_mpa_from_bar(in1);
            let t = to_si_t_k_from_c(in2);
            match region_from_p_t(p, t)? {
                Region::R1 => cv1_p_t(p, t),
                Region::R2 => cv2_p_t(p, t),
                Region::R3 => cv3_p_t(p, t),
                Region::R5 => cv5_p_t(p, t),
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
            }
        }
        "cv_ph" => {
            let p = to_si_p_mpa_from_bar(in1);
            let h = in2;
            match region_from_p_h(p, h)? {
                Region::R1 => cv1_p_t(p, t1_p_h(p, h)?),
                Region::R2 => cv2_p_t(p, t2_p_h(p, h)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_h(p, h)?;
                    let t = t3_p_h(p, h)?;
                    cv3_rho_t(rho, t)
                }
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
                Region::R5 => cv5_p_t(p, t5_p_h(p, h)?),
            }
        }
        "cv_ps" => {
            let p = to_si_p_mpa_from_bar(in1);
            let s = in2;
            match region_from_p_s(p, s)? {
                Region::R1 => cv1_p_t(p, t1_p_s(p, s)?),
                Region::R2 => cv2_p_t(p, t2_p_s(p, s)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_s(p, s)?;
                    let t = t3_p_s(p, s)?;
                    cv3_rho_t(rho, t)
                }
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
                Region::R5 => cv5_p_t(p, t5_p_s(p, s)?),
            }
        }
        "wvv_p" | "wv_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                w2_p_t(p, ts)
            } else {
                let rho_v = 1.0 / v3_p_h(p, h4v_p_mpa(p)?)?;
                w3_rho_t(rho_v, ts)
            }
        }
        "wvl_p" | "wl_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let ts = tsat_k_from_p_mpa(p)?;
            if p < 16.529 {
                w1_p_t(p, ts)
            } else {
                let rho_l = 1.0 / v3_p_h(p, h4l_p_mpa(p)?)?;
                w3_rho_t(rho_l, ts)
            }
        }
        "wv_t" => props(
            "wV_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "wl_t" => props(
            "wL_p",
            from_si_p_bar_from_mpa(psat_mpa_from_t_k(to_si_t_k_from_c(in1))?),
            0.0,
        ),
        "w_pt" => {
            let p = to_si_p_mpa_from_bar(in1);
            let t = to_si_t_k_from_c(in2);
            match region_from_p_t(p, t)? {
                Region::R1 => w1_p_t(p, t),
                Region::R2 => w2_p_t(p, t),
                Region::R3 => w3_p_t(p, t),
                Region::R5 => w5_p_t(p, t),
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
            }
        }
        "w_ph" => {
            let p = to_si_p_mpa_from_bar(in1);
            let h = in2;
            match region_from_p_h(p, h)? {
                Region::R1 => w1_p_t(p, t1_p_h(p, h)?),
                Region::R2 => w2_p_t(p, t2_p_h(p, h)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_h(p, h)?;
                    let t = t3_p_h(p, h)?;
                    w3_rho_t(rho, t)
                }
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
                Region::R5 => w5_p_t(p, t5_p_h(p, h)?),
            }
        }
        "w_ps" => {
            let p = to_si_p_mpa_from_bar(in1);
            let s = in2;
            match region_from_p_s(p, s)? {
                Region::R1 => w1_p_t(p, t1_p_s(p, s)?),
                Region::R2 => w2_p_t(p, t2_p_s(p, s)?),
                Region::R3 => {
                    let rho = 1.0 / v3_p_s(p, s)?;
                    let t = t3_p_s(p, s)?;
                    w3_rho_t(rho, t)
                }
                Region::R4 => Err(If97Error::InvalidIntermediateValue),
                Region::R5 => w5_p_t(p, t5_p_s(p, s)?),
            }
        }
        "my_pt" => my_allregions_p_t(to_si_p_mpa_from_bar(in1), to_si_t_k_from_c(in2)),
        "my_ph" => my_allregions_p_h(to_si_p_mpa_from_bar(in1), in2),
        "my_ps" => {
            let h = props("h_ps", in1, in2)?;
            my_allregions_p_h(to_si_p_mpa_from_bar(in1), h)
        }
        "tcl_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let t = tsat_k_from_p_mpa(p)?;
            let rho = 1.0 / v4l_p_mpa(p)?;
            tc_ptrho(p, t, rho)
        }
        "tcv_p" => {
            let p = to_si_p_mpa_from_bar(in1);
            let t = tsat_k_from_p_mpa(p)?;
            let rho = 1.0 / v4v_p_mpa(p)?;
            tc_ptrho(p, t, rho)
        }
        "tcl_t" => {
            let t = to_si_t_k_from_c(in1);
            let p = psat_mpa_from_t_k(t)?;
            let rho = 1.0 / v4l_p_mpa(p)?;
            tc_ptrho(p, t, rho)
        }
        "tcv_t" => {
            let t = to_si_t_k_from_c(in1);
            let p = psat_mpa_from_t_k(t)?;
            let rho = 1.0 / v4v_p_mpa(p)?;
            tc_ptrho(p, t, rho)
        }
        "tc_pt" => {
            let v = props("v_pT", in1, in2)?;
            let p = to_si_p_mpa_from_bar(in1);
            let t = to_si_t_k_from_c(in2);
            tc_ptrho(p, t, 1.0 / v)
        }
        "tc_ph" => {
            let v = props("v_ph", in1, in2)?;
            let t_c = props("T_ph", in1, in2)?;
            let p = to_si_p_mpa_from_bar(in1);
            let t = to_si_t_k_from_c(t_c);
            tc_ptrho(p, t, 1.0 / v)
        }
        "tc_hs" => {
            let p_bar = props("p_hs", in1, in2)?;
            let v = props("v_ph", p_bar, in1)?;
            let t_c = props("T_ph", p_bar, in1)?;
            tc_ptrho(to_si_p_mpa_from_bar(p_bar), to_si_t_k_from_c(t_c), 1.0 / v)
        }
        "st_t" => surface_tension_t(to_si_t_k_from_c(in1)),
        "st_p" => {
            let t_c = props("Tsat_p", in1, 0.0)?;
            surface_tension_t(to_si_t_k_from_c(t_c))
        }
        "x_ph" => x4_ph(to_si_p_mpa_from_bar(in1), in2),
        "x_ps" => x4_ps(to_si_p_mpa_from_bar(in1), in2),
        "vx_ph" => {
            let p = to_si_p_mpa_from_bar(in1);
            let h = in2;
            let v_l = v4l_p_mpa(p)?;
            let v_v = v4v_p_mpa(p)?;
            let x = x4_ph(p, h)?;
            Ok(x * v_v / (x * v_v + (1.0 - x) * v_l))
        }
        "vx_ps" => {
            let p = to_si_p_mpa_from_bar(in1);
            let s = in2;
            let v_l = v4l_p_mpa(p)?;
            let v_v = v4v_p_mpa(p)?;
            let x = x4_ps(p, s)?;
            Ok(x * v_v / (x * v_v + (1.0 - x) * v_l))
        }
        _ => Err(If97Error::InvalidIntermediateValue),
    }
}

pub fn props_si(fun: &str, in1: f64, in2: f64) -> If97Result<f64> {
    let f = fun.trim().to_ascii_lowercase();
    match f.as_str() {
        "tsat_p" => Ok(to_si_t_k_from_c(props(
            "tsat_p",
            from_si_p_bar_from_pa(in1),
            0.0,
        )?)),
        "t_ph" => Ok(to_si_t_k_from_c(props(
            "t_ph",
            from_si_p_bar_from_pa(in1),
            in2,
        )?)),
        "t_ps" => Ok(to_si_t_k_from_c(props(
            "t_ps",
            from_si_p_bar_from_pa(in1),
            in2,
        )?)),
        "t_hs" => Ok(to_si_t_k_from_c(props("t_hs", in1, in2)?)),
        "psat_t" => Ok(to_si_p_pa_from_bar(props(
            "psat_t",
            from_si_t_c_from_k(in1),
            0.0,
        )?)),
        "p_hs" => Ok(to_si_p_pa_from_bar(props("p_hs", in1, in2)?)),
        "hv_p" => props("hv_p", from_si_p_bar_from_pa(in1), 0.0),
        "hl_p" => props("hl_p", from_si_p_bar_from_pa(in1), 0.0),
        "hv_t" => props("hv_t", from_si_t_c_from_k(in1), 0.0),
        "hl_t" => props("hl_t", from_si_t_c_from_k(in1), 0.0),
        "h_pt" => props("h_pt", from_si_p_bar_from_pa(in1), from_si_t_c_from_k(in2)),
        "h_ps" => props("h_ps", from_si_p_bar_from_pa(in1), in2),
        "h_px" => props("h_px", from_si_p_bar_from_pa(in1), in2),
        "h_prho" => props("h_prho", from_si_p_bar_from_pa(in1), in2),
        "h_tx" => props("h_tx", from_si_t_c_from_k(in1), in2),
        "vv_p" => props("vv_p", from_si_p_bar_from_pa(in1), 0.0),
        "vl_p" => props("vl_p", from_si_p_bar_from_pa(in1), 0.0),
        "vv_t" => props("vv_t", from_si_t_c_from_k(in1), 0.0),
        "vl_t" => props("vl_t", from_si_t_c_from_k(in1), 0.0),
        "v_pt" => props("v_pt", from_si_p_bar_from_pa(in1), from_si_t_c_from_k(in2)),
        "v_ph" => props("v_ph", from_si_p_bar_from_pa(in1), in2),
        "v_ps" => props("v_ps", from_si_p_bar_from_pa(in1), in2),
        "rhov_p" => props("rhov_p", from_si_p_bar_from_pa(in1), 0.0),
        "rhol_p" => props("rhol_p", from_si_p_bar_from_pa(in1), 0.0),
        "rhov_t" => props("rhov_t", from_si_t_c_from_k(in1), 0.0),
        "rhol_t" => props("rhol_t", from_si_t_c_from_k(in1), 0.0),
        "rho_pt" => props(
            "rho_pt",
            from_si_p_bar_from_pa(in1),
            from_si_t_c_from_k(in2),
        ),
        "rho_ph" => props("rho_ph", from_si_p_bar_from_pa(in1), in2),
        "rho_ps" => props("rho_ps", from_si_p_bar_from_pa(in1), in2),
        "sv_p" => props("sv_p", from_si_p_bar_from_pa(in1), 0.0),
        "sl_p" => props("sl_p", from_si_p_bar_from_pa(in1), 0.0),
        "sv_t" => props("sv_t", from_si_t_c_from_k(in1), 0.0),
        "sl_t" => props("sl_t", from_si_t_c_from_k(in1), 0.0),
        "s_pt" => props("s_pt", from_si_p_bar_from_pa(in1), from_si_t_c_from_k(in2)),
        "s_ph" => props("s_ph", from_si_p_bar_from_pa(in1), in2),
        "uv_p" => props("uv_p", from_si_p_bar_from_pa(in1), 0.0),
        "ul_p" => props("ul_p", from_si_p_bar_from_pa(in1), 0.0),
        "uv_t" => props("uv_t", from_si_t_c_from_k(in1), 0.0),
        "ul_t" => props("ul_t", from_si_t_c_from_k(in1), 0.0),
        "u_pt" => props("u_pt", from_si_p_bar_from_pa(in1), from_si_t_c_from_k(in2)),
        "u_ph" => props("u_ph", from_si_p_bar_from_pa(in1), in2),
        "u_ps" => props("u_ps", from_si_p_bar_from_pa(in1), in2),
        "cpv_p" => props("cpv_p", from_si_p_bar_from_pa(in1), 0.0),
        "cpl_p" => props("cpl_p", from_si_p_bar_from_pa(in1), 0.0),
        "cpv_t" => props("cpv_t", from_si_t_c_from_k(in1), 0.0),
        "cpl_t" => props("cpl_t", from_si_t_c_from_k(in1), 0.0),
        "cp_pt" => props("cp_pt", from_si_p_bar_from_pa(in1), from_si_t_c_from_k(in2)),
        "cp_ph" => props("cp_ph", from_si_p_bar_from_pa(in1), in2),
        "cp_ps" => props("cp_ps", from_si_p_bar_from_pa(in1), in2),
        "cvv_p" => props("cvv_p", from_si_p_bar_from_pa(in1), 0.0),
        "cvl_p" => props("cvl_p", from_si_p_bar_from_pa(in1), 0.0),
        "cvv_t" => props("cvv_t", from_si_t_c_from_k(in1), 0.0),
        "cvl_t" => props("cvl_t", from_si_t_c_from_k(in1), 0.0),
        "cv_pt" => props("cv_pt", from_si_p_bar_from_pa(in1), from_si_t_c_from_k(in2)),
        "cv_ph" => props("cv_ph", from_si_p_bar_from_pa(in1), in2),
        "cv_ps" => props("cv_ps", from_si_p_bar_from_pa(in1), in2),
        "wv_p" | "wvv_p" => props("wv_p", from_si_p_bar_from_pa(in1), 0.0),
        "wl_p" | "wvl_p" => props("wl_p", from_si_p_bar_from_pa(in1), 0.0),
        "wv_t" => props("wv_t", from_si_t_c_from_k(in1), 0.0),
        "wl_t" => props("wl_t", from_si_t_c_from_k(in1), 0.0),
        "w_pt" => props("w_pt", from_si_p_bar_from_pa(in1), from_si_t_c_from_k(in2)),
        "w_ph" => props("w_ph", from_si_p_bar_from_pa(in1), in2),
        "w_ps" => props("w_ps", from_si_p_bar_from_pa(in1), in2),
        "my_pt" => props("my_pt", from_si_p_bar_from_pa(in1), from_si_t_c_from_k(in2)),
        "my_ph" => props("my_ph", from_si_p_bar_from_pa(in1), in2),
        "my_ps" => props("my_ps", from_si_p_bar_from_pa(in1), in2),
        "tcl_p" => props("tcl_p", from_si_p_bar_from_pa(in1), 0.0),
        "tcv_p" => props("tcv_p", from_si_p_bar_from_pa(in1), 0.0),
        "tcl_t" => props("tcl_t", from_si_t_c_from_k(in1), 0.0),
        "tcv_t" => props("tcv_t", from_si_t_c_from_k(in1), 0.0),
        "tc_pt" => props("tc_pt", from_si_p_bar_from_pa(in1), from_si_t_c_from_k(in2)),
        "tc_ph" => props("tc_ph", from_si_p_bar_from_pa(in1), in2),
        "tc_hs" => props("tc_hs", in1, in2),
        "st_t" => props("st_t", from_si_t_c_from_k(in1), 0.0),
        "st_p" => props("st_p", from_si_p_bar_from_pa(in1), 0.0),
        "x_ph" => props("x_ph", from_si_p_bar_from_pa(in1), in2),
        "x_ps" => props("x_ps", from_si_p_bar_from_pa(in1), in2),
        "vx_ph" => props("vx_ph", from_si_p_bar_from_pa(in1), in2),
        "vx_ps" => props("vx_ps", from_si_p_bar_from_pa(in1), in2),
        _ => Err(If97Error::InvalidIntermediateValue),
    }
}

const R1_I: [i32; 34] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29,
    30, 31, 32,
];
const R1_J: [i32; 34] = [
    -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11,
    -6, -29, -31, -38, -39, -40, -41,
];
const R1_N: [f64; 34] = [
    0.146_329_712_131_67,
    -0.845_481_871_691_14,
    -3.756_360_367_204,
    3.385_516_916_838_5,
    -0.957_919_633_878_72,
    0.157_720_385_132_28,
    -0.016_616_417_199_501,
    8.121_462_998_356_8e-4,
    2.831_908_012_380_4e-4,
    -6.070_630_156_587_4e-4,
    -0.018_990_068_218_419,
    -0.032_529_748_770_505,
    -0.021_841_717_175_414,
    -5.283_835_796_993e-5,
    -4.718_432_107_326_7e-4,
    -3.000_178_079_302_6e-4,
    4.766_139_390_698_7e-5,
    -4.414_184_533_084_6e-6,
    -7.269_499_629_759_4e-16,
    -3.167_964_484_505_4e-5,
    -2.827_079_798_531_2e-6,
    -8.520_512_812_010_3e-10,
    -2.242_528_190_8e-6,
    -6.517_122_289_560_1e-7,
    -1.434_172_993_792_4e-13,
    -4.051_699_686_011_7e-7,
    -1.273_430_174_164_1e-9,
    -1.742_487_123_063_4e-10,
    -6.876_213_129_553_1e-19,
    1.447_830_782_852_1e-20,
    2.633_578_166_279_5e-23,
    -1.194_762_264_007_1e-23,
    1.822_809_458_140_4e-24,
    -9.353_708_729_245_8e-26,
];

fn r1_common(p_mpa: f64, t_k: f64) -> If97Result<(f64, f64)> {
    if !p_mpa.is_finite() || !t_k.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if p_mpa <= 0.0 || t_k <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let pi = p_mpa / 16.53;
    let tau = 1386.0 / t_k;
    Ok((pi, tau))
}

pub fn v1_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r1_common(p_mpa, t_k)?;
    let mut gamma_der_pi = 0.0;
    for i in 0..34 {
        gamma_der_pi -=
            R1_N[i] * (R1_I[i] as f64) * (7.1 - pi).powi(R1_I[i] - 1) * (tau - 1.222).powi(R1_J[i]);
    }
    let v = R_KJ_KG_K * t_k / p_mpa * pi * gamma_der_pi / 1000.0;
    if !v.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(v)
}

pub fn h1_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r1_common(p_mpa, t_k)?;
    let mut gamma_der_tau = 0.0;
    for i in 0..34 {
        gamma_der_tau +=
            R1_N[i] * (7.1 - pi).powi(R1_I[i]) * (R1_J[i] as f64) * (tau - 1.222).powi(R1_J[i] - 1);
    }
    let h = R_KJ_KG_K * t_k * tau * gamma_der_tau;
    if !h.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(h)
}

pub fn u1_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r1_common(p_mpa, t_k)?;
    let mut gamma_der_tau = 0.0;
    let mut gamma_der_pi = 0.0;
    for i in 0..34 {
        gamma_der_pi -=
            R1_N[i] * (R1_I[i] as f64) * (7.1 - pi).powi(R1_I[i] - 1) * (tau - 1.222).powi(R1_J[i]);
        gamma_der_tau +=
            R1_N[i] * (7.1 - pi).powi(R1_I[i]) * (R1_J[i] as f64) * (tau - 1.222).powi(R1_J[i] - 1);
    }
    let u = R_KJ_KG_K * t_k * (tau * gamma_der_tau - pi * gamma_der_pi);
    if !u.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(u)
}

pub fn s1_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r1_common(p_mpa, t_k)?;
    let mut gamma = 0.0;
    let mut gamma_der_tau = 0.0;
    for i in 0..34 {
        gamma_der_tau +=
            R1_N[i] * (7.1 - pi).powi(R1_I[i]) * (R1_J[i] as f64) * (tau - 1.222).powi(R1_J[i] - 1);
        gamma += R1_N[i] * (7.1 - pi).powi(R1_I[i]) * (tau - 1.222).powi(R1_J[i]);
    }
    let s = R_KJ_KG_K * tau * gamma_der_tau - R_KJ_KG_K * gamma;
    if !s.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(s)
}

pub fn cp1_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r1_common(p_mpa, t_k)?;
    let mut gamma_der_tautau = 0.0;
    for i in 0..34 {
        gamma_der_tautau += R1_N[i]
            * (7.1 - pi).powi(R1_I[i])
            * (R1_J[i] as f64)
            * ((R1_J[i] - 1) as f64)
            * (tau - 1.222).powi(R1_J[i] - 2);
    }
    let cp = -R_KJ_KG_K * tau * tau * gamma_der_tautau;
    if !cp.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(cp)
}

pub fn cv1_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r1_common(p_mpa, t_k)?;
    let mut gamma_der_pi = 0.0;
    let mut gamma_der_pipi = 0.0;
    let mut gamma_der_pitau = 0.0;
    let mut gamma_der_tautau = 0.0;
    for i in 0..34 {
        gamma_der_pi -=
            R1_N[i] * (R1_I[i] as f64) * (7.1 - pi).powi(R1_I[i] - 1) * (tau - 1.222).powi(R1_J[i]);
        gamma_der_pipi += R1_N[i]
            * (R1_I[i] as f64)
            * ((R1_I[i] - 1) as f64)
            * (7.1 - pi).powi(R1_I[i] - 2)
            * (tau - 1.222).powi(R1_J[i]);
        gamma_der_pitau -= R1_N[i]
            * (R1_I[i] as f64)
            * (7.1 - pi).powi(R1_I[i] - 1)
            * (R1_J[i] as f64)
            * (tau - 1.222).powi(R1_J[i] - 1);
        gamma_der_tautau += R1_N[i]
            * (7.1 - pi).powi(R1_I[i])
            * (R1_J[i] as f64)
            * ((R1_J[i] - 1) as f64)
            * (tau - 1.222).powi(R1_J[i] - 2);
    }
    if gamma_der_pipi == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let cv = R_KJ_KG_K
        * (-tau * tau * gamma_der_tautau
            + (gamma_der_pi - tau * gamma_der_pitau).powi(2) / gamma_der_pipi);
    if !cv.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(cv)
}

pub fn w1_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r1_common(p_mpa, t_k)?;
    let mut gamma_der_pi = 0.0;
    let mut gamma_der_pipi = 0.0;
    let mut gamma_der_pitau = 0.0;
    let mut gamma_der_tautau = 0.0;
    for i in 0..34 {
        gamma_der_pi -=
            R1_N[i] * (R1_I[i] as f64) * (7.1 - pi).powi(R1_I[i] - 1) * (tau - 1.222).powi(R1_J[i]);
        gamma_der_pipi += R1_N[i]
            * (R1_I[i] as f64)
            * ((R1_I[i] - 1) as f64)
            * (7.1 - pi).powi(R1_I[i] - 2)
            * (tau - 1.222).powi(R1_J[i]);
        gamma_der_pitau -= R1_N[i]
            * (R1_I[i] as f64)
            * (7.1 - pi).powi(R1_I[i] - 1)
            * (R1_J[i] as f64)
            * (tau - 1.222).powi(R1_J[i] - 1);
        gamma_der_tautau += R1_N[i]
            * (7.1 - pi).powi(R1_I[i])
            * (R1_J[i] as f64)
            * ((R1_J[i] - 1) as f64)
            * (tau - 1.222).powi(R1_J[i] - 2);
    }
    let denom = (gamma_der_pi - tau * gamma_der_pitau).powi(2) / (tau * tau * gamma_der_tautau)
        - gamma_der_pipi;
    if denom == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let w2 = 1000.0 * R_KJ_KG_K * t_k * gamma_der_pi.powi(2) / denom;
    let w = sqrt_nonneg(w2)?;
    if !w.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(w)
}

const T1_PH_I: [i32; 20] = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6];
const T1_PH_J: [i32; 20] = [
    0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32,
];
const T1_PH_N: [f64; 20] = [
    -238.724_899_245_21,
    404.211_886_379_45,
    113.497_468_817_18,
    -5.845_761_604_803_9,
    -1.528_548_241_314e-4,
    -1.086_670_769_537_7e-6,
    -13.391_744_872_602,
    43.211_039_183_559,
    -54.010_067_170_506,
    30.535_892_203_916,
    -6.596_474_942_363_8,
    9.396_540_087_836_3e-3,
    1.157_364_750_534e-7,
    -2.585_864_128_207_3e-5,
    -4.064_436_308_479_9e-9,
    6.645_618_619_163_5e-8,
    8.067_073_410_302_7e-11,
    -9.347_777_121_394_7e-13,
    5.826_544_202_060_1e-15,
    -1.502_018_595_350_3e-17,
];

pub fn t1_p_h(p_mpa: f64, h: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !h.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let pi = p_mpa;
    let eta = h / 2500.0;
    let mut t = 0.0;
    for i in 0..20 {
        t += T1_PH_N[i] * pi.powi(T1_PH_I[i]) * (eta + 1.0).powi(T1_PH_J[i]);
    }
    if !t.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(t)
}

const T1_PS_I: [i32; 20] = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4];
const T1_PS_J: [i32; 20] = [
    0, 1, 2, 3, 11, 31, 0, 1, 2, 3, 12, 31, 0, 1, 2, 9, 31, 10, 32, 32,
];
const T1_PS_N: [f64; 20] = [
    174.782_680_583_07,
    34.806_930_892_873,
    6.529_258_497_845_5,
    0.330_399_817_754_89,
    -1.928_138_292_319_6e-7,
    -2.490_919_724_457_3e-23,
    -0.261_076_364_893_32,
    0.225_929_659_815_86,
    -0.064_256_463_395_226,
    7.887_628_927_052_6e-3,
    3.567_211_060_736_6e-10,
    1.733_249_699_489_5e-24,
    5.660_890_065_483_7e-4,
    -3.263_548_313_971_7e-4,
    4.477_828_669_063_2e-5,
    -5.132_215_690_850_7e-10,
    -4.252_265_704_220_7e-26,
    2.640_044_136_068_9e-13,
    7.812_460_045_972_3e-29,
    -3.073_219_990_366_8e-31,
];

pub fn t1_p_s(p_mpa: f64, s: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let pi = p_mpa;
    let sigma = s;
    let mut t = 0.0;
    for i in 0..20 {
        t += T1_PS_N[i] * pi.powi(T1_PS_I[i]) * (sigma + 2.0).powi(T1_PS_J[i]);
    }
    if !t.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(t)
}

const P1_HS_I: [i32; 19] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 5];
const P1_HS_J: [i32; 19] = [0, 1, 2, 4, 5, 6, 8, 14, 0, 1, 4, 6, 0, 1, 10, 4, 1, 4, 0];
const P1_HS_N: [f64; 19] = [
    -0.691_997_014_660_582,
    -18.361_254_878_756,
    -9.283_324_092_973_35,
    65.963_956_990_990_6,
    -16.206_038_891_202_4,
    450.620_017_338_667,
    854.680_678_224_17,
    6_075.232_140_011_62,
    32.648_768_262_185_6,
    -26.940_884_458_293_1,
    -319.947_848_334_3,
    -928.354_307_043_32,
    30.363_453_745_524_9,
    -65.054_042_244_414_6,
    -4_309.913_165_161_3,
    -747.512_324_096_068,
    730.000_345_529_245,
    1_142.840_325_690_21,
    -436.407_041_874_559,
];

pub fn p1_h_s(h: f64, s: f64) -> If97Result<f64> {
    if !h.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let eta = h / 3400.0;
    let sigma = s / 7.6;
    let mut p = 0.0;
    for i in 0..19 {
        p += P1_HS_N[i] * (eta + 0.05).powi(P1_HS_I[i]) * (sigma + 0.05).powi(P1_HS_J[i]);
    }
    let p = p * 100.0;
    if !p.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(p)
}

pub fn t1_p_rho(p_mpa: f64, rho: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !rho.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if rho <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let mut low = 273.15;
    let mut high = tsat_k_from_p_mpa(p_mpa)?;
    for _ in 0..200 {
        let ts = (low + high) / 2.0;
        let rhos = 1.0 / v1_p_t(p_mpa, ts)?;
        if (rho - rhos).abs() <= 1e-5 {
            return Ok(ts);
        }
        if rhos < rho {
            high = ts;
        } else {
            low = ts;
        }
    }
    Err(If97Error::InvalidIntermediateValue)
}

const R2_J0: [i32; 9] = [0, 1, -5, -4, -3, -2, -1, 2, 3];
const R2_N0: [f64; 9] = [
    -9.692_768_650_021_7,
    10.086_655_968_018,
    -0.005_608_791_128_302,
    0.071_452_738_081_455,
    -0.407_104_982_239_28,
    1.424_081_917_144_4,
    -4.383_951_131_945,
    -0.284_086_324_607_72,
    0.021_268_463_753_307,
];
const R2_IR: [i32; 43] = [
    1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10,
    16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24,
];
const R2_JR: [i32; 43] = [
    0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4,
    10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58,
];
const R2_NR: [f64; 43] = [
    -1.773_174_247_321_3e-3,
    -0.017_834_862_292_358,
    -0.045_996_013_696_365,
    -0.057_581_259_083_432,
    -0.050_325_278_727_93,
    -3.303_264_167_020_3e-5,
    -1.894_898_751_631_5e-4,
    -3.939_277_724_335_5e-3,
    -0.043_797_295_650_573,
    -2.667_454_791_408_7e-5,
    2.048_173_769_230_9e-8,
    4.387_066_728_443_5e-7,
    -3.227_767_723_857e-5,
    -1.503_392_454_214_8e-3,
    -0.040_668_253_562_649,
    -7.884_730_955_936_7e-10,
    1.279_071_785_228_5e-8,
    4.822_537_271_850_7e-7,
    2.292_207_633_766_1e-6,
    -1.671_476_645_106_1e-11,
    -2.117_147_232_135_5e-3,
    -23.895_741_934_104,
    -5.905_956_432_427e-18,
    -1.262_180_889_910_1e-6,
    -0.038_946_842_435_739,
    1.125_621_136_045_9e-11,
    -8.231_134_089_799_8,
    1.980_971_280_208_8e-8,
    1.040_696_521_017_4e-19,
    -1.023_474_709_592_9e-13,
    -1.001_817_937_951_1e-9,
    -8.088_290_864_698_5e-11,
    0.106_930_318_794_09,
    -0.336_622_505_741_71,
    8.918_584_535_542_1e-25,
    3.062_931_687_623_2e-13,
    -4.200_246_769_820_8e-6,
    -5.905_602_968_563_9e-26,
    3.782_694_761_345_7e-6,
    -1.276_860_893_468_1e-15,
    7.308_761_059_506_1e-29,
    5.541_471_535_077_8e-17,
    -9.436_970_724_121e-7,
];

fn r2_common(p_mpa: f64, t_k: f64) -> If97Result<(f64, f64)> {
    if !p_mpa.is_finite() || !t_k.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if p_mpa <= 0.0 || t_k <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let pi = p_mpa;
    let tau = 540.0 / t_k;
    Ok((pi, tau))
}

pub fn v2_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r2_common(p_mpa, t_k)?;
    let g0_pi = 1.0 / pi;
    let mut gr_pi = 0.0;
    for i in 0..43 {
        gr_pi += R2_NR[i] * (R2_IR[i] as f64) * pi.powi(R2_IR[i] - 1) * (tau - 0.5).powi(R2_JR[i]);
    }
    let v = R_KJ_KG_K * t_k / p_mpa * pi * (g0_pi + gr_pi) / 1000.0;
    if !v.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(v)
}

pub fn h2_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r2_common(p_mpa, t_k)?;
    let mut g0_tau = 0.0;
    for i in 0..9 {
        g0_tau += R2_N0[i] * (R2_J0[i] as f64) * tau.powi(R2_J0[i] - 1);
    }
    let mut gr_tau = 0.0;
    for i in 0..43 {
        gr_tau += R2_NR[i] * pi.powi(R2_IR[i]) * (R2_JR[i] as f64) * (tau - 0.5).powi(R2_JR[i] - 1);
    }
    let h = R_KJ_KG_K * t_k * tau * (g0_tau + gr_tau);
    if !h.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(h)
}

pub fn u2_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r2_common(p_mpa, t_k)?;
    let g0_pi = 1.0 / pi;
    let mut g0_tau = 0.0;
    for i in 0..9 {
        g0_tau += R2_N0[i] * (R2_J0[i] as f64) * tau.powi(R2_J0[i] - 1);
    }
    let mut gr_pi = 0.0;
    let mut gr_tau = 0.0;
    for i in 0..43 {
        gr_pi += R2_NR[i] * (R2_IR[i] as f64) * pi.powi(R2_IR[i] - 1) * (tau - 0.5).powi(R2_JR[i]);
        gr_tau += R2_NR[i] * pi.powi(R2_IR[i]) * (R2_JR[i] as f64) * (tau - 0.5).powi(R2_JR[i] - 1);
    }
    let u = R_KJ_KG_K * t_k * (tau * (g0_tau + gr_tau) - pi * (g0_pi + gr_pi));
    if !u.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(u)
}

pub fn s2_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r2_common(p_mpa, t_k)?;
    let mut g0 = pi.ln();
    let mut g0_tau = 0.0;
    for i in 0..9 {
        g0 += R2_N0[i] * tau.powi(R2_J0[i]);
        g0_tau += R2_N0[i] * (R2_J0[i] as f64) * tau.powi(R2_J0[i] - 1);
    }
    let mut gr = 0.0;
    let mut gr_tau = 0.0;
    for i in 0..43 {
        gr += R2_NR[i] * pi.powi(R2_IR[i]) * (tau - 0.5).powi(R2_JR[i]);
        gr_tau += R2_NR[i] * pi.powi(R2_IR[i]) * (R2_JR[i] as f64) * (tau - 0.5).powi(R2_JR[i] - 1);
    }
    let s = R_KJ_KG_K * (tau * (g0_tau + gr_tau) - (g0 + gr));
    if !s.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(s)
}

pub fn cp2_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r2_common(p_mpa, t_k)?;
    let mut g0_tautau = 0.0;
    for i in 0..9 {
        g0_tautau +=
            R2_N0[i] * (R2_J0[i] as f64) * ((R2_J0[i] - 1) as f64) * tau.powi(R2_J0[i] - 2);
    }
    let mut gr_tautau = 0.0;
    for i in 0..43 {
        gr_tautau += R2_NR[i]
            * pi.powi(R2_IR[i])
            * (R2_JR[i] as f64)
            * ((R2_JR[i] - 1) as f64)
            * (tau - 0.5).powi(R2_JR[i] - 2);
    }
    let cp = -R_KJ_KG_K * tau * tau * (g0_tautau + gr_tautau);
    if !cp.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(cp)
}

pub fn cv2_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r2_common(p_mpa, t_k)?;
    let mut g0_tautau = 0.0;
    for i in 0..9 {
        g0_tautau +=
            R2_N0[i] * (R2_J0[i] as f64) * ((R2_J0[i] - 1) as f64) * tau.powi(R2_J0[i] - 2);
    }
    let mut gr_pi = 0.0;
    let mut gr_pitau = 0.0;
    let mut gr_pipi = 0.0;
    let mut gr_tautau = 0.0;
    for i in 0..43 {
        gr_pi += R2_NR[i] * (R2_IR[i] as f64) * pi.powi(R2_IR[i] - 1) * (tau - 0.5).powi(R2_JR[i]);
        gr_pipi += R2_NR[i]
            * (R2_IR[i] as f64)
            * ((R2_IR[i] - 1) as f64)
            * pi.powi(R2_IR[i] - 2)
            * (tau - 0.5).powi(R2_JR[i]);
        gr_pitau += R2_NR[i]
            * (R2_IR[i] as f64)
            * pi.powi(R2_IR[i] - 1)
            * (R2_JR[i] as f64)
            * (tau - 0.5).powi(R2_JR[i] - 1);
        gr_tautau += R2_NR[i]
            * pi.powi(R2_IR[i])
            * (R2_JR[i] as f64)
            * ((R2_JR[i] - 1) as f64)
            * (tau - 0.5).powi(R2_JR[i] - 2);
    }
    let denom = 1.0 - pi * pi * gr_pipi;
    if denom == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let cv = R_KJ_KG_K
        * (-tau * tau * (g0_tautau + gr_tautau)
            - (1.0 + pi * gr_pi - tau * pi * gr_pitau).powi(2) / denom);
    if !cv.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(cv)
}

pub fn w2_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let (pi, tau) = r2_common(p_mpa, t_k)?;
    let mut g0_tautau = 0.0;
    for i in 0..9 {
        g0_tautau +=
            R2_N0[i] * (R2_J0[i] as f64) * ((R2_J0[i] - 1) as f64) * tau.powi(R2_J0[i] - 2);
    }
    let mut gr_pi = 0.0;
    let mut gr_pitau = 0.0;
    let mut gr_pipi = 0.0;
    let mut gr_tautau = 0.0;
    for i in 0..43 {
        gr_pi += R2_NR[i] * (R2_IR[i] as f64) * pi.powi(R2_IR[i] - 1) * (tau - 0.5).powi(R2_JR[i]);
        gr_pipi += R2_NR[i]
            * (R2_IR[i] as f64)
            * ((R2_IR[i] - 1) as f64)
            * pi.powi(R2_IR[i] - 2)
            * (tau - 0.5).powi(R2_JR[i]);
        gr_pitau += R2_NR[i]
            * (R2_IR[i] as f64)
            * pi.powi(R2_IR[i] - 1)
            * (R2_JR[i] as f64)
            * (tau - 0.5).powi(R2_JR[i] - 1);
        gr_tautau += R2_NR[i]
            * pi.powi(R2_IR[i])
            * (R2_JR[i] as f64)
            * ((R2_JR[i] - 1) as f64)
            * (tau - 0.5).powi(R2_JR[i] - 2);
    }
    let term1 = 1.0 + 2.0 * pi * gr_pi + pi * pi * gr_pi * gr_pi;
    let denom = (1.0 - pi * pi * gr_pipi)
        + (1.0 + pi * gr_pi - tau * pi * gr_pitau).powi(2) / (tau * tau * (g0_tautau + gr_tautau));
    if denom == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let w2 = 1000.0 * R_KJ_KG_K * t_k * term1 / denom;
    let w = sqrt_nonneg(w2)?;
    if !w.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(w)
}

pub fn t2_p_h(p_mpa: f64, h: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !h.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let sub_reg = if p_mpa < 4.0 {
        1
    } else if p_mpa
        < (905.842_785_147_23 - 0.679_557_863_992_41 * h + 1.280_900_273_013_6e-4 * h * h)
    {
        2
    } else {
        3
    };

    match sub_reg {
        1 => {
            let ji: [i32; 34] = [
                0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40, 42, 44, 24,
                44, 12, 32, 44, 32, 36, 42, 34, 44, 28,
            ];
            let ii: [i32; 34] = [
                0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4,
                5, 5, 5, 6, 6, 7,
            ];
            let ni: [f64; 34] = [
                1_089.895_231_828_8,
                849.516_544_955_35,
                -107.817_480_918_26,
                33.153_654_801_263,
                -7.423_201_679_024_8,
                11.765_048_724_356,
                1.844_574_935_579,
                -4.179_270_054_962_4,
                6.247_819_693_581_2,
                -17.344_563_108_114,
                -200.581_768_620_96,
                271.960_654_737_96,
                -455.113_182_858_18,
                3_091.968_860_475_5,
                252_266.403_578_72,
                -6.170_742_286_833_9e-3,
                -0.310_780_466_295_83,
                11.670_873_077_107,
                128_127_984.040_46,
                -985_549_096.232_76,
                2_822_454_697.300_2,
                -3_594_897_141.070_3,
                1_722_734_991.319_7,
                -13_551.334_240_775,
                12_848_734.664_65,
                1.386_572_428_322_6,
                235_988.325_565_14,
                -13_105_236.545_054,
                7_399.983_547_476_6,
                -551_966.970_300_6,
                3_715_408.599_623_3,
                19_127.729_239_66,
                -415_351.648_356_34,
                -62.459_855_192_507,
            ];
            let hs = h / 2000.0;
            let mut t = 0.0;
            for i in 0..34 {
                t += ni[i] * p_mpa.powi(ii[i]) * (hs - 2.1).powi(ji[i]);
            }
            Ok(t)
        }
        2 => {
            let ji: [i32; 38] = [
                0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18, 40, 1, 2, 12,
                24, 2, 12, 18, 24, 28, 40, 18, 24, 40, 28, 2, 28, 1, 40,
            ];
            let ii: [i32; 38] = [
                0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4,
                4, 4, 5, 5, 5, 6, 7, 7, 9, 9,
            ];
            let ni: [f64; 38] = [
                1_489.504_107_951_6,
                743.077_983_140_34,
                -97.708_318_797_837,
                2.474_246_470_567_4,
                -0.632_813_200_160_26,
                1.138_595_212_965_8,
                -0.478_118_636_486_25,
                8.520_812_343_154_4e-3,
                0.937_471_473_779_32,
                3.359_311_860_491_6,
                3.380_935_560_145_4,
                0.168_445_396_719_04,
                0.738_757_452_366_95,
                -0.471_287_374_361_86,
                0.150_202_731_397_07,
                -2.176_411_421_975e-3,
                -0.021_810_755_324_761,
                -0.108_297_844_036_77,
                -0.046_333_324_635_812,
                7.128_035_195_955_1e-5,
                1.103_283_178_999_9e-4,
                1.895_524_838_790_2e-4,
                3.089_154_116_053_7e-3,
                1.355_550_455_494_9e-3,
                2.864_023_747_745_6e-7,
                -1.077_985_735_751_2e-5,
                -7.646_271_245_481_4e-5,
                1.405_239_281_831_6e-5,
                -3.108_381_433_143_4e-5,
                -1.030_273_821_210_3e-6,
                2.821_728_163_504e-7,
                1.270_490_227_194_5e-6,
                7.380_335_346_829_2e-8,
                -1.103_013_923_890_9e-8,
                -8.145_636_520_783_3e-14,
                -2.518_054_568_296_2e-11,
                -1.756_523_396_940_7e-18,
                8.693_415_634_416_3e-15,
            ];
            let hs = h / 2000.0;
            let mut t = 0.0;
            for i in 0..38 {
                t += ni[i] * (p_mpa - 2.0).powi(ii[i]) * (hs - 2.6).powi(ji[i]);
            }
            Ok(t)
        }
        _ => {
            let ji: [i32; 23] = [
                0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22,
            ];
            let ii: [i32; 23] = [
                -7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6,
            ];
            let ni: [f64; 23] = [
                -3_236_839_855_524.2,
                7_326_335_090_218.1,
                358_250_899_454.47,
                -583_401_318_515.9,
                -10_783_068_217.47,
                20_825_544_563.171,
                610_747.835_645_16,
                859_777.225_355_8,
                -25_745.723_604_17,
                31_081.088_422_714,
                1_208.231_586_593_6,
                482.197_551_092_55,
                3.796_600_127_248_6,
                -10.842_984_880_077,
                -0.045_364_172_676_66,
                1.455_911_565_869_8e-13,
                1.126_159_740_723e-12,
                -1.780_498_224_068_6e-11,
                1.232_457_969_083_2e-7,
                -1.160_692_113_098_4e-6,
                2.784_636_708_855_4e-5,
                -5.927_003_847_417_6e-4,
                1.291_858_299_187_8e-3,
            ];
            let hs = h / 2000.0;
            let mut t = 0.0;
            for i in 0..23 {
                t += ni[i] * (p_mpa + 25.0).powi(ii[i]) * (hs - 1.8).powi(ji[i]);
            }
            Ok(t)
        }
    }
}

#[allow(clippy::approx_constant)]
pub fn t2_p_s(p_mpa: f64, s: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let sub_reg = if p_mpa < 4.0 {
        1
    } else if s < 5.85 {
        3
    } else {
        2
    };

    match sub_reg {
        1 => {
            let ii: [f64; 46] = [
                -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.25, -1.25, -1.25, -1.0, -1.0, -1.0, -1.0,
                -1.0, -1.0, -0.75, -0.75, -0.5, -0.5, -0.5, -0.5, -0.25, -0.25, -0.25, -0.25, 0.25,
                0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 1.0,
                1.0, 1.25, 1.25, 1.5, 1.5,
            ];
            let ji: [i32; 46] = [
                -24, -23, -19, -13, -11, -10, -19, -15, -6, -26, -21, -17, -16, -9, -8, -15, -14,
                -26, -13, -9, -7, -27, -25, -11, -6, 1, 4, 8, 11, 0, 1, 5, 6, 10, 14, 16, 0, 4, 9,
                17, 7, 18, 3, 15, 5, 18,
            ];
            let ni: [f64; 46] = [
                -392_359.838_619_84,
                515_265.738_272_7,
                40_482.443_161_048,
                -321.937_909_239_02,
                96.961_424_218_694,
                -22.867_846_371_773,
                -449_429.141_243_57,
                -5_011.833_602_016_6,
                0.356_844_635_600_15,
                44_235.335_848_19,
                -13_673.388_811_708,
                421_632.602_078_64,
                22_516.925_837_475,
                474.421_448_656_46,
                -149.311_307_976_47,
                -197_811.263_204_52,
                -23_554.399_470_76,
                -19_070.616_302_076,
                55_375.669_883_164,
                3_829.369_143_736_3,
                -603.918_605_805_67,
                1_936.310_262_033_1,
                4_266.064_369_861,
                -5_978.063_887_271_8,
                -704.014_639_268_62,
                338.367_841_075_53,
                20.862_786_635_187,
                0.033_834_172_656_196,
                -4.312_442_841_489_3e-5,
                166.537_913_564_12,
                -139.862_920_558_98,
                -0.788_495_479_998_72,
                0.072_132_411_753_872,
                -5.975_483_939_828_3e-3,
                -1.214_135_895_390_4e-5,
                2.322_709_673_387_1e-7,
                -10.538_463_566_194,
                2.071_892_549_650_2,
                -0.072_193_155_260_427,
                2.074_988_708_112e-7,
                -0.018_340_657_911_379,
                2.903_627_234_869_6e-7,
                0.210_375_278_936_19,
                2.568_123_972_999_9e-4,
                -0.012_799_002_933_781,
                -8.219_810_265_201_8e-6,
            ];
            let pi = p_mpa;
            let sigma = s / 2.0;
            let mut teta = 0.0;
            for i in 0..46 {
                teta += ni[i] * pi.powf(ii[i]) * (sigma - 2.0).powi(ji[i]);
            }
            Ok(teta)
        }
        2 => {
            let ii: [i32; 44] = [
                -6, -6, -5, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, 0,
                0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5,
            ];
            let ji: [i32; 44] = [
                0, 11, 0, 11, 0, 1, 11, 0, 1, 11, 12, 0, 1, 6, 10, 0, 1, 5, 8, 9, 0, 1, 2, 4, 5, 6,
                9, 0, 1, 2, 3, 7, 8, 0, 1, 5, 0, 1, 3, 0, 1, 0, 1, 2,
            ];
            let ni: [f64; 44] = [
                316_876.650_834_97,
                20.864_175_881_858,
                -398_593.998_035_99,
                -21.816_058_518_877,
                223_697.851_942_42,
                -2_784.170_344_581_7,
                9.920_743_607_148,
                -75_197.512_299_157,
                2_970.860_595_115_8,
                -3.440_687_854_852_6,
                0.388_155_642_491_15,
                17_511.295_085_75,
                -1_423.711_285_444_9,
                1.094_380_336_416_7,
                0.899_716_193_084_95,
                -3_375.974_009_895_8,
                471.628_858_183_55,
                -1.918_824_199_367_9,
                0.410_785_804_921_96,
                -0.334_653_781_720_97,
                1_387.003_477_750_5,
                -406.633_261_958_38,
                41.727_347_159_61,
                2.193_254_943_453_2,
                -1.032_005_000_907_7,
                0.358_829_435_167_03,
                5.251_145_372_606_6e-3,
                12.838_916_450_705,
                -2.864_243_721_938_1,
                0.569_126_836_648_55,
                -0.099_962_954_584_931,
                -3.263_203_777_845_9e-3,
                2.332_092_257_672_3e-4,
                -0.153_348_098_574_5,
                0.029_072_288_239_902,
                3.753_470_274_116_7e-4,
                1.729_669_170_241_1e-3,
                -3.855_605_084_450_4e-4,
                -3.501_771_229_260_8e-5,
                -1.456_639_363_149_2e-5,
                5.642_085_726_726_9e-6,
                4.128_615_007_460_5e-8,
                -2.068_467_111_882_4e-8,
                1.640_939_367_472_5e-9,
            ];
            let pi = p_mpa;
            let sigma = s / 0.7853;
            let mut teta = 0.0;
            for i in 0..44 {
                teta += ni[i] * pi.powi(ii[i]) * (10.0 - sigma).powi(ji[i]);
            }
            Ok(teta)
        }
        _ => {
            let ii: [i32; 30] = [
                -2, -2, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7,
                7, 7, 7,
            ];
            let ji: [i32; 30] = [
                0, 1, 0, 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 2, 0, 1, 5, 0, 1, 4, 0, 1, 2, 0, 1, 0, 1, 3,
                4, 5,
            ];
            let ni: [f64; 30] = [
                909.685_010_053_65,
                2_404.566_708_842,
                -591.623_263_871_3,
                541.454_041_280_74,
                -270.983_084_111_92,
                979.765_250_979_26,
                -469.667_729_594_35,
                14.399_274_604_723,
                -19.104_204_230_429,
                5.329_916_711_197_1,
                -21.252_975_375_934,
                -0.311_473_344_137_6,
                0.603_348_408_946_23,
                -0.042_764_839_702_509,
                5.818_559_725_525_9e-3,
                -0.014_597_008_284_753,
                5.663_117_563_102_7e-3,
                -7.615_586_458_457_7e-5,
                2.244_034_291_933_2e-4,
                -1.256_109_501_341_3e-5,
                6.332_313_266_093_4e-7,
                -2.054_198_967_537_5e-6,
                3.640_537_039_008_2e-8,
                -2.975_989_778_921_5e-9,
                1.013_661_852_976_3e-8,
                5.992_571_969_235_1e-12,
                -2.067_787_010_516_4e-11,
                -2.087_427_818_188_6e-11,
                1.016_216_682_508_9e-10,
                -1.642_982_828_134_7e-10,
            ];
            let pi = p_mpa;
            let sigma = s / 2.9251;
            let mut teta = 0.0;
            for i in 0..30 {
                teta += ni[i] * pi.powi(ii[i]) * (2.0 - sigma).powi(ji[i]);
            }
            Ok(teta)
        }
    }
}

pub fn p2_h_s(h: f64, s: f64) -> If97Result<f64> {
    if !h.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let limit = -3_498.980_834_321_39 + 2_575.607_169_058_76 * s - 421.073_558_227_969 * s * s
        + 27.634_906_379_994_4 * s * s * s;
    let sub_reg = if h < limit {
        1
    } else if s < 5.85 {
        3
    } else {
        2
    };

    match sub_reg {
        1 => {
            let ii: [i32; 29] = [
                0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 6,
                7,
            ];
            let ji: [i32; 29] = [
                1, 3, 6, 16, 20, 22, 0, 1, 2, 3, 5, 6, 10, 16, 20, 22, 3, 16, 20, 0, 2, 3, 6, 16,
                16, 3, 16, 3, 1,
            ];
            let ni: [f64; 29] = [
                -1.825_753_619_230_32e-2,
                -0.125_229_548_799_536,
                0.592_290_437_320_145,
                6.047_697_061_851_22,
                238.624_965_444_474,
                -298.639_090_222_922,
                0.051_225_081_304_075,
                -0.437_266_515_606_486,
                0.413_336_902_999_504,
                -5.164_682_545_747_73,
                -5.570_148_384_457_11,
                12.855_503_782_447_8,
                11.414_410_895_329,
                -119.504_225_652_714,
                -2_847.779_859_615_6,
                4_317.578_464_080_06,
                1.128_940_408_026_5,
                1_974.091_862_063_19,
                1_516.124_447_060_87,
                1.413_244_514_212_35e-2,
                0.585_501_282_219_601,
                -2.972_580_758_630_12,
                5.945_673_148_473_19,
                -6_236.565_657_989_05,
                9_659.862_351_333_32,
                6.815_009_349_481_34,
                -6_332.072_868_244_89,
                -5.589_192_244_657_6,
                4.006_457_984_720_63e-2,
            ];
            let eta = h / 4200.0;
            let sigma = s / 12.0;
            let mut pi = 0.0;
            for i in 0..29 {
                pi += ni[i] * (eta - 0.5).powi(ii[i]) * (sigma - 1.2).powi(ji[i]);
            }
            Ok(pi.powi(4) * 4.0)
        }
        2 => {
            let ii: [i32; 33] = [
                0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8,
                8, 8, 8, 12, 14,
            ];
            let ji: [i32; 33] = [
                0, 1, 2, 4, 8, 0, 1, 2, 3, 5, 12, 1, 6, 18, 0, 1, 7, 12, 1, 16, 1, 12, 1, 8, 18, 1,
                16, 1, 3, 14, 18, 10, 16,
            ];
            let ni: [f64; 33] = [
                8.014_969_899_294_95e-2,
                -0.543_862_807_146_111,
                0.337_455_597_421_283,
                8.905_554_511_574_5,
                313.840_736_431_485,
                0.797_367_065_977_789,
                -1.216_169_735_562_4,
                8.728_033_869_374_77,
                -16.976_978_175_760_2,
                -186.552_827_328_416,
                95_115.927_434_423_7,
                -18.916_851_012_049_4,
                -4_334.070_371_948_4,
                543_212_633.012_715,
                0.144_793_408_386_013,
                128.024_559_637_516,
                -67_230.953_407_126_8,
                33_697_238.009_528_7,
                -586.634_196_762_72,
                -22_140_322_476.988_9,
                1_716.066_687_083_89,
                -570_817_595.806_302,
                -3_121.096_931_784_82,
                -2_078_413.846_330_1,
                3_056_059_461_577.86,
                3_221.570_043_143_33,
                326_810_259_797.295,
                -1_441.041_589_344_87,
                410.694_867_802_691,
                109_077_066_873.024,
                -24_796_465_425_889.3,
                1_888_019_068.651_34,
                -123_651_009_018_773.0,
            ];
            let eta = h / 4100.0;
            let sigma = s / 7.9;
            let mut pi = 0.0;
            for i in 0..33 {
                pi += ni[i] * (eta - 0.6).powi(ii[i]) * (sigma - 1.01).powi(ji[i]);
            }
            Ok(pi.powi(4) * 100.0)
        }
        _ => {
            let ii: [i32; 31] = [
                0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 6, 6,
                10, 12, 16,
            ];
            let ji: [i32; 31] = [
                0, 1, 2, 3, 4, 8, 0, 2, 5, 8, 14, 2, 3, 7, 10, 18, 0, 5, 8, 16, 18, 18, 1, 4, 6,
                14, 8, 18, 7, 7, 10,
            ];
            let ni: [f64; 31] = [
                0.112_225_607_199_012,
                -3.390_059_536_067_12,
                -32.050_391_173_009_4,
                -197.597_305_104_9,
                -407.693_861_553_446,
                13_294.377_522_233_1,
                1.708_468_397_740_07,
                37.369_419_814_224_5,
                3_581.443_658_154_34,
                423_014.446_424_664,
                -751_071_025.760_063,
                52.344_612_760_789_8,
                -228.351_290_812_417,
                -960_652.417_056_937,
                -80_705_929.252_607_4,
                1_626_980_172_256.69,
                0.772_465_073_604_171,
                46_392.997_383_774_6,
                -13_731_788.513_412_8,
                1_704_703_926_305.12,
                -25_110_462_818_730.8,
                31_774_883_083_552.0,
                53.868_562_367_531_2,
                -55_308.909_462_516_9,
                -1_028_615.224_214_05,
                2_042_494_187_562.34,
                273_918_446.626_977,
                -2.639_631_463_126_85e15,
                -1_078_908_541.080_88,
                -29_649_262_098.012_4,
                -1.117_549_073_234_24e15,
            ];
            let eta = h / 3500.0;
            let sigma = s / 5.9;
            let mut pi = 0.0;
            for i in 0..31 {
                pi += ni[i] * (eta - 0.7).powi(ii[i]) * (sigma - 1.1).powi(ji[i]);
            }
            Ok(pi.powi(4) * 100.0)
        }
    }
}

pub fn t2_p_rho(p_mpa: f64, rho: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !rho.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if rho <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let mut low = if p_mpa < 16.5292 {
        tsat_k_from_p_mpa(p_mpa)?
    } else {
        b23_t_k_from_p_mpa(p_mpa)?
    };
    let mut high = 1073.15;
    for _ in 0..250 {
        let ts = (low + high) / 2.0;
        let rhos = 1.0 / v2_p_t(p_mpa, ts)?;
        if (rho - rhos).abs() <= 1e-6 {
            return Ok(ts);
        }
        if rhos < rho {
            high = ts;
        } else {
            low = ts;
        }
    }
    Err(If97Error::InvalidIntermediateValue)
}

const R3_I: [i32; 40] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6,
    6, 7, 8, 9, 9, 10, 10, 11,
];
const R3_J: [i32; 40] = [
    0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1,
    3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26,
];
const R3_N: [f64; 40] = [
    1.065_807_002_851_3,
    -15.732_845_290_239,
    20.944_396_974_307,
    -7.686_770_787_871_6,
    2.618_594_778_795_4,
    -2.808_078_114_862,
    1.205_336_969_651_7,
    -8.456_681_281_250_2e-3,
    -1.265_431_547_771_4,
    -1.152_440_780_668_1,
    0.885_210_439_843_18,
    -0.642_077_651_816_07,
    0.384_934_601_866_71,
    -0.852_147_088_242_06,
    4.897_228_154_187_7,
    -3.050_261_725_696_5,
    0.039_420_536_879_154,
    0.125_584_084_243_08,
    -0.279_993_296_987_1,
    1.389_979_956_946,
    -2.018_991_502_357,
    -8.214_763_717_396_3e-3,
    -0.475_960_357_349_23,
    0.043_984_074_473_5,
    -0.444_764_354_287_39,
    0.905_720_707_197_33,
    0.705_224_500_879_67,
    0.107_705_126_263_32,
    -0.329_136_232_589_54,
    -0.508_710_620_411_58,
    -0.022_175_400_873_096,
    0.094_260_751_665_092,
    0.164_362_784_479_61,
    -0.013_503_372_241_348,
    -0.014_834_345_352_472,
    5.792_295_362_808_4e-4,
    3.230_890_470_371_1e-3,
    8.096_480_299_621_5e-5,
    -1.655_767_979_503_7e-4,
    -4.492_389_906_181_5e-5,
];

fn r3_common(rho: f64, t_k: f64) -> If97Result<(f64, f64)> {
    if !rho.is_finite() || !t_k.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if rho <= 0.0 || t_k <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let delta = rho / 322.0;
    let tau = 647.096 / t_k;
    Ok((delta, tau))
}

fn r3_phi(delta: f64, tau: f64) -> If97Result<f64> {
    if delta <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let mut phi = R3_N[0] * delta.ln();
    for i in 1..40 {
        phi += R3_N[i] * delta.powi(R3_I[i]) * tau.powi(R3_J[i]);
    }
    if !phi.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(phi)
}

fn r3_phi_delta(delta: f64, tau: f64) -> If97Result<f64> {
    if delta <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let mut phi_delta = R3_N[0] / delta;
    for i in 1..40 {
        let ii = R3_I[i];
        if ii == 0 {
            continue;
        }
        phi_delta += R3_N[i] * (ii as f64) * delta.powi(ii - 1) * tau.powi(R3_J[i]);
    }
    if !phi_delta.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(phi_delta)
}

fn r3_phi_tau(delta: f64, tau: f64) -> If97Result<f64> {
    let mut phi_tau = 0.0;
    for i in 1..40 {
        let jj = R3_J[i];
        if jj == 0 {
            continue;
        }
        phi_tau += R3_N[i] * delta.powi(R3_I[i]) * (jj as f64) * tau.powi(jj - 1);
    }
    if !phi_tau.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(phi_tau)
}

fn r3_phi_tautau(delta: f64, tau: f64) -> If97Result<f64> {
    let mut phi_tautau = 0.0;
    for i in 0..40 {
        let jj = R3_J[i];
        if jj < 2 {
            continue;
        }
        phi_tautau +=
            R3_N[i] * delta.powi(R3_I[i]) * (jj as f64) * ((jj - 1) as f64) * tau.powi(jj - 2);
    }
    if !phi_tautau.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(phi_tautau)
}

fn r3_phi_deltadelta(delta: f64, tau: f64) -> If97Result<f64> {
    if delta <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let mut phi_deltadelta = -R3_N[0] / (delta * delta);
    for i in 1..40 {
        let ii = R3_I[i];
        if ii < 2 {
            continue;
        }
        phi_deltadelta +=
            R3_N[i] * (ii as f64) * ((ii - 1) as f64) * delta.powi(ii - 2) * tau.powi(R3_J[i]);
    }
    if !phi_deltadelta.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(phi_deltadelta)
}

fn r3_phi_deltatau(delta: f64, tau: f64) -> If97Result<f64> {
    if delta <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let mut phi_deltatau = 0.0;
    for i in 1..40 {
        let ii = R3_I[i];
        let jj = R3_J[i];
        if ii == 0 || jj == 0 {
            continue;
        }
        phi_deltatau += R3_N[i] * (ii as f64) * delta.powi(ii - 1) * (jj as f64) * tau.powi(jj - 1);
    }
    if !phi_deltatau.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(phi_deltatau)
}

pub fn p3_rho_t(rho: f64, t_k: f64) -> If97Result<f64> {
    let (delta, tau) = r3_common(rho, t_k)?;
    let phi_delta = r3_phi_delta(delta, tau)?;
    let p = rho * R_KJ_KG_K * t_k * delta * phi_delta / 1000.0;
    if !p.is_finite() || p <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(p)
}

pub fn u3_rho_t(rho: f64, t_k: f64) -> If97Result<f64> {
    let (delta, tau) = r3_common(rho, t_k)?;
    let phi_tau = r3_phi_tau(delta, tau)?;
    let u = R_KJ_KG_K * t_k * (tau * phi_tau);
    if !u.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(u)
}

pub fn h3_rho_t(rho: f64, t_k: f64) -> If97Result<f64> {
    let (delta, tau) = r3_common(rho, t_k)?;
    let phi_delta = r3_phi_delta(delta, tau)?;
    let phi_tau = r3_phi_tau(delta, tau)?;
    let h = R_KJ_KG_K * t_k * (tau * phi_tau + delta * phi_delta);
    if !h.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(h)
}

pub fn s3_rho_t(rho: f64, t_k: f64) -> If97Result<f64> {
    let (delta, tau) = r3_common(rho, t_k)?;
    let phi = r3_phi(delta, tau)?;
    let phi_tau = r3_phi_tau(delta, tau)?;
    let s = R_KJ_KG_K * (tau * phi_tau - phi);
    if !s.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(s)
}

pub fn cp3_rho_t(rho: f64, t_k: f64) -> If97Result<f64> {
    let (delta, tau) = r3_common(rho, t_k)?;
    let phi_tautau = r3_phi_tautau(delta, tau)?;
    let phi_delta = r3_phi_delta(delta, tau)?;
    let phi_deltatau = r3_phi_deltatau(delta, tau)?;
    let phi_deltadelta = r3_phi_deltadelta(delta, tau)?;
    let denom = 2.0 * delta * phi_delta + delta * delta * phi_deltadelta;
    if denom == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let num = (delta * phi_delta - delta * tau * phi_deltatau).powi(2);
    let cp = R_KJ_KG_K * (-tau * tau * phi_tautau + num / denom);
    if !cp.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(cp)
}

pub fn cv3_rho_t(rho: f64, t_k: f64) -> If97Result<f64> {
    let (delta, tau) = r3_common(rho, t_k)?;
    let phi_tautau = r3_phi_tautau(delta, tau)?;
    let cv = -R_KJ_KG_K * tau * tau * phi_tautau;
    if !cv.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(cv)
}

pub fn w3_rho_t(rho: f64, t_k: f64) -> If97Result<f64> {
    let (delta, tau) = r3_common(rho, t_k)?;
    let phi_tautau = r3_phi_tautau(delta, tau)?;
    let phi_delta = r3_phi_delta(delta, tau)?;
    let phi_deltatau = r3_phi_deltatau(delta, tau)?;
    let phi_deltadelta = r3_phi_deltadelta(delta, tau)?;
    let denom = tau * tau * phi_tautau;
    if denom == 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let term = (delta * phi_delta - delta * tau * phi_deltatau).powi(2) / denom;
    let w2 = 1000.0
        * R_KJ_KG_K
        * t_k
        * (2.0 * delta * phi_delta + delta * delta * phi_deltadelta - term);
    let w = sqrt_nonneg(w2)?;
    if !w.is_finite() {
        return Err(If97Error::InvalidIntermediateValue);
    }
    Ok(w)
}

pub fn rho3_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !t_k.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if p_mpa <= 0.0 || t_k <= 0.0 {
        return Err(If97Error::InvalidIntermediateValue);
    }
    let mut best_rho = None;
    let mut best_abs_err = f64::INFINITY;
    let mut bracket = None;
    let mut prev = None;

    let mut sample_rhos = Vec::with_capacity(260);
    sample_rhos.extend_from_slice(&[1e-6, 1e-4, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]);
    let mut rho = 20.0;
    while rho <= 2_000.0 {
        sample_rhos.push(rho);
        rho += 10.0;
    }

    for &rho in &sample_rhos {
        let p = match p3_rho_t(rho, t_k) {
            Ok(p) => p,
            Err(_) => continue,
        };
        let f = p - p_mpa;
        let abs_f = f.abs();
        if abs_f < best_abs_err {
            best_abs_err = abs_f;
            best_rho = Some(rho);
        }
        if let Some((rho_prev, f_prev)) = prev {
            if f_prev == 0.0 {
                return Ok(rho_prev);
            }
            if f == 0.0 {
                return Ok(rho);
            }
            if (f_prev < 0.0 && f > 0.0) || (f_prev > 0.0 && f < 0.0) {
                bracket = Some((rho_prev, f_prev, rho));
                break;
            }
        }
        prev = Some((rho, f));
    }

    let tol_p = 1e-12 * (1.0 + p_mpa.abs());
    if let Some((mut low, mut f_low, mut high)) = bracket {
        for _ in 0..250 {
            let mid = (low + high) / 2.0;
            let p_mid = p3_rho_t(mid, t_k)?;
            let f_mid = p_mid - p_mpa;
            if f_mid.abs() <= tol_p {
                return Ok(mid);
            }
            if (f_low < 0.0 && f_mid > 0.0) || (f_low > 0.0 && f_mid < 0.0) {
                high = mid;
            } else {
                low = mid;
                f_low = f_mid;
            }
            if (high - low).abs() <= 1e-12 * (1.0 + mid.abs()) {
                return Ok((low + high) / 2.0);
            }
        }
        return Err(If97Error::InvalidIntermediateValue);
    }

    let mut rho = best_rho.ok_or(If97Error::InvalidIntermediateValue)?;
    for _ in 0..60 {
        let (delta, tau) = r3_common(rho, t_k)?;
        let phi_delta = r3_phi_delta(delta, tau)?;
        let phi_deltadelta = r3_phi_deltadelta(delta, tau)?;
        let p = rho * R_KJ_KG_K * t_k * delta * phi_delta / 1000.0;
        let f = p - p_mpa;
        if f.abs() <= tol_p {
            return Ok(rho);
        }
        let dp_drho =
            (R_KJ_KG_K * t_k / 1000.0) * (2.0 * delta * phi_delta + delta * delta * phi_deltadelta);
        if !dp_drho.is_finite() || dp_drho == 0.0 {
            break;
        }
        let mut step = f / dp_drho;
        if step.abs() > 0.5 * rho {
            step = step.signum() * 0.5 * rho;
        }
        let next_rho = rho - step;
        if !next_rho.is_finite() || next_rho <= 0.0 {
            break;
        }
        rho = next_rho;
    }

    Err(If97Error::InvalidIntermediateValue)
}

pub fn v3_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let rho = rho3_p_t(p_mpa, t_k)?;
    Ok(1.0 / rho)
}

pub fn h3_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let rho = rho3_p_t(p_mpa, t_k)?;
    h3_rho_t(rho, t_k)
}

pub fn u3_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let rho = rho3_p_t(p_mpa, t_k)?;
    u3_rho_t(rho, t_k)
}

pub fn s3_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let rho = rho3_p_t(p_mpa, t_k)?;
    s3_rho_t(rho, t_k)
}

pub fn cp3_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let rho = rho3_p_t(p_mpa, t_k)?;
    cp3_rho_t(rho, t_k)
}

pub fn cv3_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let rho = rho3_p_t(p_mpa, t_k)?;
    cv3_rho_t(rho, t_k)
}

pub fn w3_p_t(p_mpa: f64, t_k: f64) -> If97Result<f64> {
    let rho = rho3_p_t(p_mpa, t_k)?;
    w3_rho_t(rho, t_k)
}

pub fn t3_p_h(p_mpa: f64, h: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !h.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let h3ab = 2_014.640_042_068_75 + 3.746_965_501_369_83 * p_mpa
        - 2.199_219_010_541_87e-2 * p_mpa * p_mpa
        + 8.751_316_860_099_5e-5 * p_mpa * p_mpa * p_mpa;
    if h < h3ab {
        let ii: [i32; 31] = [
            -12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -8, -8, -8, -8, -5, -3, -2, -2,
            -2, -1, -1, 0, 0, 1, 3, 3, 4, 4, 10, 12,
        ];
        let ji: [i32; 31] = [
            0, 1, 2, 6, 14, 16, 20, 22, 1, 5, 12, 0, 2, 4, 10, 2, 0, 1, 3, 4, 0, 2, 0, 1, 1, 0, 1,
            0, 3, 4, 5,
        ];
        let ni: [f64; 31] = [
            -1.336_456_678_112_15e-7,
            4.559_126_568_029_78e-6,
            -1.462_946_407_009_79e-5,
            6.393_413_129_700_8e-3,
            372.783_927_268_847,
            -7_186.543_774_604_47,
            573_494.752_103_4,
            -2_675_693.291_114_39,
            -3.340_662_833_026_14e-5,
            -2.454_792_140_695_97e-2,
            47.808_784_776_499_6,
            7.646_641_318_189_04e-6,
            1.283_506_276_769_72e-3,
            1.712_190_813_773_31e-2,
            -8.510_073_045_832_13,
            -1.365_134_616_297_81e-2,
            -3.844_609_975_966_57e-6,
            3.374_238_079_116_55e-3,
            -0.551_624_873_066_791,
            0.729_202_277_107_47,
            -9.925_227_573_760_41e-3,
            -0.119_308_831_407_288,
            0.793_929_190_615_421,
            0.454_270_731_799_386,
            0.209_998_591_259_91,
            -6.421_098_239_047_38e-3,
            -0.023_515_586_860_454,
            2.522_331_083_416_12e-3,
            -7.648_851_333_681_19e-3,
            1.361_764_275_742_91e-2,
            -1.330_278_835_756_69e-2,
        ];
        let ps = p_mpa / 100.0;
        let hs = h / 2300.0;
        let mut ts = 0.0;
        for i in 0..31 {
            ts += ni[i] * (ps + 0.24).powi(ii[i]) * (hs - 0.615).powi(ji[i]);
        }
        let t = ts * 760.0;
        if !t.is_finite() {
            return Err(If97Error::InvalidIntermediateValue);
        }
        Ok(t)
    } else {
        let ii: [i32; 33] = [
            -12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2,
            -1, -1, -1, -1, -1, -1, 0, 0, 1, 3, 5, 6, 8,
        ];
        let ji: [i32; 33] = [
            0, 1, 0, 1, 5, 10, 12, 0, 1, 2, 4, 10, 0, 1, 2, 0, 1, 5, 0, 4, 2, 4, 6, 10, 14, 16, 0,
            2, 1, 1, 1, 1, 1,
        ];
        let ni: [f64; 33] = [
            3.232_545_736_449_2e-5,
            -1.275_755_565_871_81e-4,
            -4.758_518_773_560_68e-4,
            1.561_830_141_816_02e-3,
            0.105_724_860_113_781,
            -85.851_422_113_253_4,
            724.140_095_480_911,
            2.964_758_102_732_57e-3,
            -5.927_219_833_659_88e-3,
            -1.263_054_228_186_66e-2,
            -0.115_716_196_364_853,
            84.900_096_973_959_5,
            -1.086_022_600_866_15e-2,
            1.543_044_753_288_51e-2,
            7.504_554_415_244_66e-2,
            2.525_209_736_129_82e-2,
            -6.025_079_012_329_96e-2,
            -3.076_222_213_505_01,
            -5.740_119_598_648_79e-2,
            5.034_713_609_398_49,
            -0.925_081_888_584_834,
            3.917_338_829_175_46,
            -77.314_600_713_019,
            9_493.087_620_985_87,
            -1_410_437.196_794_09,
            8_491_662.308_190_26,
            0.861_095_729_446_704,
            0.323_346_442_811_72,
            0.873_281_936_020_439,
            -0.436_653_048_526_683,
            0.286_596_714_529_479,
            -0.131_778_331_276_228,
            6.766_820_643_302_75e-3,
        ];
        let hs = h / 2800.0;
        let ps = p_mpa / 100.0;
        let mut ts = 0.0;
        for i in 0..33 {
            ts += ni[i] * (ps + 0.298).powi(ii[i]) * (hs - 0.72).powi(ji[i]);
        }
        let t = ts * 860.0;
        if !t.is_finite() {
            return Err(If97Error::InvalidIntermediateValue);
        }
        Ok(t)
    }
}

pub fn v3_p_h(p_mpa: f64, h: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !h.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    let h3ab = 2_014.640_042_068_75 + 3.746_965_501_369_83 * p_mpa
        - 2.199_219_010_541_87e-2 * p_mpa * p_mpa
        + 8.751_316_860_099_5e-5 * p_mpa * p_mpa * p_mpa;
    if h < h3ab {
        let ii: [i32; 32] = [
            -12, -12, -12, -12, -10, -10, -10, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1,
            -1, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8,
        ];
        let ji: [i32; 32] = [
            6, 8, 12, 18, 4, 7, 10, 5, 12, 3, 4, 22, 2, 3, 7, 3, 16, 0, 1, 2, 3, 0, 1, 0, 1, 2, 0,
            2, 0, 2, 2, 2,
        ];
        let ni: [f64; 32] = [
            5.299_440_629_660_28e-3,
            -0.170_099_690_234_461,
            11.132_381_431_292_7,
            -2_178.981_231_451_25,
            -5.060_618_279_808_75e-4,
            0.556_495_239_685_324,
            -9.436_727_260_940_16,
            -0.297_856_807_561_527,
            93.935_394_371_718_6,
            1.929_449_394_659_81e-2,
            0.421_740_664_704_763,
            -3_689_141.262_823_3,
            -7.375_668_476_006_39e-3,
            -0.354_753_242_424_366,
            -1.997_681_693_387_27,
            1.154_562_970_590_49,
            5_683.668_758_159_6,
            8.081_695_401_246_68e-3,
            0.172_416_341_519_307,
            1.042_701_752_929_27,
            -0.297_691_372_792_847,
            0.560_394_465_163_593,
            0.275_234_661_176_914,
            -0.148_347_894_866_012,
            -6.511_425_134_785_15e-2,
            -2.924_687_153_863_02,
            6.648_760_969_526_65e-2,
            3.523_350_142_638_44,
            -1.463_407_923_133_32e-2,
            -2.245_034_866_681_84,
            1.105_334_647_061_42,
            -4.087_573_444_956_12e-2,
        ];
        let ps = p_mpa / 100.0;
        let hs = h / 2100.0;
        let mut vs = 0.0;
        for i in 0..32 {
            vs += ni[i] * (ps + 0.128).powi(ii[i]) * (hs - 0.727).powi(ji[i]);
        }
        let v = vs * 0.0028;
        if !v.is_finite() || v <= 0.0 {
            return Err(If97Error::InvalidIntermediateValue);
        }
        Ok(v)
    } else {
        let ii: [i32; 30] = [
            -12, -12, -8, -8, -8, -8, -8, -8, -6, -6, -6, -6, -6, -6, -4, -4, -4, -3, -3, -2, -2,
            -1, -1, -1, -1, 0, 1, 1, 2, 2,
        ];
        let ji: [i32; 30] = [
            0, 1, 0, 1, 3, 6, 7, 8, 0, 1, 2, 5, 6, 10, 3, 6, 10, 0, 2, 1, 2, 0, 1, 4, 5, 0, 0, 1,
            2, 6,
        ];
        let ni: [f64; 30] = [
            -2.251_969_343_363_18e-9,
            1.406_743_633_134_86e-8,
            2.337_840_852_805_6e-6,
            -3.318_337_152_290_01e-5,
            1.079_567_785_143_18e-3,
            -0.271_382_067_378_863,
            1.072_022_624_903_33,
            -0.853_821_329_075_382,
            -2.152_141_943_405_26e-5,
            7.696_560_882_227_3e-4,
            -4.311_365_804_338_64e-3,
            0.453_342_167_309_331,
            -0.507_749_535_873_652,
            -100.475_154_528_389,
            -0.219_201_924_648_793,
            -3.210_879_656_689_17,
            607.567_815_637_771,
            5.576_864_506_859_32e-4,
            0.187_499_040_029_55,
            9.053_680_304_481_07e-3,
            0.285_417_173_048_685,
            3.299_240_309_960_98e-2,
            0.239_897_419_685_483,
            4.827_549_959_513_94,
            -11.803_575_370_223_1,
            0.169_490_044_091_791,
            -1.799_672_225_077_87e-2,
            3.718_101_163_326_74e-2,
            -5.362_883_350_650_96e-2,
            1.606_971_010_925_2,
        ];
        let ps = p_mpa / 100.0;
        let hs = h / 2800.0;
        let mut vs = 0.0;
        for i in 0..30 {
            vs += ni[i] * (ps + 0.0661).powi(ii[i]) * (hs - 0.72).powi(ji[i]);
        }
        let v = vs * 0.0088;
        if !v.is_finite() || v <= 0.0 {
            return Err(If97Error::InvalidIntermediateValue);
        }
        Ok(v)
    }
}

pub fn t3_p_s(p_mpa: f64, s: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if s <= 4.412_021_482_234_76 {
        let ii: [i32; 33] = [
            -12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -6, -6, -5, -5, -5, -4, -4, -4, -2,
            -2, -1, -1, 0, 0, 0, 1, 2, 2, 3, 8, 8, 10,
        ];
        let ji: [i32; 33] = [
            28, 32, 4, 10, 12, 14, 5, 7, 8, 28, 2, 6, 32, 0, 14, 32, 6, 10, 36, 1, 4, 1, 6, 0, 1,
            4, 0, 0, 3, 2, 0, 1, 2,
        ];
        let ni: [f64; 33] = [
            1_500_420_082.638_75,
            -159_397_258_480.424,
            5.021_811_402_179_75e-4,
            -67.205_776_785_546_6,
            1_450.585_454_044_56,
            -8_238.895_348_888_9,
            -0.154_852_214_233_853,
            11.230_504_674_669_5,
            -29.700_021_348_282_2,
            43_856_513_263.549_5,
            1.378_378_386_354_64e-3,
            -2.974_785_271_574_62,
            9_717_779_473_494.13,
            -5.715_277_670_523_98e-5,
            28_830.794_977_842,
            -74_442_828_926_270.3,
            12.801_732_484_892_1,
            -368.275_545_889_071,
            6.647_689_047_791_77e15,
            0.044_935_925_195_888,
            -4.228_978_360_996_55,
            -0.240_614_376_434_179,
            -4.743_413_652_549_24,
            0.724_093_999_126_11,
            0.923_874_349_695_897,
            3.990_436_552_810_15,
            3.840_666_518_680_09e-2,
            -3.593_443_655_718_48e-3,
            -0.735_196_448_821_653,
            0.188_367_048_396_131,
            1.410_642_668_187_04e-4,
            -2.574_185_014_963_37e-3,
            1.232_200_248_515_55e-3,
        ];
        let sigma = s / 4.4;
        let pi = p_mpa / 100.0;
        let mut teta = 0.0;
        for i in 0..33 {
            teta += ni[i] * (pi + 0.24).powi(ii[i]) * (sigma - 0.703).powi(ji[i]);
        }
        Ok(teta * 760.0)
    } else {
        let ii: [i32; 28] = [
            -12, -12, -12, -12, -8, -8, -8, -6, -6, -6, -5, -5, -5, -5, -5, -4, -3, -3, -2, 0, 2,
            3, 4, 5, 6, 8, 12, 14,
        ];
        let ji: [i32; 28] = [
            1, 3, 4, 7, 0, 1, 3, 0, 2, 4, 0, 1, 2, 4, 6, 12, 1, 6, 2, 0, 1, 1, 0, 24, 0, 3, 1, 2,
        ];
        let ni: [f64; 28] = [
            0.527_111_701_601_66,
            -40.131_783_005_274_2,
            153.020_073_134_484,
            -2_247.993_982_188_27,
            -0.193_993_484_669_048,
            -1.404_675_578_937_68,
            42.679_987_811_402_4,
            0.752_810_643_416_743,
            22.665_723_861_641_7,
            -622.873_556_909_932,
            -0.660_823_667_935_396,
            0.841_267_087_271_658,
            -25.371_750_176_439_7,
            485.708_963_532_948,
            880.531_517_490_555,
            2_650_155.927_946_26,
            -0.359_287_150_025_783,
            -656.991_567_673_753,
            2.417_681_491_853_67,
            0.856_873_461_222_588,
            0.655_143_675_313_458,
            -0.213_535_213_206_406,
            5.629_749_576_063_48e-3,
            -316_955_725_450_471.0,
            -6.999_970_001_524_57e-4,
            1.198_458_032_107_67e-2,
            1.938_481_220_220_95e-5,
            -2.150_957_491_823_09e-5,
        ];
        let sigma = s / 5.3;
        let pi = p_mpa / 100.0;
        let mut teta = 0.0;
        for i in 0..28 {
            teta += ni[i] * (pi + 0.76).powi(ii[i]) * (sigma - 0.818).powi(ji[i]);
        }
        Ok(teta * 860.0)
    }
}

pub fn v3_p_s(p_mpa: f64, s: f64) -> If97Result<f64> {
    if !p_mpa.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if s <= 4.412_021_482_234_76 {
        let ii: [i32; 28] = [
            -12, -12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -5, -4, -3, -3, -2, -2, -1, -1,
            0, 0, 0, 1, 2, 4, 5, 6,
        ];
        let ji: [i32; 28] = [
            10, 12, 14, 4, 8, 10, 20, 5, 6, 14, 16, 28, 1, 5, 2, 4, 3, 8, 1, 2, 0, 1, 3, 0, 0, 2,
            2, 0,
        ];
        let ni: [f64; 28] = [
            79.554_407_409_397_5,
            -2_382.612_429_845_9,
            17_681.310_061_778_7,
            -1.105_247_270_803_79e-3,
            -15.321_383_365_532_6,
            297.544_599_376_982,
            -35_031_520.687_124_2,
            0.277_513_761_062_119,
            -0.523_964_271_036_888,
            -148_011.182_995_403,
            1_600_148.993_742_66,
            1_708_023_226_634.27,
            2.468_669_960_064_94e-4,
            1.653_260_847_979_8,
            -0.118_008_384_666_987,
            2.537_986_423_559,
            0.965_127_704_669_424,
            -28.217_242_053_282_6,
            0.203_224_612_353_823,
            1.106_481_860_635_13,
            0.526_127_948_451_28,
            0.277_000_018_736_321,
            1.081_533_405_011_32,
            -7.441_278_853_578_93e-2,
            1.640_944_435_413_84e-2,
            -6.804_682_753_010_65e-2,
            0.025_798_857_610_164,
            -1.457_498_619_444_16e-4,
        ];
        let pi = p_mpa / 100.0;
        let sigma = s / 4.4;
        let mut omega = 0.0;
        for i in 0..28 {
            omega += ni[i] * (pi + 0.187).powi(ii[i]) * (sigma - 0.755).powi(ji[i]);
        }
        Ok(omega * 0.0028)
    } else {
        let ii: [i32; 31] = [
            -12, -12, -12, -12, -12, -12, -10, -10, -10, -10, -8, -5, -5, -5, -4, -4, -4, -4, -3,
            -2, -2, -2, -2, -2, -2, 0, 0, 0, 1, 1, 2,
        ];
        let ji: [i32; 31] = [
            0, 1, 2, 3, 5, 6, 0, 1, 2, 4, 0, 1, 2, 3, 0, 1, 2, 3, 1, 0, 1, 2, 3, 4, 12, 0, 1, 2, 0,
            2, 2,
        ];
        let ni: [f64; 31] = [
            5.915_997_803_222_38e-5,
            -1.854_659_971_378_56e-3,
            1.041_905_104_800_13e-2,
            5.986_473_020_385_9e-3,
            -0.771_391_189_901_699,
            1.725_497_655_570_36,
            -4.670_760_798_465_26e-4,
            1.345_338_233_844_39e-2,
            -8.080_943_368_054_95e-2,
            0.508_139_374_365_767,
            1.285_846_433_616_83e-3,
            -1.638_993_539_154_35,
            5.869_381_993_180_63,
            -2.924_666_679_186_13,
            -6.140_763_014_995_37e-3,
            5.761_990_140_491_72,
            -12.161_332_060_678_8,
            1.676_375_409_579_44,
            -7.441_358_387_734_63,
            3.781_680_914_376_59e-2,
            4.014_322_030_276_88,
            16.027_983_747_918_5,
            3.178_487_793_477_28,
            -3.583_623_103_048_53,
            -1_159_952.604_468_27,
            0.199_256_573_577_909,
            -0.122_270_624_794_624,
            -19.144_914_371_658_6,
            -1.504_480_029_052_84e-2,
            14.640_790_016_215_4,
            -3.274_777_871_882_3,
        ];
        let pi = p_mpa / 100.0;
        let sigma = s / 5.3;
        let mut omega = 0.0;
        for i in 0..31 {
            omega += ni[i] * (pi + 0.298).powi(ii[i]) * (sigma - 0.816).powi(ji[i]);
        }
        Ok(omega * 0.0088)
    }
}

pub fn p3_h_s(h: f64, s: f64) -> If97Result<f64> {
    if !h.is_finite() || !s.is_finite() {
        return Err(If97Error::NonFiniteInput);
    }
    if s < 4.412_021_482_234_76 {
        let ii: [i32; 33] = [
            0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 6, 7, 8, 10, 10, 14, 18, 20, 22,
            22, 24, 28, 28, 32, 32,
        ];
        let ji: [i32; 33] = [
            0, 1, 5, 0, 3, 4, 8, 14, 6, 16, 0, 2, 3, 0, 1, 4, 5, 28, 28, 24, 1, 32, 36, 22, 28, 36,
            16, 28, 36, 16, 36, 10, 28,
        ];
        let ni: [f64; 33] = [
            7.708_898_283_269_34,
            -26.083_500_912_868_8,
            267.416_218_930_389,
            17.222_108_949_684_4,
            -293.542_332_145_97,
            614.135_601_882_478,
            -61_056.275_772_567_4,
            -65_127_225.111_821_9,
            73_591.931_352_193_7,
            -11_664_650_591.419_1,
            35.526_708_643_446_1,
            -596.144_543_825_955,
            -475.842_430_145_708,
            69.678_196_535_950_3,
            335.674_250_377_312,
            25_052.680_913_088_2,
            146_997.380_630_766,
            5.380_693_150_915_34e19,
            1.436_198_272_913_46e21,
            3.649_858_661_659_94e19,
            -2_547.415_611_567_75,
            2.401_201_970_965_63e27,
            -3.938_474_646_794_96e29,
            1.470_734_070_248_52e24,
            -4.263_912_504_320_59e31,
            1.945_093_406_210_77e38,
            6.662_121_321_148_96e23,
            7.067_770_165_528_58e33,
            1.755_636_219_755_76e41,
            1.084_086_074_291_24e28,
            7.308_727_051_751_51e43,
            1.591_458_473_988_7e24,
            3.771_216_059_433_24e40,
        ];
        let sigma = s / 4.4;
        let eta = h / 2300.0;
        let mut pi = 0.0;
        for i in 0..33 {
            pi += ni[i] * (eta - 1.01).powi(ii[i]) * (sigma - 0.75).powi(ji[i]);
        }
        Ok(pi * 99.0)
    } else {
        let ii: [i32; 35] = [
            -12, -12, -12, -12, -12, -10, -10, -10, -10, -8, -8, -6, -6, -6, -6, -5, -4, -4, -4,
            -3, -3, -3, -3, -2, -2, -1, 0, 2, 2, 5, 6, 8, 10, 14, 14,
        ];
        let ji: [i32; 35] = [
            2, 10, 12, 14, 20, 2, 10, 14, 18, 2, 8, 2, 6, 7, 8, 10, 4, 5, 8, 1, 3, 5, 6, 0, 1, 0,
            3, 0, 1, 0, 1, 1, 1, 3, 7,
        ];
        let ni: [f64; 35] = [
            1.252_443_607_179_79e-13,
            -1.265_993_225_537_13e-2,
            5.068_780_301_406_26,
            31.784_717_115_420_2,
            -391_041.161_399_932,
            -9.757_334_063_920_44e-11,
            -18.631_241_948_827_9,
            510.973_543_414_101,
            373_847.005_822_362,
            2.998_040_246_665_72e-8,
            20.054_439_382_034_2,
            -4.980_304_876_628_29e-6,
            -10.230_180_636_003,
            55.281_912_699_032_5,
            -206.211_367_510_878,
            -7_940.122_323_248_23,
            7.822_484_720_281_53,
            -58.654_432_690_246_8,
            3_550.736_476_964_81,
            -1.153_031_072_901_62e-4,
            -1.750_924_031_718_02,
            257.981_687_748_16,
            -727.048_374_179_467,
            1.216_448_226_091_98e-4,
            3.931_378_717_626_92e-2,
            7.041_810_059_092_96e-3,
            -82.910_820_069_811,
            -0.265_178_818_131_25,
            13.753_168_245_399_1,
            -52.239_409_075_304_6,
            2_405.562_989_410_48,
            -22_736.163_126_892_9,
            89_074.634_393_256_7,
            -23_923_456.582_248_6,
            5_687_958_081.297_14,
        ];
        let sigma = s / 5.3;
        let eta = h / 2800.0;
        let mut pi = 0.0;
        for i in 0..35 {
            pi += ni[i] * (eta - 0.681).powi(ii[i]) * (sigma - 0.792).powi(ji[i]);
        }
        Ok(16.6 / pi)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64, tol: f64) {
        let err = (actual - expected).abs();
        assert!(
            err <= tol,
            "actual={actual}, expected={expected}, abs_err={err}, tol={tol}"
        );
    }

    fn assert_rel_close(actual: f64, expected: f64, rel: f64) {
        if expected == 0.0 {
            assert_close(actual, expected, rel);
            return;
        }
        let err = ((actual - expected) / expected).abs();
        assert!(
            err <= rel,
            "actual={actual}, expected={expected}, rel_err={err}, rel_tol={rel}"
        );
    }

    #[test]
    fn region4_psat_reference_points() {
        assert_close(psat_mpa_from_t_k(300.0).unwrap(), 0.003_536_589_41, 5e-12);
        assert_close(psat_mpa_from_t_k(500.0).unwrap(), 2.638_897_76, 1e-8);
        assert_close(psat_mpa_from_t_k(600.0).unwrap(), 12.344_314_6, 1e-7);
    }

    #[test]
    fn region4_tsat_reference_points() {
        assert_close(tsat_k_from_p_mpa(0.1).unwrap(), 372.755_919, 1e-6);
        assert_close(tsat_k_from_p_mpa(1.0).unwrap(), 453.035_632, 1e-6);
        assert_close(tsat_k_from_p_mpa(10.0).unwrap(), 584.149_488, 1e-6);
    }

    #[test]
    fn region4_roundtrip() {
        let temps = [300.0, 350.0, 373.15, 500.0, 640.0];
        for &t in &temps {
            let p = psat_mpa_from_t_k(t).unwrap();
            let t2 = tsat_k_from_p_mpa(p).unwrap();
            assert_close(t2, t, 2e-9);
        }
    }

    #[test]
    fn region1_table5() {
        let p = [30.0 / 10.0, 800.0 / 10.0, 30.0 / 10.0];
        let t = [300.0, 300.0, 500.0];
        let if97_v = [0.001_002_151_68, 0.000_971_180_894, 0.001_202_418];
        let if97_h = [115.331_273, 184.142_828, 975.542_239];
        let if97_u = [112.324_818, 106.448_356, 971.934_985];
        let if97_s = [0.392_294_792, 0.368_563_852, 2.580_419_12];
        let if97_cp = [4.173_012_18, 4.010_089_87, 4.655_806_82];
        let if97_w = [1_507.739_21, 1_634.690_54, 1_240.713_37];

        for i in 0..3 {
            assert_rel_close(v1_p_t(p[i], t[i]).unwrap(), if97_v[i], 1e-8);
            assert_rel_close(h1_p_t(p[i], t[i]).unwrap(), if97_h[i], 1e-8);
            assert_rel_close(u1_p_t(p[i], t[i]).unwrap(), if97_u[i], 1e-8);
            assert_rel_close(s1_p_t(p[i], t[i]).unwrap(), if97_s[i], 1e-8);
            assert_rel_close(cp1_p_t(p[i], t[i]).unwrap(), if97_cp[i], 1e-8);
            assert_rel_close(w1_p_t(p[i], t[i]).unwrap(), if97_w[i], 1e-8);
        }
    }

    #[test]
    fn region1_backward_tables() {
        let p = [30.0 / 10.0, 800.0 / 10.0, 800.0 / 10.0];
        let h = [500.0, 500.0, 1500.0];
        let if97_t_ph = [391.798_509, 378.108_626, 611.041_229];
        for i in 0..3 {
            assert_rel_close(t1_p_h(p[i], h[i]).unwrap(), if97_t_ph[i], 1e-8);
        }

        let s = [0.5, 0.5, 3.0];
        let if97_t_ps = [307.842_258, 309.979_785, 565.899_909];
        for i in 0..3 {
            assert_rel_close(t1_p_s(p[i], s[i]).unwrap(), if97_t_ps[i], 1e-8);
        }

        let h2 = [0.001, 90.0, 1500.0];
        let s2 = [0.0, 0.0, 3.4];
        let if97_p_hs = [0.000_980_098_061_2, 91.929_547_272, 58.682_944_23];
        for i in 0..3 {
            assert_rel_close(p1_h_s(h2[i], s2[i]).unwrap(), if97_p_hs[i], 1e-8);
        }
    }

    #[test]
    fn region2_table15() {
        let p = [0.035 / 10.0, 0.035 / 10.0, 300.0 / 10.0];
        let t = [300.0, 700.0, 700.0];
        let if97_v = [39.491_386_6, 92.301_589_8, 0.005_429_466_19];
        let if97_h = [2_549.911_45, 3_335.683_75, 2_631.494_74];
        let if97_u = [2_411.691_6, 3_012.628_19, 2_468.610_76];
        let if97_s = [8.522_389_67, 10.174_999_6, 5.175_402_98];
        let if97_cp = [1.913_001_62, 2.081_412_74, 10.350_509_2];
        let if97_w = [427.920_172, 644.289_068, 480.386_523];

        for i in 0..3 {
            assert_rel_close(v2_p_t(p[i], t[i]).unwrap(), if97_v[i], 1e-8);
            assert_rel_close(h2_p_t(p[i], t[i]).unwrap(), if97_h[i], 1e-8);
            assert_rel_close(u2_p_t(p[i], t[i]).unwrap(), if97_u[i], 1e-8);
            assert_rel_close(s2_p_t(p[i], t[i]).unwrap(), if97_s[i], 1e-8);
            assert_rel_close(cp2_p_t(p[i], t[i]).unwrap(), if97_cp[i], 1e-8);
            assert_rel_close(w2_p_t(p[i], t[i]).unwrap(), if97_w[i], 1e-8);
        }
    }

    #[test]
    fn region2_backward_tables() {
        let p = [
            0.01 / 10.0,
            30.0 / 10.0,
            30.0 / 10.0,
            50.0 / 10.0,
            50.0 / 10.0,
            250.0 / 10.0,
            400.0 / 10.0,
            600.0 / 10.0,
            600.0 / 10.0,
        ];
        let h = [
            3000.0, 3000.0, 4000.0, 3500.0, 4000.0, 3500.0, 2700.0, 2700.0, 3200.0,
        ];
        let if97_t_ph = [
            534.433_241,
            575.373_37,
            1_010.775_77,
            801.299_102,
            1_015.315_83,
            875.279_054,
            743.056_411,
            791.137_067,
            882.756_86,
        ];
        for i in 0..9 {
            assert_rel_close(t2_p_h(p[i], h[i]).unwrap(), if97_t_ph[i], 1e-8);
        }

        let p2 = [
            1.0 / 10.0,
            1.0 / 10.0,
            25.0 / 10.0,
            80.0 / 10.0,
            80.0 / 10.0,
            900.0 / 10.0,
            200.0 / 10.0,
            800.0 / 10.0,
            800.0 / 10.0,
        ];
        let s2 = [7.5, 8.0, 8.0, 6.0, 7.5, 6.0, 5.75, 5.25, 5.75];
        let if97_t_ps = [
            399.517_097,
            514.127_081,
            1_039.849_17,
            600.484_04,
            1_064.955_56,
            1_038.011_26,
            697.992_849,
            854.011_484,
            949.017_998,
        ];
        for i in 0..9 {
            assert_rel_close(t2_p_s(p2[i], s2[i]).unwrap(), if97_t_ps[i], 1e-8);
        }

        let h3 = [
            2800.0, 2800.0, 4100.0, 2800.0, 3600.0, 3600.0, 2800.0, 2800.0, 3400.0,
        ];
        let s3 = [6.5, 9.5, 9.5, 6.0, 6.0, 7.0, 5.1, 5.8, 5.8];
        let if97_p_hs = [
            1.371_012_767,
            0.001_879_743_844,
            0.102_478_899_7,
            4.793_911_442,
            83.955_192_09,
            7.527_161_441,
            94.392_020_6,
            8.414_574_124,
            83.769_038_79,
        ];
        for i in 0..9 {
            assert_rel_close(p2_h_s(h3[i], s3[i]).unwrap(), if97_p_hs[i], 1e-8);
        }
    }

    #[test]
    fn region3_table33_rho_t() {
        let t = [650.0, 650.0, 750.0];
        let rho = [500.0, 200.0, 500.0];
        let if97_p = [25.583_701_8, 22.293_064_3, 78.309_563_9];
        let if97_h = [1_863.430_19, 2_375.124_01, 2_258.688_45];
        let if97_u = [1_812.262_79, 2_263.658_68, 2_102.069_32];
        let if97_s = [4.054_272_73, 4.854_387_92, 4.469_719_06];
        let if97_cp = [13.893_571_7, 44.657_934_2, 6.341_653_59];
        let if97_w = [502.005_554, 383.444_594, 760.696_041];

        for i in 0..3 {
            let p_calc = p3_rho_t(rho[i], t[i]).unwrap();
            assert_rel_close(p_calc, if97_p[i], 1e-8);
            assert_rel_close(h3_rho_t(rho[i], t[i]).unwrap(), if97_h[i], 1e-8);
            assert_rel_close(u3_rho_t(rho[i], t[i]).unwrap(), if97_u[i], 1e-8);
            assert_rel_close(s3_rho_t(rho[i], t[i]).unwrap(), if97_s[i], 1e-8);
            assert_rel_close(cp3_rho_t(rho[i], t[i]).unwrap(), if97_cp[i], 1e-8);
            assert_rel_close(w3_rho_t(rho[i], t[i]).unwrap(), if97_w[i], 1e-8);

            let rho2 = rho3_p_t(if97_p[i], t[i]).unwrap_or_else(|e| {
                panic!(
                    "rho3_p_t failed: i={i} p_mpa={} t_k={} err={e:?}",
                    if97_p[i], t[i]
                )
            });
            assert_rel_close(rho2, rho[i], 2e-8);
        }
    }

    #[test]
    fn region3_backward_t_p_h() {
        let p = [
            200.0 / 10.0,
            500.0 / 10.0,
            1000.0 / 10.0,
            200.0 / 10.0,
            500.0 / 10.0,
            1000.0 / 10.0,
        ];
        let h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0];
        let if97_t = [
            629.308_389_2,
            690.571_833_8,
            733.616_301_4,
            641.841_805_3,
            735.184_861_8,
            842.046_087_6,
        ];
        for i in 0..6 {
            assert_rel_close(t3_p_h(p[i], h[i]).unwrap(), if97_t[i], 1e-8);
        }
    }

    #[test]
    fn region3_backward_v_p_h() {
        let p = [
            200.0 / 10.0,
            500.0 / 10.0,
            1000.0 / 10.0,
            200.0 / 10.0,
            500.0 / 10.0,
            1000.0 / 10.0,
        ];
        let h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0];
        let if97_v = [
            0.001_749_903_962,
            0.001_908_139_035,
            0.001_676_229_776,
            0.006_670_547_043,
            0.002_801_244_5,
            0.002_404_234_998,
        ];
        for i in 0..6 {
            assert_rel_close(v3_p_h(p[i], h[i]).unwrap(), if97_v[i], 1e-7);
        }
    }

    #[test]
    fn region3_backward_t_p_s() {
        let p = [
            200.0 / 10.0,
            500.0 / 10.0,
            1000.0 / 10.0,
            200.0 / 10.0,
            500.0 / 10.0,
            1000.0 / 10.0,
        ];
        let s = [3.7, 3.5, 4.0, 5.0, 4.5, 5.0];
        let if97_t = [
            620.884_156_3,
            618.154_902_9,
            705.688_023_7,
            640.117_644_3,
            716.368_751_7,
            847.433_282_5,
        ];
        for i in 0..6 {
            assert_rel_close(t3_p_s(p[i], s[i]).unwrap(), if97_t[i], 1e-8);
        }
    }

    #[test]
    fn region3_backward_v_p_s() {
        let p = [
            200.0 / 10.0,
            500.0 / 10.0,
            1000.0 / 10.0,
            200.0 / 10.0,
            500.0 / 10.0,
            1000.0 / 10.0,
        ];
        let s = [3.7, 3.5, 4.0, 5.0, 4.5, 5.0];
        let if97_v = [
            0.001_639_890_984,
            0.001_423_030_205,
            0.001_555_893_131,
            0.006_262_101_987,
            0.002_332_634_294,
            0.002_449_610_757,
        ];
        for i in 0..6 {
            assert_rel_close(v3_p_s(p[i], s[i]).unwrap(), if97_v[i], 1e-8);
        }
    }

    #[test]
    fn region3_backward_p_h_s() {
        let h = [1700.0, 2000.0, 2100.0, 2500.0, 2400.0, 2700.0];
        let s = [3.8, 4.2, 4.3, 5.1, 4.7, 5.0];
        let if97_p = [
            25.557_032_46,
            45.408_734_68,
            60.781_233_4,
            17.206_124_13,
            63.639_248_87,
            88.390_432_81,
        ];
        for i in 0..6 {
            assert_rel_close(p3_h_s(h[i], s[i]).unwrap(), if97_p[i], 1e-8);
        }
    }

    #[test]
    fn region5_table42() {
        let t = [1500.0, 1500.0, 2000.0];
        let p = [0.5, 8.0, 8.0];
        let if97_v = [1.384_553_54, 0.086_515_661_6, 0.115_743_146];
        let if97_h = [5_219.763_32, 5_206.096_34, 6_583.802_91];
        let if97_u = [4_527.486_54, 4_513.971_05, 5_657.857_74];
        let if97_s = [9.654_084_31, 8.365_467_24, 9.156_710_44];
        let if97_cp = [2.616_102_28, 2.644_538_66, 2.853_067_5];
        let if97_w = [917.071_933, 919.708_859, 1_054.358_06];
        for i in 0..3 {
            assert_rel_close(v5_p_t(p[i], t[i]).unwrap(), if97_v[i], 1e-8);
            assert_rel_close(h5_p_t(p[i], t[i]).unwrap(), if97_h[i], 1e-8);
            assert_rel_close(u5_p_t(p[i], t[i]).unwrap(), if97_u[i], 1e-8);
            assert_rel_close(s5_p_t(p[i], t[i]).unwrap(), if97_s[i], 1e-8);
            assert_rel_close(cp5_p_t(p[i], t[i]).unwrap(), if97_cp[i], 1e-8);
            assert_rel_close(w5_p_t(p[i], t[i]).unwrap(), if97_w[i], 1e-8);
        }
    }

    #[test]
    fn xsteam_verification_calling_functions() {
        const FUN: &[&str] = &[
            "Tsat_p", "T_ph", "T_ps", "T_hs", "psat_T", "p_hs", "hV_p", "hL_p", "hV_T", "hL_T",
            "h_pT", "h_ps", "h_px", "h_prho", "h_Tx", "vV_p", "vL_p", "vV_T", "vL_T", "v_pT",
            "v_ph", "v_ps", "rhoV_p", "rhoL_p", "rhoV_T", "rhoL_T", "rho_pT", "rho_ph", "rho_ps",
            "sV_p", "sL_p", "sV_T", "sL_T", "s_pT", "s_ph", "uV_p", "uL_p", "uV_T", "uL_T", "u_pT",
            "u_ph", "u_ps", "CpV_p", "CpL_p", "CpV_T", "CpL_T", "Cp_pT", "Cp_ph", "Cp_ps", "CvV_p",
            "CvL_p", "CvV_T", "CvL_T", "Cv_pT", "Cv_ph", "Cv_ps", "wV_p", "wL_p", "wV_T", "wL_T",
            "w_pT", "w_ph", "w_ps", "my_pT", "my_ph", "my_ps", "tcL_p", "tcV_p", "tcL_T", "tcV_T",
            "tc_pT", "tc_ph", "tc_hs", "st_T", "st_p", "x_ph", "x_ps", "vx_ph", "vx_ps",
        ];
        const IN1: &[f64] = &[
            1.0, 1.0, 1.0, 100.0, 100.0, 84.0, 1.0, 1.0, 100.0, 100.0, 1.0, 1.0, 1.0, 1.0, 100.0,
            1.0, 1.0, 100.0, 100.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0, 100.0, 1.0, 1.0, 1.0, 0.006117,
            0.0061171, 0.0001, 100.0, 1.0, 1.0, 1.0, 1.0, 100.0, 100.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            100.0, 100.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0, 100.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0,
            100.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 25.0, 25.0, 1.0, 1.0, 100.0, 100.0, 1.0,
            1.0, 1.0, 1.0, 1.0,
        ];
        const IN2: &[f64] = &[
            0.0,
            100.0,
            1.0,
            0.2,
            0.0,
            0.296,
            0.0,
            0.0,
            0.0,
            0.0,
            20.0,
            1.0,
            0.5,
            2.0,
            0.5,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            1000.0,
            5.0,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            1000.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            20.0,
            84.01181117,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            1000.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            200.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            200.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            200.0,
            1.0,
            100.0,
            100.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            25.0,
            100.0,
            0.34,
            0.0,
            0.0,
            1000.0,
            4.0,
            418.0,
            4.0,
        ];
        const CONTROL: &[f64] = &[
            99.60591861,
            23.84481908,
            73.70859421,
            13.84933511,
            1.014179779,
            2.295498269,
            2674.949641,
            417.4364858,
            2675.572029,
            419.099155,
            84.01181117,
            308.6107171,
            1546.193063,
            1082.773391,
            1547.33559210927,
            1.694022523,
            0.001043148,
            1.671860601,
            0.001043455,
            1.695959407,
            0.437925658,
            1.03463539,
            0.590310924,
            958.6368897,
            0.598135993,
            958.3542773,
            0.589636754,
            2.283492601,
            975.6236788,
            9.155465556,
            1.8359e-05,
            9.155756716,
            1.307014328,
            0.296482921,
            0.296813845,
            2505.547389,
            417.332171,
            2506.015308,
            418.9933299,
            2506.171426,
            956.2074342,
            308.5082185,
            2.075938025,
            4.216149431,
            2.077491868,
            4.216645119,
            2.074108555,
            4.17913573168802,
            4.190607038,
            1.552696979,
            3.769699683,
            1.553698696,
            3.76770022,
            1.551397249,
            4.035176364,
            3.902919468,
            472.0541571,
            1545.451948,
            472.2559492,
            1545.092249,
            472.3375235,
            1542.682475,
            1557.8585,
            1.22704e-05,
            0.000914003770302108,
            0.000384222,
            0.677593822,
            0.024753668,
            0.607458162,
            0.018326723,
            0.607509806,
            0.605710062,
            0.606283124,
            0.0589118685876641,
            0.058987784,
            0.258055424,
            0.445397961,
            0.288493093,
            0.999233827,
        ];

        assert_eq!(FUN.len(), IN1.len(), "fun/in1 length mismatch");
        assert_eq!(FUN.len(), IN2.len(), "fun/in2 length mismatch");
        assert_eq!(FUN.len(), CONTROL.len(), "fun/control length mismatch");

        for i in 0..FUN.len() {
            let actual = props(FUN[i], IN1[i], IN2[i]).unwrap_or_else(|e| {
                panic!(
                    "props failed: i={i} fun={} in1={} in2={} err={e:?}",
                    FUN[i], IN1[i], IN2[i]
                )
            });
            assert_rel_close(actual, CONTROL[i], 1e-5);
        }
    }

    #[test]
    fn props_si_pressure_unit_is_pa() {
        let t_k = 373.15;
        let p_pa = props_si("psat_t", t_k, 0.0).unwrap();
        let expected_pa = psat_mpa_from_t_k(t_k).unwrap() * 1_000_000.0;
        assert_rel_close(p_pa, expected_pa, 1e-12);

        let t1 = props_si("tsat_p", 100_000.0, 0.0).unwrap();
        let t2 = to_si_t_k_from_c(props("tsat_p", 1.0, 0.0).unwrap());
        assert_rel_close(t1, t2, 1e-12);

        let p_pa2 = props_si("p_hs", 84.0, 0.296).unwrap();
        let p_bar2 = props("p_hs", 84.0, 0.296).unwrap();
        assert_rel_close(p_pa2, to_si_p_pa_from_bar(p_bar2), 1e-12);
    }
}
