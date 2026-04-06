pub mod error;
pub use error::{Error, Result};

use std::{collections::HashMap, f64::consts::PI, ops::Range, sync::OnceLock, thread, time::Instant};

use ndarray::parallel::prelude::*;
use ndarray::{Array2, s};
use rand::prelude::*;
use ruviz::prelude::*;
use strum::EnumCount;

// use uom::si::f64::*;
// use uom::si::{length::{millimeter, centimeter},reciprocal_length::{reciprocal_millimeter, reciprocal_centimeter}};

/// An `ndarray::Array2<f64>` of a simulation accumulated `Photon`s `weight` losses, as a slice of `Voxels` along the XZ axis
pub type VoxelsSliceXZ = Array2<f64>;

/// A `u16` wrapper to unambiguously pass a `wavelength` around
#[derive(Debug, Clone, Copy)]
pub struct Wavelength(pub u16);

/// Prepares, runs and plots a Monte Carlo simulation
///
/// ### Parameters
///
/// - `k_photon` - the number of 1000s of photons simulated
/// - `wavelength` - the wavelength of the simulated photons, in nm
/// - `verbose` - verbosity of the run
///
/// ### Errors
///
/// May return `Err` from `std::thread::available_parallelism()` or the `ruviz` crate `Plot` command
pub fn run(k_photon: u16, wavelength: Wavelength, verbose: bool) -> Result<()> {
    let available_parallelism = thread::available_parallelism()?.get();

    let nb_photons = 1_000_usize * k_photon as usize;
    let chunk_size = (nb_photons / available_parallelism) + 1;
    if verbose {
        dbg!(&SKIN_DELTA_Z_IN_UM);
    }

    let start = Instant::now();

    let op_vxls = monte_carlo(nb_photons, wavelength, chunk_size, verbose);
    println!(
        "monte_carlo of {k_photon}k photon at {}nm finished on {available_parallelism} cores in {:?}",
        wavelength.0,
        start.elapsed()
    );

    if let Some(mut vxls) = op_vxls {
        plot(&mut vxls, k_photon, wavelength, verbose)?;
    } else {
        println!("wavelenght of {}nm not supported.", &wavelength.0);
    }

    Ok(())
}

/// Plots the `VoxelsSliceXZ` returned by `monte_carlo()`
///
/// ### Parameters
///
/// - `vxl_xz_slice` - the `ndarray::Array2<f64>` of a simulation accumulated `Photon`s `weight` losses, as a slice of
///   `Voxels` along the XZ axis
/// - `k_photon` - the number of 1000s of `Photon`s simulated
/// - `wavelength` - the `wavelength` of the simulated `Photon`s
/// - `verbose` - verbosity of the plotting
///
/// ### Errors
///
/// May return `Err` from the `ruviz` crate `Plot` command
pub fn plot(vxl_xz_slice: &mut VoxelsSliceXZ, k_photon: u16, wavelength: Wavelength, verbose: bool) -> Result<()> {
    let binding = vxl_xz_slice.t();
    let view = binding.slice(s![..;-1, ..]);

    let x_axis_scaling = f64::from(VXLS_X_SIZE) * MM_PER_VXL;
    let z_axis_scaling = f64::from(VXLS_Z_SIZE) * MM_PER_VXL;

    let start = Instant::now();

    Plot::new()
        .heatmap(
            &view,
            Some(
                HeatmapConfig::new()
                     .value_scale(AxisScale::SymLog { linthresh: (0.000_001) })
                     .colorbar(true)
                     .colorbar_label("Absorbed Energy")
            ),
        )
        .xlim(0., x_axis_scaling)
        .xlabel("Position in x(mm)")
        .ylim(z_axis_scaling, 0.)
        .ylabel("Depth (mm)")
        .title(format!("{k_photon}k photons at {}nm", wavelength.0))
        .save(format!("heatmap-{}nm{k_photon}kphoton.png", wavelength.0))?;

    if verbose {
        println!("Plot finished in {:?}", start.elapsed());
    }

    Ok(())
}

/// Number of mm in a cm
const MM_PER_CM: f64 = 10.;
/// Number of µm in a single Voxel
const UM_PER_VXL: u16 = 10;
/// Number of mm in a single Voxel
const MM_PER_VXL: f64 = UM_PER_VXL as f64 / 1_000.;

/// The different kinds of skin layers
#[derive(Debug, EnumCount, Clone, Copy, PartialEq)]
enum SkinLayerKind {
    StratumCorneum,
    Epiderme,
    PapileDermique,
    DermeSuperieur,
    DermeReticulaire,
    DermeProfond,
}
use SkinLayerKind::{DermeProfond, DermeReticulaire, DermeSuperieur, Epiderme, PapileDermique, StratumCorneum};

/// The static model for a `Skin` layer
#[derive(Debug, Clone, Copy)]
struct SkinLayer {
    /// The kind of skin layer
    kind: SkinLayerKind,
    /// The depth of a skin layer in µm
    dz_in_um: u16,
    /// The indice de refraction of a skin layer  
    indice_refrac: f64,
    /// The `v_b` of a skin layer  
    v_b: f64,
    /// The `v_w` of a skin layer  
    v_w: f64,
}

/// Complete static model of the `SkinLayer`s of a simulation
#[rustfmt::skip]
const SKIN_LAYERS: [SkinLayer; SkinLayerKind::COUNT] = [
    SkinLayer { dz_in_um:  20, indice_refrac: 1.42, v_b: 0.00, v_w: 0.05, kind: StratumCorneum },
    SkinLayer { dz_in_um: 250, indice_refrac: 1.42, v_b: 0.00, v_w: 0.20, kind: Epiderme },
    SkinLayer { dz_in_um: 100, indice_refrac: 1.39, v_b: 0.04, v_w: 0.50, kind: PapileDermique },
    SkinLayer { dz_in_um:  80, indice_refrac: 1.39, v_b: 0.30, v_w: 0.60, kind: DermeSuperieur },
    SkinLayer { dz_in_um: 200, indice_refrac: 1.39, v_b: 0.04, v_w: 0.70, kind: DermeReticulaire },
    SkinLayer { dz_in_um: 300, indice_refrac: 1.39, v_b: 0.01, v_w: 0.70, kind: DermeProfond },
];

/// Sum of the `.dz_un_um` of all `SkinLayer`s in `SKIN_LAYERS`
const SKIN_DELTA_Z_IN_UM: u16 = SKIN_LAYERS[StratumCorneum as usize].dz_in_um
    + SKIN_LAYERS[Epiderme as usize].dz_in_um
    + SKIN_LAYERS[PapileDermique as usize].dz_in_um
    + SKIN_LAYERS[DermeSuperieur as usize].dz_in_um
    + SKIN_LAYERS[DermeReticulaire as usize].dz_in_um
    + SKIN_LAYERS[DermeProfond as usize].dz_in_um;

/// Number of `Voxels` in the `x` dimension
const VXLS_X_SIZE: i32 = 800; //  8.0 mm
/// Number of `Voxels` in the `y` dimension
const VXLS_Y_SIZE: i32 = 800; //  8.0 mm
/// Number of `Voxels` in the `z` dimension
const VXLS_Z_SIZE: i32 = (SKIN_DELTA_Z_IN_UM / UM_PER_VXL) as i32;
/// `(VXLS_X_SIZE as usize, VXLS_Z_SIZE as usize)`
const VXL_XZ_DIMS: (usize, usize) = (VXLS_X_SIZE as usize, VXLS_Z_SIZE as usize);

/// A `f64` wrapper to unambiguously pass a `delta_w_coeff` around
#[derive(Debug, Clone, Copy, Default)]
struct DeltaWCoeff(f64);

/// The complete skin layer model used in the simulation, with wavelength-related computed members
#[derive(Debug)]
struct SkinLayerAtWl {
    /// Reference to the static skin layer model it extends
    layer: &'static SkinLayer,

    /// Ratio of this layer `indice_refrac` to that of the skin layer above in the model
    pub indice_refrac_ratio: f64,
    /// Range of voxel `z` coordinates contained in this skin ayer, considering all those above it
    pub vxl_z_range: Range<i32>,

    /// Computed `mu_t` of this skin layer, in `reciprocal_vxl` or vxl^-1
    pub mu_t: f64, //
    /// Computed `mu_a/mu_t` of this skin layer
    pub mu_a_on_mu_t: DeltaWCoeff,
}

impl SkinLayerAtWl {
    /// Returns this skin layer `SkinLayerKind`
    const fn kind(&self) -> SkinLayerKind {
        self.layer.kind
    }

    /// Returns this skin layer `.indice_refrac`
    const fn indice_refrac(&self) -> f64 {
        self.layer.indice_refrac
    }

    /// Recomputes the `mu_t` and `mu_a_on_mu_t` member of the skin layer according to some `wavelength` parameters
    fn model_at(&mut self, wavelength_f: f64, val_hb_0: f64, val_hb_1: f64, val_water: f64) {
        //  "J'ai assumé que l'output était en cm^-1, donc divisé par 10 pour mm^-1"
        const MAGIG_NUM_0: f64 = 2.303 * 150. / 64_500. / MM_PER_CM;
        const V_MEL: f64 = 0.1;
        /// Saturation en oxygène artériel de 98%
        const SA_O2: f64 = 0.98;
        /// Saturation en oxygène veineux est 10% moindre que l'arteriel
        const SV_O2: f64 = 0.9 * SA_O2;
        const MU_S: f64 = 20.1 * MM_PER_VXL;

        //  "J'ai assumé que l'output était en cm^-1, donc divisé par 10 pour mm^-1"
        let mu_a_baseline = 7.84e-7 * wavelength_f.powf(-3.255) / MM_PER_CM; //  reciprocal_mm

        let mu_absorption_oxyhemoglobine = val_hb_0 * MAGIG_NUM_0; // reciprocal_mm
        let mu_a_hemoglobine = val_hb_1 * MAGIG_NUM_0; //  reciprocal_mm

        let mu_w = val_water / MM_PER_CM; //  reciprocal_mm

        //  Coeff artere & veine
        let (mu_art, mu_vei) = (
            SA_O2.mul_add(mu_absorption_oxyhemoglobine, (1. - SA_O2) * mu_a_hemoglobine),
            SV_O2.mul_add(mu_absorption_oxyhemoglobine, (1. - SV_O2) * mu_a_hemoglobine),
        );

        let SkinLayer { kind, v_w, v_b, .. } = *self.layer;

        let mu_a = MM_PER_VXL
            * match kind {
                StratumCorneum | Epiderme => (1. - (V_MEL + v_w)).mul_add(
                    mu_a_baseline,
                    V_MEL.mul_add(6.6e+10 * wavelength_f.powf(-3.33) / MM_PER_CM, v_w * mu_w),
                ),
                _ => {
                    //  Tout sous l'Epiderme
                    let v_a = v_b / 2.; //  split 50-50 with v_v
                    let v_v = v_a;
                    (1. - (v_a + v_v + v_w)).mul_add(mu_a_baseline, v_w.mul_add(mu_w, v_a * mu_art + v_v * mu_vei))
                }
            };
        self.mu_t = mu_a + MU_S;
        self.mu_a_on_mu_t = DeltaWCoeff(mu_a / self.mu_t);
    }
}

/// Lowest photon `Wavelength` allowed by this simulator.
const LOWEST_WAVELENGTH: u16 = 600; // nanometer
/// Number of even `Wavelength` supported starting at `LOWEST_WAVELENGTH`
const HB_SIZE: u16 = 201;
///  Source Scott Prahl, (lambda en nm, le reste en cm^-1/M)
const HB_ARRAY: [(f64, f64); HB_SIZE as usize] = [
    (3_200.0, 14_677.2),
    (2_664.0, 13_622.4),
    (2_128.0, 12_567.6),
    (1_789.2, 11_513.2),
    (1_647.6, 10_477.6),
    (1_506.0, 9_443.6),
    (1_364.4, 8_591.2),
    (1_222.8, 7_762.0),
    (1_110.0, 7_344.8),
    (1_026.0, 6_927.2),
    (942.0, 6_509.6),
    (858.0, 6_193.2),
    (774.0, 5_906.8),
    (707.6, 5_620.0),
    (658.8, 5_366.8),
    (610.0, 5_148.8),
    (561.2, 4_930.8),
    (512.4, 4_730.8),
    (478.8, 4_602.4),
    (460.4, 4_473.6),
    (442.0, 4_345.2),
    (423.6, 4_216.8),
    (405.2, 4_088.4),
    (390.4, 3_965.08),
    (379.2, 3_857.6),
    (368.0, 3_750.12),
    (356.8, 3_642.64),
    (345.6, 3_535.16),
    (335.2, 3_427.68),
    (325.6, 3_320.2),
    (319.6, 3_226.56),
    (314.0, 3_140.28),
    (308.4, 3_053.96),
    (302.8, 2_967.68),
    (298.0, 2_881.4),
    (294.0, 2_795.12),
    (290.0, 2_708.84),
    (285.6, 2_627.64),
    (282.0, 2_554.4),
    (279.2, 2_481.16),
    (277.6, 2_407.92),
    (276.0, 2_334.68),
    (274.4, 2_261.48),
    (272.8, 2_188.24),
    (274.4, 2_115.0),
    (276.0, 2_051.96),
    (277.6, 2_000.48),
    (279.2, 1_949.04),
    (282.0, 1_897.56),
    (286.0, 1_846.08),
    (290.0, 1_794.28),
    (294.0, 1_741.0),
    (298.0, 1_687.76),
    (302.8, 1_634.48),
    (308.4, 1_583.52),
    (314.0, 1_540.48),
    (319.6, 1_497.4),
    (325.2, 1_454.36),
    (332.0, 1_411.32),
    (340.0, 1_368.28),
    (348.0, 1_325.88),
    (356.0, 1_285.16),
    (364.0, 1_244.44),
    (372.4, 1_203.68),
    (381.2, 1_152.8),
    (390.0, 1_102.2),
    (398.8, 1_102.2),
    (407.6, 1_102.2),
    (418.8, 1_101.76),
    (432.4, 1_100.48),
    (446.0, 1_115.88),
    (459.6, 1_161.64),
    (473.2, 1_207.4),
    (487.6, 1_266.04),
    (502.8, 1_333.24),
    (518.0, 1_405.24),
    (533.2, 1_515.32),
    (548.4, 1_541.76),
    (562.0, 1_560.48),
    (574.0, 1_560.48),
    (586.0, 1_548.52),
    (598.0, 1_508.44),
    (610.0, 1_459.56),
    (622.8, 1_410.52),
    (636.4, 1_361.32),
    (650.0, 1_311.88),
    (663.6, 1_262.44),
    (677.2, 1_213.0),
    (689.2, 1_163.56),
    (699.6, 1_114.8),
    (710.0, 1_075.44),
    (720.4, 1_036.08),
    (730.8, 996.72),
    (740.0, 957.36),
    (748.0, 921.8),
    (756.0, 890.8),
    (764.0, 859.8),
    (772.0, 828.8),
    (786.4, 802.96),
    (807.2, 782.36),
    (816.0, 761.72),
    (828.0, 743.84),
    (836.0, 737.08),
    (844.0, 730.28),
    (856.0, 723.52),
    (864.0, 717.08),
    (872.0, 711.84),
    (880.0, 706.6),
    (887.2, 701.32),
    (901.6, 696.08),
    (916.0, 693.76),
    (930.4, 693.6),
    (944.8, 693.48),
    (956.4, 693.32),
    (965.2, 693.2),
    (974.0, 693.04),
    (982.8, 692.92),
    (991.6, 692.76),
    (1_001.2, 692.64),
    (1_011.6, 692.48),
    (1_022.0, 692.36),
    (1_032.4, 692.2),
    (1_042.8, 691.96),
    (1_050.0, 691.76),
    (1_054.0, 691.52),
    (1_058.0, 691.32),
    (1_062.0, 691.08),
    (1_066.0, 690.88),
    (1_072.8, 690.64),
    (1_082.4, 692.44),
    (1_092.0, 694.32),
    (1_101.6, 696.2),
    (1_111.2, 698.04),
    (1_118.4, 699.92),
    (1_123.2, 701.8),
    (1_128.0, 705.84),
    (1_132.8, 709.96),
    (1_137.6, 714.08),
    (1_142.8, 718.2),
    (1_148.4, 722.32),
    (1_154.0, 726.44),
    (1_159.6, 729.84),
    (1_165.2, 733.2),
    (1_170.0, 736.6),
    (1_174.0, 739.96),
    (1_178.0, 743.6),
    (1_182.0, 747.24),
    (1_186.0, 750.88),
    (1_190.0, 754.52),
    (1_194.0, 758.16),
    (1_198.0, 761.84),
    (1_202.0, 765.04),
    (1_206.0, 767.44),
    (1_209.2, 769.8),
    (1_211.6, 772.16),
    (1_214.0, 774.56),
    (1_216.4, 776.92),
    (1_218.8, 778.4),
    (1_220.8, 778.04),
    (1_222.4, 777.72),
    (1_224.0, 777.36),
    (1_225.6, 777.04),
    (1_227.2, 776.64),
    (1_226.8, 772.36),
    (1_224.4, 768.08),
    (1_222.0, 763.84),
    (1_219.6, 752.28),
    (1_217.2, 737.56),
    (1_215.6, 722.88),
    (1_214.8, 708.16),
    (1_214.0, 693.44),
    (1_213.2, 678.72),
    (1_212.4, 660.52),
    (1_210.4, 641.08),
    (1_207.2, 621.64),
    (1_204.0, 602.24),
    (1_200.8, 583.4),
    (1_197.6, 568.92),
    (1_194.0, 554.48),
    (1_190.0, 540.04),
    (1_186.0, 525.56),
    (1_182.0, 511.12),
    (1_178.0, 495.36),
    (1_173.2, 473.32),
    (1_167.6, 451.32),
    (1_162.0, 429.32),
    (1_156.4, 415.28),
    (1_150.8, 402.28),
    (1_144.0, 389.288),
    (1_136.0, 374.944),
    (1_128.0, 359.656),
    (1_120.0, 344.372),
    (1_112.0, 329.084),
    (1_102.4, 313.796),
    (1_091.2, 298.508),
    (1_080.0, 283.22),
    (1_068.8, 267.932),
    (1_057.6, 252.648),
    (1_046.4, 237.36),
    (1_035.2, 222.072),
    (1_024.0, 206.784),
];

/// Returns `Some` `hb` specs of a `Wavelength` if it's known, `None` otherwise
fn try_get_hb_for(wavelength: Wavelength) -> Option<(f64, f64)> {
    if (LOWEST_WAVELENGTH..LOWEST_WAVELENGTH + HB_SIZE * 2).contains(&wavelength.0) && wavelength.0.is_multiple_of(2) {
        Some(HB_ARRAY[usize::from((wavelength.0 - LOWEST_WAVELENGTH) / 2)])
    } else {
        None
    }
}
/// Source Segelstein, cm^-1
static WATER_HASHMAP: OnceLock<HashMap<u16, f64>> = OnceLock::new();

/// Returns `Some` `water` specs of a `Wavelength` if it's known, `None` otherwise
#[allow(clippy::too_many_lines)]
fn try_get_water_for(wavelength: Wavelength) -> Option<f64> {
    WATER_HASHMAP
        .get_or_init(|| {
            [
                (600, 0.002_018_4),
                (605, 0.002_350_1),
                (610, 0.002_552_4),
                (615, 0.002_716_7),
                (619, 0.002_838_3),
                (625, 0.002_958_7),
                (629, 0.002_998_4),
                (635, 0.003_069_9),
                (640, 0.003_084_1),
                (646, 0.003_125_5),
                (650, 0.003_235_8),
                (655, 0.003_411_3),
                (661, 0.003_689_8),
                (665, 0.003_836_2),
                (670, 0.003_935_5),
                (675, 0.004_055_9),
                (681, 0.004_245_4),
                (685, 0.004_529_8),
                (690, 0.004_830_3),
                (695, 0.005_357_4),
                (700, 0.006_012_0),
                (705, 0.007_311_2),
                (710, 0.008_851_0),
                (715, 0.010_544),
                (719, 0.012_736),
                (724, 0.015_850),
                (729, 0.019_810),
                (735, 0.023_063),
                (740, 0.024_773),
                (745, 0.025_818),
                (750, 0.026_125),
                (755, 0.026_294),
                (760, 0.026_114),
                (766, 0.025_770),
                (769, 0.024_950),
                (775, 0.023_981),
                (780, 0.022_706),
                (785, 0.021_429),
                (791, 0.020_374),
                (794, 0.019_902),
                (800, 0.019_640),
                (805, 0.019_815),
                (809, 0.020_657),
                (815, 0.022_335),
                (820, 0.024_829),
                (824, 0.027_737),
                (830, 0.030_905),
                (836, 0.033_732),
                (840, 0.036_808),
                (845, 0.039_990),
                (849, 0.043_343),
                (855, 0.046_336),
                (859, 0.048_978),
                (865, 0.051_515),
                (871, 0.054_074),
                (875, 0.056_111),
                (879, 0.057_942),
                (885, 0.060_113),
                (889, 0.062_224),
                (895, 0.064_867),
                (900, 0.067_924),
                (906, 0.071_455),
                (910, 0.078_707),
                (914, 0.092_052),
                (920, 0.113_38),
                (925, 0.144_05),
                (929, 0.185_05),
                (935, 0.237_92),
                (940, 0.290_05),
                (944, 0.340_35),
                (951, 0.387_59),
                (955, 0.419_76),
                (959, 0.439_84),
                (966, 0.450_57),
                (971, 0.453_45),
                (975, 0.448_52),
                (980, 0.438_51),
                (984, 0.426_03),
                (991, 0.412_58),
                (995, 0.395_27),
                (1_000, 0.376_99),
                (1_009, 0.334_77),
                (1_021, 0.289_48),
                (1_030, 0.244_13),
                (1_040, 0.204_20),
                (1_050, 0.169_83),
                (1_059, 0.154_14),
                (1_069, 0.148_00),
                (1_079, 0.154_78),
                (1_089, 0.172_97),
                (1_099, 0.195_30),
                (1_109, 0.230_93),
                (1_119, 0.295_12),
                (1_130, 0.430_26),
                (1_140, 0.655_99),
                (1_151, 1.016_0),
                (1_159, 1.159_1),
                (1_169, 1.204_0),
            ]
            .into_iter()
            .collect()
        })
        .get(&wavelength.0)
        .copied()
}

/// Skin model used in the Monte Carlo simulation
#[derive(Debug)]
struct Skin {
    layers: [SkinLayerAtWl; SkinLayerKind::COUNT],
    wavelength: Wavelength,
}

/// Converts µm to voxel discrete unit.
const fn um_to_vxl(um: u16) -> i32 {
    (um / UM_PER_VXL) as i32
}

impl Default for Skin {
    fn default() -> Self {
        const INDICE_REFRAC_AIR: f64 = 1.000_293; // à T.P.N. Unused, but for isometric solution.

        let mut dz_in_vxl_acc = 0;
        let mut indice_refrac_denominateur = INDICE_REFRAC_AIR;

        Self {
            layers: SKIN_LAYERS.iter().map(|skin_layer| {
                let SkinLayer { dz_in_um, indice_refrac, .. } = skin_layer;

                let dz_in_vxl_start = dz_in_vxl_acc;
                dz_in_vxl_acc += um_to_vxl(*dz_in_um);

                let indice_refrac_numerateur = indice_refrac_denominateur;  // precedent
                indice_refrac_denominateur = *indice_refrac;

                SkinLayerAtWl {
                    layer: skin_layer,

                    indice_refrac_ratio: indice_refrac_numerateur / indice_refrac_denominateur,
                    vxl_z_range: dz_in_vxl_start..dz_in_vxl_acc,

                    mu_t: 0.,
                    mu_a_on_mu_t: DeltaWCoeff::default(),
                }
            }).collect::<Vec<_>>().try_into().expect("SKIN_LAYERS [SkinLayer; SkinLayerKind::COUNT] to SkinLayerAtWl[SkinLayerAtWl; SkinLayerKind::COUNT]"),
            wavelength: Wavelength(LOWEST_WAVELENGTH),
        }
    }
}

impl Skin {
    /// Returns `Some` version of this `Skin` model if it can be built for a given `Wavelength`, or `None` otherwise
    fn try_model_at(&mut self, wavelength: Wavelength) -> Option<&mut Self> {
        if let Some((val_hb_0, val_hb_1)) = try_get_hb_for(wavelength)
            && let Some(val_water) = try_get_water_for(wavelength)
        {
            let wavelength_f = f64::from(wavelength.0);

            for skin_layer_at_wl in &mut self.layers {
                skin_layer_at_wl.model_at(wavelength_f, val_hb_0, val_hb_1, val_water);
            }
            self.wavelength = wavelength;
            Some(self)
        } else {
            None
        }
    }

    /// Returns the `SkinLayerAtWl` for a `SkinLayerKind` passed in parameter
    const fn skin_layer(&self, skin_layer_kind: SkinLayerKind) -> &SkinLayerAtWl {
        &self.layers[skin_layer_kind as usize]
    }

    /// Randomly generates a new `path_seg.len_in_vxl` according to a `SkinLayerKind` passed in parameter
    fn new_path_seg_len_in_vxl(&self, src_skin_layer_kind: SkinLayerKind, rng: &mut ThreadRng) -> f64 {
        //  parcours du groupe de photon pour cette itération
        let mut xsi_1: f64 = rng.random::<f64>();

        //  Exclude 0 from the interval because ln(0) = ∞, or NaN.
        #[allow(clippy::while_float)]
        while 0. == xsi_1 {
            xsi_1 = rng.random::<f64>();
        }
        //  -xsi_1.ln() runs from 0. to very large (f64::MAX ?).
        -xsi_1.ln() / self.skin_layer(src_skin_layer_kind).mu_t //  in vxl
    }

    /// Returns `Some` `SkinLayerAtWl` if the photon depth passed in parameter is within the `Voxels`, `None` otherwise
    fn try_skin_layer(&self, photon_pos_z: i32) -> Option<&SkinLayerAtWl> {
        self.layers
            .iter()
            .find(|skin_layer| skin_layer.vxl_z_range.contains(&photon_pos_z))
    }

    /// Simulate the absorption of a `Photon`, moving it into the `Voxels`
    ///
    /// Returns `Some` pair of new `PhotonPathSeg` direction and `delta_w_coeff` to be applied,
    /// or `None` if the `Photon` displacement moves it out of the `Voxels`
    fn absorption(&self, photon: &mut Photon, rng: &mut ThreadRng) -> Option<(UnitVec, DeltaWCoeff)> {
        // let Photon { pos, path_seg, .. } = photon;
        let SkinLayer {
            kind: src_skin_layer_kind,
            indice_refrac: n1,
            ..
        } = *photon.skin_layer(self);
        let path_seg = photon.path_seg;
        let src_pos_z = photon.pos.z;

        // With `len_in_vxl: f64` possibly extremely large, `photon.path_seg.dz()` conversion
        // to `i32` may resut in `i32::MIN` or `i32:MAX`. We insure the addition here doesn't
        // panic (or wrap-around in `release` mode) by using `i32::saturating_add_signed()`.
        if let Some(dst_skin_layer) = self.try_skin_layer(src_pos_z.saturating_add(photon.path_seg.dz())) {
            // Some(dst_skin_layer) means within voxels in z

            let mu_a_on_mu_t = dst_skin_layer.mu_a_on_mu_t;
            let dst_skin_layer_kind = dst_skin_layer.kind();

            //  Trans skin layers
            if src_skin_layer_kind == dst_skin_layer_kind {
                //  photon stays in the same skin layer
                //  returns photon.pos.is_within_xy_of_vxls()
                if photon.move_of(path_seg.len_in_vxl, 0., dst_skin_layer_kind) {
                    // photon is staying within voxels in z
                    Some((photon.path_seg.dir.iso_refrac(rng), mu_a_on_mu_t))
                } else {
                    None
                } //  not within voxels.
            } else {
                //  Stop at skin layer transition for this iteration.

                let moving_len = path_seg.len_in_vxl_from(dst_skin_layer.vxl_z_range.start - src_pos_z);
                //  returns photon.pos.is_within_xy_of_vxls()
                if photon.move_of(moving_len, path_seg.len_in_vxl - moving_len, dst_skin_layer_kind) {
                    // skin layer transition is always within voxels in z
                    // skin layers with different indice_refrac: the path_seg.dir will change
                    #[allow(clippy::float_cmp)] //  as we do want to filter the exact know ratios = 1.0
                    if 1.0 == dst_skin_layer.indice_refrac_ratio {
                        Some((photon.path_seg.dir.iso_refrac(rng), mu_a_on_mu_t))
                    } else {
                        Some(photon.path_seg.dir.aniso_refrac(
                            n1,
                            dst_skin_layer.indice_refrac(),
                            dst_skin_layer.indice_refrac_ratio,
                            rng,
                        ))
                    }
                } else {
                    None
                } //  not within voxels.
            }
        }
        // ! Some(dst_skin_layer) means not within voxels in z
        else {
            None
        }
    }
}

/// A trio of `i32` Ranges along the `x`, `y` and `z` axis allowing to find if a `Photon` is located within the `Voxels`
struct VxlRanges {
    x: Range<i32>,
    y: Range<i32>,
    z: Range<i32>,
}

impl VxlRanges {
    /// Instantiate a `VxlRanges` from its `x`. `y`, and `z` dimensions
    pub const fn new(x_size: i32, y_size: i32, z_size: i32) -> Self {
        Self {
            x: 0..x_size,
            y: 0..y_size,
            z: 0..z_size,
        }
    }
    /// Return true if a `Photon` `x` position is contained within the `Voxels`
    fn x_contains(&self, x: i32) -> bool {
        self.x.contains(&x)
    }
    /// Return true if a `Photon` `y` position is contained within the `Voxels`
    fn y_contains(&self, y: i32) -> bool {
        self.y.contains(&y)
    }
    /// Return true if a `Photon` `z` position is contained within the `Voxels`
    fn z_contains(&self, z: i32) -> bool {
        self.z.contains(&z)
    }
}

/// The single instance of `VxlRanges` used in this simulation
const VXL_RANGES: VxlRanges = VxlRanges::new(VXLS_X_SIZE, VXLS_Y_SIZE, VXLS_Z_SIZE);

/// `A i32`-based 3D position of a `Photon`, in or out of the `Voxels`
#[derive(Debug, Clone, Copy)]
struct VoxelPos {
    ///  Position of a `Photon` along the `x` axis, in or out of the `Voxels`
    x: i32,
    ///  Position of a `Photon` along the `y` axis, in or out of the `Voxels`
    y: i32,
    ///  Position of a `Photon` along the `z` axis, in or out of the `Voxels`
    z: i32,
}

impl VoxelPos {
    ///  Returns `x as usize`.
    ///
    /// `x` is validated as never negative by `skin.try_skin_layer()` at all times, so usize conversion is *always* valid.
    const fn x(&self) -> usize {
        #![allow(clippy::cast_sign_loss)]
        self.x as usize
    }

    ///  Returns `z as usize`.
    ///
    ///  `z` is validated as never negative  by `is_within_xy_of_vxls()` at all times, so usize conversion is *always* valid.
    const fn z(&self) -> usize {
        #![allow(clippy::cast_sign_loss)]
        self.z as usize
    }

    /// Test that this `VoxelPos` is within the Voxels x, y and z ranges.
    fn is_within_vxls(&self) -> bool {
        VXL_RANGES.x_contains(self.x) && VXL_RANGES.y_contains(self.y) && VXL_RANGES.z_contains(self.z)
    }

    /// Test that this `VoxelPos` is within the Voxels x and y ranges.
    fn is_within_xy_of_vxls(&self) -> bool {
        VXL_RANGES.x_contains(self.x) && VXL_RANGES.y_contains(self.y)
    }

    /// Adds the `delta_pos` to `self` and returns .`is_within_xy_of_vxls()`.
    ///
    /// Since it is likely that `dx`, `dy` and or `dz` ends up being `i32::MIN` or `i32:MAX`,
    /// on overflow by `f64` conversion to `i32`, it's important to prevent panic (or worst: silent
    /// wrap-around in release mode). So, we cap the `x`, `y`, `z` values to `i32::MIN` and `i32:MAX`
    /// by using `i32::saturating_add_signed()`.
    ///
    /// `i32::MIN` and `i32:MAX` values will then definitely fail the final `is_within_xy_of_vxls()`.
    fn move_of(&mut self, delta_pos: DeltaVoxelPos) -> bool {
        self.x = self.x.saturating_add(delta_pos.dx);
        self.y = self.y.saturating_add(delta_pos.dy);
        self.z = self.z.saturating_add(delta_pos.dz);
        self.is_within_xy_of_vxls()
    }
}

/// A `i32`-based 3D distance between two positions of a `Photon`, in or out of the `Voxels`
#[derive(Debug, Clone, Copy)]
struct DeltaVoxelPos {
    ///  Distance between two positions of a `Photon` along the `x` axis, in or out of the `Voxels`
    dx: i32,
    ///  Distance between two positions of a `Photon` along the `y` axis, in or out of the `Voxels`
    dy: i32,
    ///  Distance between two positions of a `Photon` along the `z` axis, in or out of the `Voxels`
    dz: i32,
}

/// A `f64`-based 3D unitary vector used as the direction of a `PhotonPathSeg` random-length displacement in the simulation
#[derive(Debug, Clone, Copy)]
struct UnitVec {
    /// Direction, along the `x` axis, of a `PhotonPathSeg` random-length displacement in the simulation
    pub ux: f64,
    /// Direction, along the `y` axis, of a `PhotonPathSeg` random-length displacement in the simulation
    pub uy: f64,
    /// Direction, along the `z` axis, of a `PhotonPathSeg` random-length displacement in the simulation
    pub uz: f64,
}

impl Default for UnitVec {
    /// The default `{ ux: 0., uy: 0., uz: 1. }` direction that all new simulatad photon start with
    fn default() -> Self {
        Self { ux: 0., uy: 0., uz: 1. }
    }
}

/// Returns sin if parameter is cos and cos if parameter is sin
fn co_sin(sin_or_cos: f64) -> f64 {
    sin_or_cos.mul_add(-sin_or_cos, 1.).sqrt()
}

impl UnitVec {
    /// Instantiates a new `UnitVec`    
    const fn new(ux: f64, uy: f64, uz: f64) -> Self {
        Self { ux, uy, uz }
    }

    /// Returns `f64` based `[ux, uy, uz] * len` converted to `i32` based `[dx, dy, dz]`.
    /// Since we *know* that 'len' can be extremely big (ln(0) = ∞, we prevent `0.` but not
    /// `f64::EPSILON`) it is entirely likely that some returned `dx`, `dy` and or `dz` ends
    /// up being `i32::MIN` or `i32:MAX` on overflow by `f64` conversion to `i32`.
    fn delta_pos(&self, len: f64) -> DeltaVoxelPos {
        //  As `i32` returns `i32:MIN` or `i32::MAX` on overflow by `f64`conversion, which
        // //  will definitely ultimately result in `pos.is_within_xy_of_vxls()` failting.
        #[allow(clippy::cast_possible_truncation)]
        DeltaVoxelPos {
            dx: (self.ux * len) as i32,
            dy: (self.uy * len) as i32,
            dz: (self.uz * len) as i32,
        }
    }

    /// Changes the direction of the `PhotonPathSeg` so that it's refelected along the `z` axis.
    fn uz_reflected(&self) -> Self {
        //  Only .uz is reflected the other unitary element ux and uy stay the same.
        Self {
            ux: self.ux,
            uy: self.uy,
            uz: -self.uz,
        }
    }

    /// Randomly changes the direction of the `PhotonPathSeg` to reflect displacement through matter with uniform `indice_refrac`
    fn iso_refrac(&self, rng: &mut ThreadRng) -> Self {
        /// Coefficient d'anisotropie
        const G: f64 = 0.9;
        const G_SQUARE: f64 = G * G;
        const ONE_PLUS_G_SQUARE: f64 = 1. + G_SQUARE;
        const ONE_MINUS_G_SQUARE: f64 = 1. - G_SQUARE;
        const ONE_MINUS_G: f64 = 1. - G;
        const TWO_G: f64 = 2. * G;
        const ONE_ON_TWO_G: f64 = 1. / (2. * G);

        // xsi_2
        let phi = 2.0 * PI * rng.random::<f64>();
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();

        //  xsi_3
        let mut cos_theta = ONE_ON_TWO_G
            * (ONE_MINUS_G_SQUARE / TWO_G.mul_add(rng.random::<f64>(), ONE_MINUS_G)).mul_add(
                -(ONE_MINUS_G_SQUARE / TWO_G.mul_add(rng.random::<f64>(), ONE_MINUS_G)),
                ONE_PLUS_G_SQUARE,
            );
        //  Ici aussi, pour éviter les erreurs arithmétiques flottantes cos doit être entre (-1, 1)
        cos_theta = cos_theta.clamp(-1., 1.);
        let sin_theta = co_sin(cos_theta);

        //  We know this is likely because of `common `default() intitialisation case.
        #[allow(clippy::float_cmp)]
        if 1. == self.uz.abs() {
            Self::new(
                sin_theta * sin_phi, //  self.ux
                sin_theta * cos_phi, //  self.uy
                cos_theta,           //  self.uz
            )
        } else {
            let not_uz = co_sin(self.uz);
            let uz_cos_phi = self.uz * cos_phi;
            Self::new(
                cos_theta.mul_add(self.ux, sin_theta * self.ux.mul_add(uz_cos_phi, -(self.uy * sin_phi)) / not_uz),
                cos_theta.mul_add(self.uy, sin_theta * self.uy.mul_add(uz_cos_phi, -(self.ux * sin_phi)) / not_uz),
                self.uz.mul_add(cos_theta, -(not_uz * sin_theta * cos_phi)),
            )
        }
    }

    /// Changes the direction of the `PhotonPathSeg` to randomly simulate either reflexion or transmission
    /// at matter interface with different `indice_refrac`.
    fn aniso_refrac(&self, n1: f64, n2: f64, n1_on_n2: f64, rng: &mut ThreadRng) -> (Self, DeltaWCoeff) {
        let cos_theta_i = self.uz.abs().clamp(0., 1.); //  cos_theta_i in [0., 1.] range
        let sin_theta_i = co_sin(cos_theta_i); //  sin_theta_i in [0., 1.] range too.

        let sin_theta_t = (n1_on_n2 * sin_theta_i).clamp(0.0, 1.0); //  Loi de Snell-Descartes

        //  Reflexion / Transmission computing
        //  We know this is likely because of known n1_on_n2 values larger than 1 and the clamping just above.
        #[allow(clippy::float_cmp)]
        if 1.0 == sin_theta_t {
            //  Complete internal reflexion
            //  rng.random::<f64>() generates a random number over 0.0..1.0, so always < 1.0.
            //  Therefore, if 1.0 == sin_theta_t, r = sin_theta_t is always > rng.random::<f64>()
            //  without even having to call the rng.
            //  Also r_square is 1.0, and delta_w_coeff = (1. - r_square) = 0.
            (self.uz_reflected(), DeltaWCoeff::default())
        } else {
            let cos_theta_t = co_sin(sin_theta_t);
            let r = n1.mul_add(cos_theta_i, -(n2 * cos_theta_t)) / n1.mul_add(cos_theta_i, n2 * cos_theta_t);
            let r_square = r.powi(2);

            if rng.random::<f64>() < r {
                //  xsi_5 <  r  : Reflexion
                (self.uz_reflected(), DeltaWCoeff(1. - r_square))
            } else {
                //  xsi_5 >= r  : Transmission
                (
                    Self::new(
                        self.ux * n1_on_n2, //  self.ux
                        self.uy * n1_on_n2, //  self.uy
                        cos_theta_t,        //  self.uz
                    ),
                    DeltaWCoeff(r_square),
                )
            }
        }
    }
}

/// A random-length displacement applied to a `Photon` for the few absorption cycles it takes to reach `0.`
#[derive(Debug, Clone, Copy, Default)]
struct PhotonPathSeg {
    /// Loss of `weight` accumulated over the simulation cycles, until the `len_in_vxl` reaches 0.
    delta_w: f64,
    /// Scalar part of the `PhotonPathSeg` displacement applied to a `Photon`
    ///
    /// Intialized with a random value by `skin.new_path_seg_len_in_vxl()` and
    /// Decreased by `photon.moved_of()` until it reaches `0.`, after which a new one is created.
    pub len_in_vxl: f64,
    /// Directional part of the `PhotonPathSeg` displacement applied to a `Photon`
    ///
    /// Intialized at `{ ux: 0., uy: 0., uz: 1. }` with every new `Photon`
    /// Preserved as new `len_in_vxl` are random-generated reaches `0.` as the `Photon` is moved through `skin.absorption()`
    /// Changed by `dir.uz_reflected()`, `dir.iso_refrac()` and `dir.aniso_refrac()` during `skin.absorption()`
    pub dir: UnitVec,
}

impl PhotonPathSeg {
    /// Instantiates a new `PhotonPathSeg` with :
    /// - `len_in_vxl` passed in parameter,
    /// - `delta_w` of `0.` and
    /// - `dir` of `UnitVec { ux: 0., uy: 0., uz: 1. }`
    #[must_use]
    pub fn new(len_in_vxl: f64) -> Self {
        Self {
            len_in_vxl,
            ..Self::default()
        }
    }

    /// Changes the `.dir` part of the `PhotonPathSeg` with a `new_dir` passed in parameter.
    const fn set_dir(&mut self, new_dir: UnitVec) {
        self.dir.ux = new_dir.ux;
        self.dir.uy = new_dir.uy;
        self.dir.uz = new_dir.uz;
    }

    /// Returns `f64` based `dir.uz * len_in_vxl` converted to `i32` based `dz`.
    /// Since we *know* that `len_in_vxl` can be extremely big (ln(0) = ∞, we prevent `0.` but not
    /// `f64::EPSILON`) it is entirely likely that some returned `dx`, `dy` and or `dz` ends up being
    /// `i32::MIN` or `i32:MAX` on overflow by `f64` conversion to `i32`.
    fn dz(&self) -> i32 {
        #![allow(clippy::cast_possible_truncation)]

        //  as i32 returns i32:MAX or i32::MIN on overflow by f64,
        //  which will definitely result in `skin.try_skin_layer()` to return `None`.
        (self.len_in_vxl * self.dir.uz) as i32
    }

    /// Returns the whole `len_in_vxl` corresponding to a `dz` displacement along the `z` axis passed in parameter.
    fn len_in_vxl_from(&self, dz: i32) -> f64 {
        f64::from(dz) / self.dir.uz
    }
}

/// The `Photon` model used in the simulation
#[derive(Debug)]
struct Photon {
    /// Weight of this `Photon`, progressively reduced until it reaches `SEUIL_DE_POIDS_CRITIQUE` at the end of its simulation
    weight: f64,
    /// Conveniently tracks in which `skin_layer` this `Photon` is now located (derived from `Photon.pos`)
    skin_layer_kind: SkinLayerKind,
    /// Position of this `Photon`, in or out of the `Voxels`
    pub pos: VoxelPos,
    /// Random-length displacement to be applied to this `Photon` for the few absorption cycles it takes to reach `0.`
    pub path_seg: PhotonPathSeg,
}

impl Photon {
    /// Instantiates a new `Photon` with :
    /// - `weight` of `0.`,
    /// - `skin_layer_kind` of `SkinLayerKind::StratumCorneum`,
    /// - `pos.x` of `VXLS_X_SIZE / 2`,
    /// - `pos.y` of `VXLS_Y_SIZE / 2`,
    /// - `pos.z` of `0.`,
    /// - `path_seg.len_in_vxl` randomly initialized,
    /// - `path_seg.delta_w` of `0.`,
    /// - `path_seg.dir` of `UnitVec { ux: 0., uy: 0., uz: 1. }`
    ///
    pub fn new(skin: &Skin, rng: &mut ThreadRng) -> Self {
        let skin_layer_kind = StratumCorneum;
        let path_seg = PhotonPathSeg::new(skin.new_path_seg_len_in_vxl(skin_layer_kind, rng));

        Self {
            weight: 1.0,
            skin_layer_kind,
            pos: VoxelPos {
                x: VXLS_X_SIZE / 2,
                y: VXLS_Y_SIZE / 2,
                z: 0,
            },
            path_seg,
        }
    }

    /// Returns `.path_seg.delta_w`
    pub const fn path_seg_delta_w(&self) -> f64 {
        self.path_seg.delta_w
    }

    /// Returns `.weight - .path_seg.delta_w`
    pub fn adjusted_weight(&self) -> f64 {
        self.weight - self.path_seg.delta_w
    }

    /// Increases the `.path_seg.delta_w` by the `delta_w_coeff` passed in parameter, times the `photon.adjusted_weight()`
    pub fn increase_path_seg_delta_w_by(&mut self, delta_w_coeff: DeltaWCoeff) {
        self.path_seg.delta_w += delta_w_coeff.0 * self.adjusted_weight();
    }

    /// Sets `.weight` to the `.adjusted_weight()`, and reset `.path_seg.delta_w` to `0.`
    pub fn apply_and_reset_path_seg_delta_w(&mut self) {
        self.weight = self.adjusted_weight();
        self.path_seg.delta_w = 0.;
    }

    /// Scales this `photon.adjusted_weight()` (so, both its `.weight` and `path_seg.delta_w`),
    /// by a factor `m` passed in parameter
    pub fn scale_adjusted_weight_by(&mut self, m: f64) {
        self.weight *= m;
        self.path_seg.delta_w *= m;
    }

    /// Randomly generate a new `.path_seg.len_in_vxl` for this `photon` next simulation cycles.
    pub fn generate_next_path_seg(&mut self, skin: &Skin, rng: &mut ThreadRng) {
        self.path_seg.len_in_vxl = skin.new_path_seg_len_in_vxl(self.skin_layer_kind, rng);
    }

    /// - Adds the `delta_pos` of `moving_len` to .pos,
    /// - Adjusts .`skin_layer_kind` and `len_in_vxl` to `remaining_len`, and
    /// - Returns .`pos.is_within_xy_of_vxls()`
    pub fn move_of(&mut self, moving_len: f64, remaining_len: f64, skin_layer_kind: SkinLayerKind) -> bool {
        self.path_seg.len_in_vxl = remaining_len;
        self.skin_layer_kind = skin_layer_kind;
        // Adds the delta_pos to .pos and returns .pos.is_within_xy_of_vxls().
        self.pos.move_of(self.path_seg.dir.delta_pos(moving_len))
    }

    /// Returns the `SkinLayer` for this `photon.skin_layer_kind`
    pub const fn skin_layer(&self, skin: &Skin) -> &SkinLayer {
        skin.skin_layer(self.skin_layer_kind).layer
    }
}

/// Runs a Monte Carlo simulation
///
/// ### Parameters
///
/// - `nb_photons` - the number of photons to simulate
/// - `wavelength` -  the wavelength of the photons to simulate, in nm
/// - `chunk_size` - the size of the chunk of photons that will be passed to parallel processing instances of this simulation
///     - a `chunk_size` equal to `nb_photons` runs the simulation without parallelism
///     - a `chunk_size` equal to `(nb_photons / std::thread::available_parallelism()?.get()) + 1` runs the simulation
///       with the maximal parallel computing ressources available
/// - `verbose` - the verbosity of the simulation
///
/// ### Returns
///
/// `Some` `ndarray::Array2<f64>` of the simulation accumulated `Photon`s `weight` losses, as a slice of `Voxels` along
/// the XZ axis, or `None` if the `Wavelength` is not supported
#[must_use]
pub fn monte_carlo(
    nb_photons: usize,
    wavelength: Wavelength,
    chunk_size: usize,
    verbose: bool,
) -> Option<VoxelsSliceXZ> {
    /// Seuil de poids critique qui établit une limite à la simulation de chaque photon.
    const SEUIL_DE_POIDS_CRITIQUE: f64 = 0.01;
    /// En dessous d’un `SEUIL_DE_POIDS_CRITIQUE` prédéterminé, le paquet a une chance sur `M` de continuer à être simulé,
    /// sinon "il disparaı̂t". Pour éviter que de l’énergie disparaisse pris dans l'ensemble, lorsque le photon survie,
    /// on multiplie son énergie par `M` afin de compenser pour les fois ou l’énergie disparait.
    const M: f64 = 10.;
    /// Probabilité qu'un photon continue à être simulé lorsque son `weight` atteint le `SEUIL_DE_POIDS_CRITIQUE`
    const M_INV: f64 = 1. / M; //  Probabilité de ne pas disparitre

    let mut skin = Skin::default();
    skin.try_model_at(wavelength).map(|skin| {
        if verbose { println!("{skin:#?}" ); }
            (0..nb_photons)
            .into_par_iter()
             // Chunking the simulations reduces the overhead of creating many arrays
            .chunks(chunk_size)
            .map(|chunk| {
                // Each thread gets a local partial grid
                let mut local_vxls = Array2::<f64>::zeros(VXL_XZ_DIMS);

                // Each thread gets its own thread-local RNG.
                // This is fast and non-blocking.
                let mut rng = rand::rng();

                for nth_photon in chunk {
                    let mut  photon = Photon::new(skin, &mut rng);
                    #[allow(clippy::while_float)]
                    while SEUIL_DE_POIDS_CRITIQUE < photon.adjusted_weight() {

                        if 0. == photon.path_seg.len_in_vxl {
                            photon.generate_next_path_seg(skin, &mut rng);
                        }
                        if let Some((new_dir, delta_w_coeff)) = skin.absorption(&mut photon, &mut rng) {
                            photon.path_seg.set_dir(new_dir);

                            photon.increase_path_seg_delta_w_by(delta_w_coeff);
                            if 0. == photon.path_seg.len_in_vxl {
                                let pos = &photon.pos;
                                local_vxls[[pos.x(), pos.z()]] += photon.path_seg_delta_w();

                                photon.apply_and_reset_path_seg_delta_w();
                            }
                            //  if the photon weight is bellow a given critical threshold (Wc), play Roulette ! 
                            //  The photon as 1/M chance to continue simulation, otherwise "it disappears". To
                            //  avoid that overall energy "disappears" too, when a photon survives, its energy 
                            //  is multiplied by M so to compensate all the other cases where energy disappears.
                            #[allow(clippy::while_float)]                     //  Roulette !  (xsi_4)
                            if photon.adjusted_weight() < SEUIL_DE_POIDS_CRITIQUE  &&  rng.random::<f64>() < M_INV {
                                photon.scale_adjusted_weight_by(M);
                            }
                        }
                        else {  //  ! within(vxls) 
                            break;
                        }
                    }
                    if verbose && nth_photon.is_multiple_of(1_000) {
                        let pos = &photon.pos;
                        if pos.is_within_vxls() {
                            println!("{nth_photon}th photon : vxls[x:{}][z:{}] = {:?}", pos.x, pos.z, local_vxls[[pos.x(), pos.z()]]);
                        }
                        else {
                            println!("{nth_photon}th photon pos : {pos:?} out of vxls[{VXLS_X_SIZE}][{VXLS_Y_SIZE}][{VXLS_Z_SIZE}] !");
                        }
                    }
                }
                local_vxls
            })
            .reduce(|| Array2::<f64>::zeros(VXL_XZ_DIMS),
                |mut total, local| {
                    total += &local;
                    total
                },
            )
    })
}
