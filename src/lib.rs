
pub mod error;
pub use error::{Error, Result};

use std::{
    thread,
    ops::Range,
    sync::OnceLock,
    collections::HashMap,
    f64::consts::PI,
    time::Instant,
};
// use tracing::debug;
// use derive_more::{Index, IndexMut};
use strum::{EnumCount};
use rand::prelude::*;
use ndarray::{Array2, s};
// use ndarray::parallel::prelude::*;
use ruviz::prelude::*;
// use rayon::prelude::*;

// use uom::si::f64::*;
// use uom::si::{length::{millimeter, centimeter},reciprocal_length::{reciprocal_millimeter, reciprocal_centimeter}};
pub fn run() -> Result<()> {
    let available_parallelism = thread::available_parallelism()?.get();
    dbg!(&available_parallelism);
    
    partie_3_1();
    partie_1_2(available_parallelism)
}


fn partie_3_1() {

}

const NB_PHOTONS_1_2:usize = 1_000_000;

const MM_PER_CM:f64 = 10.;
const CM_PER_VXL:f64 = cm(0.001);
const MM_PER_VXL:f64 = CM_PER_VXL * MM_PER_CM;
const VXL_PER_MM:f64 = 1./MM_PER_VXL;

const INDICE_REFRAC_AIR:f64 = 1.000293; // à T.P.N. Unused, but for isometric solution.

const SKIN_LAYERS: [SkinLayer; SkinLayerKind::COUNT] = [
    SkinLayer {dz_in_mm: 0.02, indice_refrac: 1.42, v_b: 0.00, v_w: 0.05, kind: StratumCorneum},
    SkinLayer {dz_in_mm: 0.25, indice_refrac: 1.42, v_b: 0.00, v_w: 0.20, kind: Epiderme},
    SkinLayer {dz_in_mm: 0.10, indice_refrac: 1.39, v_b: 0.04, v_w: 0.50, kind: PapileDermique},
    SkinLayer {dz_in_mm: 0.08, indice_refrac: 1.39, v_b: 0.30, v_w: 0.60, kind: DermeSuperieur},
    SkinLayer {dz_in_mm: 0.20, indice_refrac: 1.39, v_b: 0.04, v_w: 0.70, kind: DermeReticulaire},
    SkinLayer {dz_in_mm: 0.30, indice_refrac: 1.39, v_b: 0.01, v_w: 0.70, kind: DermeProfond},
];
const SKIN_DELTA_Z_IN_MM:f64 = SKIN_LAYERS[StratumCorneum as usize].dz_in_mm  
    + SKIN_LAYERS[Epiderme as usize].dz_in_mm
    + SKIN_LAYERS[PapileDermique as usize].dz_in_mm
    + SKIN_LAYERS[DermeSuperieur as usize].dz_in_mm
    + SKIN_LAYERS[DermeReticulaire as usize].dz_in_mm
    + SKIN_LAYERS[DermeProfond as usize].dz_in_mm;

const VXLS_X_SIZE:i32 = 800;  //  8.0 mm
const VXLS_Y_SIZE:i32 = 800;  //  8.0 mm
const VXLS_Z_SIZE:i32 = (SKIN_DELTA_Z_IN_MM * VXL_PER_MM) as i32;
const VXL_XZ_DIMS:(usize, usize) = (VXLS_X_SIZE as usize, VXLS_Z_SIZE as usize);

const MAGIG_NUM_0: f64 = 2.303 * 150. / 64_500. / MM_PER_CM; //  "J'ai assumé que l'output était en cm^-1, donc divisé par 10 pour mm^-1"  

const SA_O2:f64 = 0.98;         //  Saturation en oxygène artériel de 98%
const SV_O2:f64 = 0.9 * SA_O2;  //  Saturation en oxygène veineux est 10% moindre que l'arteriel
const SEUIL_DE_POIDS_CRITIQUE:f64 = 0.01;

//  En dessous d’un weight seuil critique (Wc) prédéterminé, le paquet a une chance sur M de continuer à être simulé, 
//  sinon "il disparaı̂t". Pour éviter que de l’énergie disparaisse pris dans l'ensemble, lorsque le photon survie, 
//  on multiplie son énergie par M afin de compenser pour les fois ou l’énergie disparait.
const M:f64 = 10.;
const M_INV:f64 = 1. / M;   //  Probabilité de ne pas disparitre

/// Coefficient d'anisotropie
const G:f64 = 0.9;

const V_MEL:f64 = 0.1;
const MU_S:f64 = 20.1 * MM_PER_VXL;

const HB_SIZE:usize = 201;  //  Source Scott Prahl, lambda:nm, le reste : cm^-1/M
const HB_ARRAY : [(f64, f64); HB_SIZE] = [
    (3200.0, 14677.2),
    (2664.0, 13622.4),
    (2128.0, 12567.6),
    (1789.2, 11513.2),
    (1647.6, 10477.6),
    (1506.0, 9443.6),
    (1364.4, 8591.2),
    (1222.8, 7762.0),
    (1110.0, 7344.8),
    (1026.0, 6927.2),
    (942.0, 6509.6),
    (858.0, 6193.2),
    (774.0, 5906.8),
    (707.6, 5620.0),
    (658.8, 5366.8),
    (610.0, 5148.8),
    (561.2, 4930.8),
    (512.4, 4730.8),
    (478.8, 4602.4),
    (460.4, 4473.6),
    (442.0, 4345.2),
    (423.6, 4216.8),
    (405.2, 4088.4),
    (390.4, 3965.08),
    (379.2, 3857.6),
    (368.0, 3750.12),
    (356.8, 3642.64),
    (345.6, 3535.16),
    (335.2, 3427.68),
    (325.6, 3320.2),
    (319.6, 3226.56),
    (314.0, 3140.28),
    (308.4, 3053.96),
    (302.8, 2967.68),
    (298.0, 2881.4),
    (294.0, 2795.12),
    (290.0, 2708.84),
    (285.6, 2627.64),
    (282.0, 2554.4),
    (279.2, 2481.16),
    (277.6, 2407.92),
    (276.0, 2334.68),
    (274.4, 2261.48),
    (272.8, 2188.24),
    (274.4, 2115.0),
    (276.0, 2051.96),
    (277.6, 2000.48),
    (279.2, 1949.04),
    (282.0, 1897.56),
    (286.0, 1846.08),
    (290.0, 1794.28),
    (294.0, 1741.0),
    (298.0, 1687.76),
    (302.8, 1634.48),
    (308.4, 1583.52),
    (314.0, 1540.48),
    (319.6, 1497.4),
    (325.2, 1454.36),
    (332.0, 1411.32),
    (340.0, 1368.28),
    (348.0, 1325.88),
    (356.0, 1285.16),
    (364.0, 1244.44),
    (372.4, 1203.68),
    (381.2, 1152.8),
    (390.0, 1102.2),
    (398.8, 1102.2),
    (407.6, 1102.2),
    (418.8, 1101.76),
    (432.4, 1100.48),
    (446.0, 1115.88),
    (459.6, 1161.64),
    (473.2, 1207.4),
    (487.6, 1266.04),
    (502.8, 1333.24),
    (518.0, 1405.24),
    (533.2, 1515.32),
    (548.4, 1541.76),
    (562.0, 1560.48),
    (574.0, 1560.48),
    (586.0, 1548.52),
    (598.0, 1508.44),
    (610.0, 1459.56),
    (622.8, 1410.52),
    (636.4, 1361.32),
    (650.0, 1311.88),
    (663.6, 1262.44),
    (677.2, 1213.0),
    (689.2, 1163.56),
    (699.6, 1114.8),
    (710.0, 1075.44),
    (720.4, 1036.08),
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
    (1001.2, 692.64),
    (1011.6, 692.48),
    (1022.0, 692.36),
    (1032.4, 692.2),
    (1042.8, 691.96),
    (1050.0, 691.76),
    (1054.0, 691.52),
    (1058.0, 691.32),
    (1062.0, 691.08),
    (1066.0, 690.88),
    (1072.8, 690.64),
    (1082.4, 692.44),
    (1092.0, 694.32),
    (1101.6, 696.2),
    (1111.2, 698.04),
    (1118.4, 699.92),
    (1123.2, 701.8),
    (1128.0, 705.84),
    (1132.8, 709.96),
    (1137.6, 714.08),
    (1142.8, 718.2),
    (1148.4, 722.32),
    (1154.0, 726.44),
    (1159.6, 729.84),
    (1165.2, 733.2),
    (1170.0, 736.6),
    (1174.0, 739.96),
    (1178.0, 743.6),
    (1182.0, 747.24),
    (1186.0, 750.88),
    (1190.0, 754.52),
    (1194.0, 758.16),
    (1198.0, 761.84),
    (1202.0, 765.04),
    (1206.0, 767.44),
    (1209.2, 769.8),
    (1211.6, 772.16),
    (1214.0, 774.56),
    (1216.4, 776.92),
    (1218.8, 778.4),
    (1220.8, 778.04),
    (1222.4, 777.72),
    (1224.0, 777.36),
    (1225.6, 777.04),
    (1227.2, 776.64),
    (1226.8, 772.36),
    (1224.4, 768.08),
    (1222.0, 763.84),
    (1219.6, 752.28),
    (1217.2, 737.56),
    (1215.6, 722.88),
    (1214.8, 708.16),
    (1214.0, 693.44),
    (1213.2, 678.72),
    (1212.4, 660.52),
    (1210.4, 641.08),
    (1207.2, 621.64),
    (1204.0, 602.24),
    (1200.8, 583.4),
    (1197.6, 568.92),
    (1194.0, 554.48),
    (1190.0, 540.04),
    (1186.0, 525.56),
    (1182.0, 511.12),
    (1178.0, 495.36),
    (1173.2, 473.32),
    (1167.6, 451.32),
    (1162.0, 429.32),
    (1156.4, 415.28),
    (1150.8, 402.28),
    (1144.0, 389.288),
    (1136.0, 374.944),
    (1128.0, 359.656),
    (1120.0, 344.372),
    (1112.0, 329.084),
    (1102.4, 313.796),
    (1091.2, 298.508),
    (1080.0, 283.22),
    (1068.8, 267.932),
    (1057.6, 252.648),
    (1046.4, 237.36),
    (1035.2, 222.072),
    (1024.0, 206.784)
];

fn hb_at_index(key: i32) -> Option<(f64, f64)> {
    let i = (key-600)/2;
    if 0<= i && i < (HB_SIZE as i32) { Some(HB_ARRAY[i as usize]) } else { None }
}
static WATER_HASHMAP: OnceLock<HashMap<i32, f64>> = OnceLock::new();

fn water_hashmap(key: i32) -> Option<f64> {
    WATER_HASHMAP.get_or_init(|| [   
        (600, 0.0020184),
        (605, 0.0023501),
        (610, 0.0025524),
        (615, 0.0027167),
        (619, 0.0028383),
        (625, 0.0029587),
        (629, 0.0029984),
        (635, 0.0030699),
        (640, 0.0030841),
        (646, 0.0031255),
        (650, 0.0032358),
        (655, 0.0034113),
        (661, 0.0036898),
        (665, 0.0038362),
        (670, 0.0039355),
        (675, 0.0040559),
        (681, 0.0042454),
        (685, 0.0045298),
        (690, 0.0048303),
        (695, 0.0053574),
        (700, 0.0060120),
        (705, 0.0073112),
        (710, 0.0088510),
        (715, 0.010544),
        (719, 0.012736),
        (724, 0.015850),
        (729, 0.019810),
        (735, 0.023063),
        (740, 0.024773),
        (745, 0.025818),
        (750, 0.026125),
        (755, 0.026294),
        (760, 0.026114),
        (766, 0.025770),
        (769, 0.024950),
        (775, 0.023981),
        (780, 0.022706),
        (785, 0.021429),
        (791, 0.020374),
        (794, 0.019902),
        (800, 0.019640),
        (805, 0.019815),
        (809, 0.020657),
        (815, 0.022335),
        (820, 0.024829),
        (824, 0.027737),
        (830, 0.030905),
        (836, 0.033732),
        (840, 0.036808),
        (845, 0.039990),
        (849, 0.043343),
        (855, 0.046336),
        (859, 0.048978),
        (865, 0.051515),
        (871, 0.054074),
        (875, 0.056111),
        (879, 0.057942),
        (885, 0.060113),
        (889, 0.062224),
        (895, 0.064867),
        (900, 0.067924),
        (906, 0.071455),
        (910, 0.078707),
        (914, 0.092052),
        (920, 0.11338),
        (925, 0.14405),
        (929, 0.18505),
        (935, 0.23792),
        (940, 0.29005),
        (944, 0.34035),
        (951, 0.38759),
        (955, 0.41976),
        (959, 0.43984),
        (966, 0.45057),
        (971, 0.45345),
        (975, 0.44852),
        (980, 0.43851),
        (984, 0.42603),
        (991, 0.41258),
        (995, 0.39527),
        (1000, 0.37699),
        (1009, 0.33477),
        (1021, 0.28948),
        (1030, 0.24413),
        (1040, 0.20420),
        (1050, 0.16983),
        (1059, 0.15414),
        (1069, 0.14800),
        (1079, 0.15478),
        (1089, 0.17297),
        (1099, 0.19530),
        (1109, 0.23093),
        (1119, 0.29512),
        (1130, 0.43026),
        (1140, 0.65599),
        (1151, 1.0160),
        (1159, 1.1591),
        (1169, 1.2040),
    ].into_iter().collect())
        .get(&key)
        .copied()
}

struct Vxls {
    x_range: Range<i32>,
    y_range: Range<i32>,
    z_range: Range<i32>,
}

impl Vxls {
    pub const fn new(x_size: i32, y_size: i32, z_size: i32) -> Self {
        Self {
            x_range : 0..x_size,
            y_range : 0..y_size,
            z_range : 0..z_size,
        }
    }
    fn x_contains(&self, x: &i32) -> bool {
        self.x_range.contains(x)
    }
    fn y_contains(&self, y: &i32) -> bool {
        self.y_range.contains(y)
    }
    fn z_contains(&self, z: &i32) -> bool {
        self.z_range.contains(z)
    }
}
const VXLS:Vxls = Vxls::new(VXLS_X_SIZE, VXLS_Y_SIZE, VXLS_Z_SIZE);

fn plot(vxls:&mut Array2<f64>, wavelength: i32, k_photon: usize) -> Result<()> {

    vxls.par_mapv_inplace(f64::log10);
    let binding = vxls.t();
    let view = binding.slice(s![..;-1, ..]);
       
    let x_axis_scaling = (VXLS_X_SIZE as f64) * MM_PER_VXL;
    let z_axis_scaling = (VXLS_Z_SIZE as f64) * MM_PER_VXL;
        
    let start = Instant::now();
    
    Plot::new()
        .heatmap(&view, Some(HeatmapConfig::new()
            .colorbar(true)
            .colorbar_label("Énergie absorbée")
            //  No way to scale the colorbar axis in log yet
            ))
        .xlim(0., x_axis_scaling)
        .xlabel("Position en x(mm)")
        .ylim(0., z_axis_scaling)    
        // .ylim(z_axis_scaling, 0.)    
        .ylabel("Profondeur (mm)")
        .title(format!("{k_photon}k photon at {wavelength}nm"))
        .save(format!(".heatmap-{wavelength}nm{k_photon}kphoton.png"))?;
    
    println!("Plot finished in {:?}", start.elapsed());

    Ok(())
}

fn partie_1_2(available_parallelism: usize) -> Result<()> {
    let chunk_size = (NB_PHOTONS_1_2 / available_parallelism) + 1;
    dbg!(&chunk_size);
    
    dbg!(&SKIN_DELTA_Z_IN_MM);

    dbg!("partie_1_2: before monte_carlo()");
    let start = Instant::now();

    let wavelength = 700;
    let k_photon = NB_PHOTONS_1_2/1_000;
    if let Some(mut vxls) = monte_carlo(wavelength, chunk_size) {
        println!("monte_carlo of {k_photon}k photon at {wavelength}nm finished in {:?}", start.elapsed());

        plot(&mut vxls, wavelength, k_photon)?;
    }
    else { 
        println!("monte_carlo of {k_photon}k photon at {wavelength}nm finished in {:?}", start.elapsed()); 
    }

    Ok(())
}

const fn cm(cm :f64) -> f64 {
    cm
}
const fn mm_to_vxl(mm: f64) -> i32 {
    (mm * VXL_PER_MM) as i32
}

#[derive(Debug, EnumCount, Clone, Copy, PartialEq)]
enum SkinLayerKind {
    StratumCorneum,
    Epiderme,
    PapileDermique,
    DermeSuperieur,
    DermeReticulaire,
    DermeProfond,
}
use SkinLayerKind::*;

#[derive(Debug, Clone, Copy)]
struct SkinLayer {
    kind: SkinLayerKind,
    dz_in_mm: f64,
    indice_refrac : f64,
    v_b : f64,
    v_w : f64,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct DeltaWCoeff(f64);

#[derive(Debug)]
struct SkinLayerAtWl {
    layer: &'static SkinLayer,
    
    pub indice_refrac_ratio: f64,
    pub vxl_z_range : Range<i32>,
    
    pub mu_t: f64,          //  in reciprocal_vxl or vxl^-1
    pub mu_a_on_mu_t: DeltaWCoeff,
}

impl SkinLayerAtWl {
    fn kind(&self) -> SkinLayerKind {
        self.layer.kind
    }

    fn indice_refrac(&self) -> f64 {
        self.layer.indice_refrac
    }
    
    fn model_at(&mut self, wavelength: f64, val_hb_0: f64, val_hb_1: f64, val_water: f64) {
        //  "J'ai assumé que l'output était en cm^-1, donc divisé par 10 pour mm^-1"
        let mu_a_baseline=  7.84e-7 * wavelength.powf(-3.255) / MM_PER_CM;   //  reciprocal_mm
        
        let mu_absorption_oxyhemoglobine = val_hb_0 * MAGIG_NUM_0; // reciprocal_mm
        let mu_a_hemoglobine = val_hb_1 * MAGIG_NUM_0;  //  reciprocal_mm
    
        let mu_w = val_water / MM_PER_CM;  //  reciprocal_mm

        //  Coeff artere & veine
        let (mu_art, mu_vei) = (
            SA_O2 * mu_absorption_oxyhemoglobine + (1. - SA_O2) * mu_a_hemoglobine,
            SV_O2 * mu_absorption_oxyhemoglobine + (1. - SV_O2) * mu_a_hemoglobine );

        let SkinLayer { kind, v_w, v_b, .. } = *self.layer;
        
        let mu_a = MM_PER_VXL * match kind {
            StratumCorneum | Epiderme =>  
                V_MEL * (6.6e+10 * wavelength.powf(-3.33) / MM_PER_CM)  +  v_w * mu_w  +  (1. - (V_MEL + v_w)) * mu_a_baseline,
            _ => {  //  Tout sous l'Epiderme
                let v_a = v_b / 2.;     //  split 50-50 with v_v
                let v_v = v_a;          //  split 50-50 with v_a
                v_a * mu_art  +  v_v * mu_vei  +  v_w * mu_w  +  (1.- (v_a + v_v + v_w)) * mu_a_baseline
            },                    
        };
        self.mu_t = mu_a + MU_S;
        self.mu_a_on_mu_t = DeltaWCoeff(mu_a / self.mu_t);
    }
    
}

#[derive(Debug, )]
struct Skin {
    g_square: f64,
    layers: [SkinLayerAtWl; SkinLayerKind::COUNT],
    wavelength_i : i32,
}

const LOWEST_WAVELENGTH:i32 = 600;  // nanometer

impl Default for Skin {
    fn default() -> Self {
        let mut dz_in_vxl_acc = 0;
        let mut indice_refrac_denominateur = INDICE_REFRAC_AIR;
        
        Self { 
            g_square: G.powi(2),
            layers: SKIN_LAYERS.iter().map(|skin_layer| { 
                let SkinLayer { dz_in_mm, indice_refrac, .. } = skin_layer;

                let dz_in_vxl_start = dz_in_vxl_acc; 
                dz_in_vxl_acc += mm_to_vxl(*dz_in_mm);
                
                let indice_refrac_numerateur = indice_refrac_denominateur /* precedent*/ ;
                indice_refrac_denominateur = *indice_refrac;
                
                SkinLayerAtWl {
                    layer: skin_layer,
                    // kind, /*dz_in_mm,*/ indice_refrac, v_b, v_w,
                    
                    indice_refrac_ratio: indice_refrac_numerateur / indice_refrac_denominateur,
                    vxl_z_range: dz_in_vxl_start..dz_in_vxl_acc,
                    
                    mu_t: 0.,
                    mu_a_on_mu_t: DeltaWCoeff::default(),
                }
            }).collect::<Vec<_>>().try_into().expect("SKIN_LAYERS [SkinLayer; SkinLayerKind::COUNT] to SkinLayerAtWl[SkinLayerAtWl; SkinLayerKind::COUNT]"),
            wavelength_i: LOWEST_WAVELENGTH,
        }
    }
}

impl Skin {
    fn try_model_at(&mut self, wavelength_i: i32) -> Option<&mut Self> {
        if let Some((val_hb_0, val_hb_1)) = hb_at_index(wavelength_i)
        && let Some(val_water) = water_hashmap(wavelength_i) {
            
            let wavelength = wavelength_i as f64;
            
            for skin_layer_at_wl in self.layers.iter_mut() {
                skin_layer_at_wl.model_at(wavelength, val_hb_0, val_hb_1, val_water);
            }
            self.wavelength_i = wavelength_i;
            Some(self)
        } else { None }
    }
    
    fn skin_layer(&self, skin_layer_kind: SkinLayerKind) -> &SkinLayerAtWl {
        &self.layers[skin_layer_kind as usize]
    }
    fn new_path_seg_len_in_vxl(&self, src_skin_layer_kind: SkinLayerKind, rng: &mut ThreadRng) -> f64 {
        //  parcours du groupe de photon pour cette itération
        let mut xsi_1:f64 = rng.random::<f64>();
        while xsi_1 == 0. {             //  Exclude 0 from the interval because ln(0) = ∞, or NaN.
            xsi_1 = rng.random::<f64>();
        }
        //  -xsi_1.ln() runs from 0. to very large (f64::MAX ?).
        -xsi_1.ln() / self.skin_layer(src_skin_layer_kind).mu_t     //  in vxl
    }
    
    fn try_skin_layer(&self, photon_pos_z: i32) -> Option<&SkinLayerAtWl> {
        self.layers.iter()
            .find(|skin_layer| skin_layer.vxl_z_range.contains(&photon_pos_z))
    }
    
    fn same_refrac(&self, dir: &UnitVec, rng: &mut ThreadRng) -> UnitVec {
        // xsi_2
        let phi = 2.0 * PI * rng.random::<f64>();  
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        
        //  xsi_3
        let mut cos_theta = 1. / (2.* G) * (1. + self.g_square - ((1. - self.g_square)/(1. - G + 2. * G * rng.random::<f64>())).powi(2)); 
        //  Ici aussi, pour éviter les erreurs arithmétiques flottantes cos doit être entre (-1, 1)
        cos_theta = cos_theta.clamp(-1., 1.);
        let sin_theta = co_sin(cos_theta);
        
        if dir.uz.abs() == 1. {
            UnitVec::new(   
                sin_theta * sin_phi,    //  dir.ux
                sin_theta * cos_phi,    //  dir.uy
                cos_theta,              //  dir.uz
            )
        } else {
            let not_uz = co_sin(dir.uz);
            let uz_cos_phi = dir.uz * cos_phi;
            UnitVec::new(   
                sin_theta * (dir.ux * uz_cos_phi  -  dir.uy * sin_phi) / not_uz  +  cos_theta * dir.ux,
                sin_theta * (dir.uy * uz_cos_phi  -  dir.ux * sin_phi) / not_uz  +  cos_theta * dir.uy,
                dir.uz * cos_theta  -  not_uz * sin_theta * cos_phi
            )
        }
    }
    
    fn diff_refrac(&self, n1: f64, n2: f64, n1_on_n2: f64, dir: &UnitVec, rng: &mut ThreadRng) -> (UnitVec, DeltaWCoeff) {
        let cos_theta_i = dir.uz.abs().clamp(0., 1.);   //  cos_theta_i in [0., 1.] range
        let sin_theta_i = co_sin(cos_theta_i);              //  sin_theta_i in [0., 1.] range too.
        
        let sin_theta_t = (n1_on_n2 * sin_theta_i).clamp(0.0, 1.0);  //  Loi de Snell-Descartes

        //  Calcul de réflexion
        if 1.0 == sin_theta_t {     //  Complete internal reflexion
            //  rng.random::<f64>() generates a random number over 0.0..1.0, so always < 1.0 
            //  therefore if 1.0 == sin_theta_t, r= sin_theta_t is guarantied > rng.random::<f64>()
            //  without even having to call it. 
            //  Also r_square is 1.0, and delta_w_coeff = (1. - r_square) = 0.
            (dir.uz_reflected(), DeltaWCoeff::default())
        } else {
            let cos_theta_t = co_sin(sin_theta_t);
            let r =  (n1 * cos_theta_i  -  n2 * cos_theta_t) / (n1 * cos_theta_i  +  n2 * cos_theta_t);
            let r_square = r.powi(2);
            
            if rng.random::<f64>() < r {                   //  xsi_5 <  r  : Reflexion
                (dir.uz_reflected(), DeltaWCoeff(1. - r_square) )
            } else {                                    //  xsi_5 >= r  : Transmission    
                (   UnitVec::new(
                        dir.ux * n1_on_n2,  //  dir.ux
                        dir.uy * n1_on_n2,  //  dir.uy
                        cos_theta_t,        //  dir.uz
                    ),              DeltaWCoeff(r_square))
            }                   
        }
    }

    fn absorption(&self, photon: &mut Photon, rng: &mut ThreadRng) -> Option<(UnitVec, DeltaWCoeff)> {
        // let Photon { pos, path_seg, .. } = photon; 
        let SkinLayer { kind:src_skin_layer_kind, indice_refrac:n1, ..} = *photon.skin_layer(self);
        let path_seg = photon.path_seg;
        let src_pos_z = photon.pos.z;
        
            // Some(dst_skin_layer) means within voxels in z 
        if let Some(dst_skin_layer) = self.try_skin_layer(src_pos_z + photon.path_seg.dz()) {
            
            let mu_a_on_mu_t = dst_skin_layer.mu_a_on_mu_t;
            let dst_skin_layer_kind = dst_skin_layer.kind();
            
            //  Trans skin layers 
            if src_skin_layer_kind != dst_skin_layer_kind {
                //  Stop at skin layer transition for this iteration.
                
                let moving_len = path_seg.new_len_in_vxl(dst_skin_layer.vxl_z_range.start - src_pos_z);
                //  returns photon.pos.is_within_xy_of_vxls()
                if photon.move_of(moving_len, path_seg.len_in_vxl - moving_len, dst_skin_layer_kind) { 
                    // skin layer transition is always within voxels in z 
                
                    //  skin layers with different indice_refrac: the path_seg.dir will change
                    if 1.0 != dst_skin_layer.indice_refrac_ratio {
                        Some(self.diff_refrac(n1, dst_skin_layer.indice_refrac(), 
                                            dst_skin_layer.indice_refrac_ratio, 
                                            &photon.path_seg.dir, rng))
                    }
                    else {
                        Some((self.same_refrac(&photon.path_seg.dir, rng), mu_a_on_mu_t))
                    }
                }
                else { None }   //  not within voxels.
            }
            else {  //  photon stays in the same skin layer
                //  returns photon.pos.is_within_xy_of_vxls()
                if photon.move_of(path_seg.len_in_vxl, 0., dst_skin_layer_kind) { 
                    // photon is staying within voxels in z

                    Some((self.same_refrac(&path_seg.dir, rng), mu_a_on_mu_t))
                }
                else { None }   //  not within voxels.
            }
        }   // ! Some(dst_skin_layer) means not within voxels in z 
        else { None }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct VoxelPos {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}
impl VoxelPos {
    fn is_within_vxls(&self) -> bool {
          VXLS.x_contains(&self.x) && VXLS.y_contains(&self.y) && VXLS.z_contains(&self.z)
    }
    fn is_within_xy_of_vxls(&self) -> bool {
          VXLS.x_contains(&self.x) && VXLS.y_contains(&self.y)
    }
    /// Adds the delta_pos to self and returns .is_within_xy_of_vxls().
    fn move_of(&mut self, delta_pos: DeltaVoxelPos) -> bool {
        self.x += delta_pos.dx;
        self.y += delta_pos.dy;
        self.z += delta_pos.dz;
        self.is_within_xy_of_vxls()
    }
}
#[derive(Debug, Clone, Copy)]
struct DeltaVoxelPos {
    dx: i32,
    dy: i32,
    dz: i32,
}

#[derive(Debug, Clone, Copy)]
pub struct UnitVec{
    pub ux: f64,
    pub uy: f64,
    pub uz: f64,
}
impl Default for UnitVec {
    fn default() -> Self {
        Self {ux: 0., uy: 0., uz: 1.,}
    }
}

impl UnitVec {
    fn new(ux: f64, uy: f64, uz: f64) -> Self {
        Self {ux, uy, uz}
    }
    
    fn delta_pos(&self, len: f64) -> DeltaVoxelPos {
        //  as i32 returns i32:MAX or i32::MIN on overflow by f64, 
        //  which will definitely result in !photon.pos.is_within(&vxls)
        DeltaVoxelPos {
            dx: (self.ux * len) as i32,
            dy: (self.uy * len) as i32,
            dz: (self.uz * len) as i32,
        }
    }
    fn uz_reflected(&self) -> Self {
        //  Only .uz is reflected the other unitary element ux and uy stay the same.
        Self {
            ux:  self.ux,
            uy:  self.uy,
            uz: -self.uz,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct PhotonPathSeg {
    delta_w : f64,
    pub len_in_vxl : f64,
    pub dir : UnitVec,
}

impl PhotonPathSeg {
    pub fn new(len_in_vxl: f64) -> Self {
        Self {
            len_in_vxl,
            .. PhotonPathSeg::default()
        }
    }
    fn set_dir(&mut self, new_dir: UnitVec) {
        self.dir.ux = new_dir.ux;
        self.dir.uy = new_dir.uy;
        self.dir.uz = new_dir.uz;
    }
    fn dz(&self) -> i32 {
        (self.len_in_vxl * self.dir.uz) as i32
    }
    fn new_len_in_vxl(&self, dz: i32) -> f64 {
        dz as f64 / self.dir.uz
    }
}
#[derive(Debug)]
struct Photon {
    // skin: &'a mut Skin,
    weight: f64,
    skin_layer_kind: SkinLayerKind,
    pub pos: VoxelPos,
    pub path_seg: PhotonPathSeg,
}

impl Photon {
    pub fn new(skin: &Skin, rng: &mut ThreadRng) -> Self {
        let skin_layer_kind= StratumCorneum;
        let path_seg = PhotonPathSeg::new(skin.new_path_seg_len_in_vxl(skin_layer_kind, rng));

        Self {
            weight: 1.0,
            pos: VoxelPos { 
                x: VXLS_X_SIZE / 2,
                y: VXLS_Y_SIZE / 2,
                z: 0,
            },
            path_seg,
            skin_layer_kind
        }
    }
    
    pub fn path_seg_delta_w(&self) -> f64 {
        self.path_seg.delta_w
    }
    pub fn adjusted_weight(&self) -> f64 {
        self.weight - self.path_seg.delta_w
    }
    pub fn increase_path_seg_delta_w_by(&mut self, delta_w_coeff: DeltaWCoeff) {
        self.path_seg.delta_w += delta_w_coeff.0 * self.adjusted_weight();
    }
    pub fn apply_and_reset_path_seg_delta_w(&mut self)  {
        self.weight = self.adjusted_weight();
        self.path_seg.delta_w = 0.;
    }
    pub fn scale_weight_by(&mut self, m: f64)  {
        self.weight *= m;
        self.path_seg.delta_w *= m;
    }

    pub fn generate_next_path_seg(&mut self, skin: &Skin, rng: &mut ThreadRng) {
        self.path_seg.len_in_vxl = skin.new_path_seg_len_in_vxl(self.skin_layer_kind, rng);
    }

    /// - Adds the delta_pos of moving_len to .pos, 
    /// - Adjusts .skin_layer_kind and len_in_vxl to remaining_len, and 
    /// - Returns .pos.is_within_xy_of_vxls()
    pub fn move_of(&mut self, moving_len:f64, remaining_len: f64, skin_layer_kind: SkinLayerKind) -> bool {
        self.path_seg.len_in_vxl = remaining_len;
        self.skin_layer_kind = skin_layer_kind;
        // Adds the delta_pos to .pos and returns .pos.is_within_xy_of_vxls().
        self.pos.move_of(self.path_seg.dir.delta_pos(moving_len))
    }
    
    pub fn skin_layer(&self, skin: &Skin) -> &SkinLayer {
        skin.skin_layer(self.skin_layer_kind).layer
    }
}

/// Return sin if cos and cos if sin
fn co_sin(sin_or_cos:f64) -> f64 {
    (1. - sin_or_cos.powi(2)).sqrt()
}

/// wavelength in nm 
fn monte_carlo(wavelength_i:i32, _chunk_size: usize) -> Option<Array2<f64>>{
    let mut skin = Skin::default();
    let mut local_vxls = Array2::<f64>::zeros(VXL_XZ_DIMS);
    skin.try_model_at(wavelength_i).map(|skin| {
        println!("{:#?}", &skin);
        // let _:i32 = (0..NB_PHOTONS_1_2).into_par_iter().map(|nth_photon| {
        let _ =
            (0..NB_PHOTONS_1_2)
            // .into_par_iter()
            //  // Chunking the simulations reduces the overhead of creating many arrays
            // .chunks(chunk_size)
            .map(|nth_photon| {
                // Each thread gets a local partial grid
                // let mut local_vxls = Array2::<f64>::zeros(VXL_XZ_DIMS);
                
                // Each thread gets its own thread-local RNG.
                // This is fast and non-blocking.
                let mut rng = rand::rng();  
                
                // for _ in chunk {
                    let mut  photon = Photon::new(skin, &mut rng);
                    // if nth_photon % 1_000 == 0 { dbg!(&nth_photon, &photon); }
                    
                    while SEUIL_DE_POIDS_CRITIQUE < photon.adjusted_weight() {
                        
                        if 0. == photon.path_seg.len_in_vxl {
                            photon.generate_next_path_seg(skin, &mut rng);    
                        }
                        if let Some((new_dir, delta_w_coeff)) = skin.absorption(&mut photon, &mut rng) {
                            photon.path_seg.set_dir(new_dir);
                            
                            photon.increase_path_seg_delta_w_by(delta_w_coeff);
                            if 0. == photon.path_seg.len_in_vxl {
                                let pos = &photon.pos;
                                local_vxls[[pos.x as usize, pos.z as usize]] += photon.path_seg_delta_w();
                                
                                photon.apply_and_reset_path_seg_delta_w(); 
                            }
                            
                            //  if the photon weight is bellow a given critical threshold (Wc), play Roulette ! 
                            //  The photon as 1/M chance to continue simulation, otherwise "it disappears". To
                            //  avoid that overall energy "disappears" too, when a photon survives, its energy 
                            //  is multiplied by M so to compensate all the other cases where energy disappears.
                                                                            //  Roulette !  (xsi_4)
                            if photon.adjusted_weight() < SEUIL_DE_POIDS_CRITIQUE  &&  rng.random::<f64>() < M_INV {
                                photon.scale_weight_by(M); 
                            }
                        }
                        else {  //  ! within(vxls) 
                            break;
                        }
                    }
                    if nth_photon % 1_000 == 0 {
                        let pos = &photon.pos;
                        if pos.is_within_vxls() {
                            println!("{nth_photon}th photon : vxls[x:{}][z:{}] = {:?}", pos.x, pos.z, local_vxls[[pos.x as usize, pos.z as usize]]);
                        }
                        else {
                            println!("{nth_photon}th photon pos : {pos:?} out of vxls[{VXLS_X_SIZE}][{VXLS_Y_SIZE}][{VXLS_Z_SIZE}] !");
                        }
                    }
                // }
                // local_vxls
                0
            })
            .sum::<i32>();
            local_vxls
            
            // .reduce(|| Array2::<f64>::zeros(VXL_XZ_DIMS), 
            //     |mut total, local| {
            //         total += &local; 
            //         total
            //     },
            // )
            
    })
}

