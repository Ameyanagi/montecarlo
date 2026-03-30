
pub mod error;
pub use error::{Error, Result};

use std::{
    ops::Range,
    sync::OnceLock,
    collections::HashMap,
    f64::consts::PI,
    // cmp::min,    //  Not available for float !!! 
};
use tracing::debug;
use derive_more::{Index, IndexMut};
use strum::{EnumCount};

// use uom::si::f64::*;
// use uom::si::{length::{millimeter, centimeter},reciprocal_length::{reciprocal_millimeter, reciprocal_centimeter}};
// use uom::si::length::millimeter;
// use uom::si::time::second;

use rand::prelude::*;
use ruviz::prelude::*;

pub fn run() -> Result<()> {
    partie_3_1();
    partie_1_2()
}


fn partie_3_1() {
    // let dz = Length::new::<centimeter>(0.01);

    // let mu_a = 20. / dz;
    // let mu_s = 0./ Length::new::<centimeter>(1.);
    // let mu_t = mu_a + mu_s;

    // let mu_x = 1.;
    // let w = 1.;

    // let m = 10.;

    // let w_c = 0.1;

}

const NB_PHOTONS_1_2:i32 = 100;

const UNIT_SCALE:f64 = 10.;
const CM_PER_VXL:f64 = cm(0.001);
const MM_PER_VXL:f64 = CM_PER_VXL * UNIT_SCALE;
const VXL_PER_MM:f64 = 1./MM_PER_VXL;

const VXLS_X_SIZE:i32 = 2_000;  //  2.0 cm
const VXLS_Y_SIZE:i32 =   800;  //  0.8 cm
const VXLS_Z_SIZE:i32 = 1_300;  //  1.3 cm

const MAGIG_NUM_0: f64 = 2.303 * 150. / 64_500. / UNIT_SCALE; //  "J'ai assumé que l'output était en cm^-1, donc divisé par 10 pour mm^-1"  

const SA_O2:f64 = 0.98;         //  Saturation en oxygène artériel de 98%
const SV_O2:f64 = 0.9 * SA_O2;  //  Saturation en oxygène veineux est 10% moindre que l'arteriel
const SEUIL_DE_POIDS_CRITIQUE:f64 = 0.01;

//  En dessous d’un poids seuil critique (Wc) prédéterminé, le paquet a une chance sur M de continuer à être simulé, 
//  sinon "il disparaı̂t". Pour éviter que de l’énergie disparaisse pris dans l'ensemble, lorsque le photon survie, 
//  on multiplie son énergie par M afin de compenser pour les fois ou l’énergie disparait.
const M:f64 = 10.;
const M_INV:f64 = 1. / M;   //  Probabilité de ne pas disparitre

const N1:f64 = 1.42;    //  indice de réfraction du stratum et de l'épiderme
const N2:f64 = 1.39;    //  indice de réfraction du derme
const N1_SUR_N2:f64 = N1/N2;    //  > 1   

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

#[derive(Index, IndexMut, Debug )]
struct VoxelsXZ {
    x_range: Range<i32>,
    y_range: Range<i32>,
    z_range: Range<i32>,
    #[index]
    #[index_mut]
    vxls:Vec<Vec<f64>>,
}
impl VoxelsXZ {
    fn new(x_size: i32, y_size: i32, z_size: i32) -> Self {
        let x_range = 0..x_size;
        let y_range = 0..y_size;
        let z_range = 0..z_size;
        Self {
            x_range,
            y_range,
            z_range,
            vxls: (0..x_size).map(|_|  
                           (0..z_size).map(|_| 0.).collect())
                        .collect(),
        }
    }
    
    fn hors_jeu(&self, pos: &VoxelPos) -> bool {
          !self.x_range.contains(&pos.x) || !self.y_range.contains(&pos.y) || !self.z_range.contains(&pos.z)
    }
}

impl NumericData2D for VoxelsXZ {
    fn shape(&self) -> (usize, usize) {
        (self.x_range.end as usize, self.z_range.end as usize)
    }
    fn try_collect_row_major_f64(&self) -> std::result::Result<Vec<f64>, ruviz::core::PlottingError> {
        let mut values = Vec::with_capacity((self.x_range.end * self.z_range.end) as usize);
        for row in &self.vxls {
            //  In Voxels the y axis is displayed in revert and we want the log of the values to be heatmaped 
            values.extend(row.into_iter().rev().map(|f| f.log10()));
        }
        Ok(values)
    }
}

fn partie_1_2() -> Result<()> {
    dbg!("partie_1_2: before Voxels allocation");
    
    let mut a = VoxelsXZ::new(VXLS_X_SIZE, VXLS_Y_SIZE, VXLS_Z_SIZE);
    dbg!("partie_1_2: before monte_carlo()");
    monte_carlo(700, &mut a);
    // dbg!(&a);
    println!("partie_1_2: before plot()");
    plot(&a)?;
    Ok(())
}

const fn cm(cm :f64) -> f64 {
    cm
}
const fn mm_to_vxl(mm: f64) -> i32 {
    (mm * VXL_PER_MM) as i32
}

struct SkinInVxl {
    stratum_corneum_lower_z: i32,
    epiderme_lower_z: i32,
    papille_dermique_lower_z: i32,
    derme_superieur_lower_z: i32,
    derme_reticulaire_lower_z: i32,
    derme_profond_lower_z: i32,
}
impl SkinInVxl {
    fn new(stratum_corneum_dz_in_mm: f64, epiderme_dz_in_mm: f64, papille_dermique_dz_in_mm: f64, derme_superieur_dz_in_mm:f64, derme_reticulaire_dz_in_mm:f64, derme_profond_dz_in_mm: f64) -> Self {
        
        let stratum_corneum_lower_z  =                             mm_to_vxl(stratum_corneum_dz_in_mm);
        let epiderme_lower_z         = stratum_corneum_lower_z   + mm_to_vxl(epiderme_dz_in_mm);
        let papille_dermique_lower_z = epiderme_lower_z          + mm_to_vxl(papille_dermique_dz_in_mm);
        let derme_superieur_lower_z  = papille_dermique_lower_z  + mm_to_vxl(derme_superieur_dz_in_mm);
        let derme_reticulaire_lower_z= derme_superieur_lower_z   + mm_to_vxl(derme_reticulaire_dz_in_mm);
        let derme_profond_lower_z    = derme_reticulaire_lower_z + mm_to_vxl(derme_profond_dz_in_mm);
        Self {
            stratum_corneum_lower_z,
            epiderme_lower_z,
            papille_dermique_lower_z,
            derme_superieur_lower_z,
            derme_reticulaire_lower_z,
            derme_profond_lower_z,
        }
    }
}

#[derive(Debug, EnumCount, Clone, Copy, PartialEq, Default)]
enum SkinLayerKind {
    #[default]
    StatumCorneum,
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
impl SkinLayer {
    fn mu_a(&self, wavelength: f64, mu_w: f64, mu_art: f64, mu_vei: f64, mu_a_baseline: f64) -> f64 {
        let Self { kind, v_w, v_b, .. } = self;
        MM_PER_VXL * match kind {
            StatumCorneum | Epiderme =>  
                V_MEL * (6.6e+10 * wavelength.powf(-3.33) / UNIT_SCALE)  +  v_w * mu_w  +  (1. - (V_MEL + v_w)) * mu_a_baseline,
            _ => {  //  Tout sous l'Epiderme
                let v_a = *v_b / 2.;     //  split 50-50 with v_v
                let v_v = v_a;          //  split 50-50 with v_a
                v_a * mu_art  +  v_v * mu_vei  +  v_w * mu_w  +  (1.- (v_a + v_v + v_w)) * mu_a_baseline
            },                    
        }
    }
    
}

const INDICE_REFRAC_AIR:f64 = 1.000293; // à T.P.N. UNused but for isometric solution.

const SKIN_LAYERS: [SkinLayer; SkinLayerKind::COUNT] = [
    SkinLayer {kind: StatumCorneum,     dz_in_mm: 0.02, indice_refrac: 1.42, v_b: 0.00, v_w: 0.05, },
    SkinLayer {kind: Epiderme,          dz_in_mm: 0.25, indice_refrac: 1.42, v_b: 0.00, v_w: 0.20, },
    SkinLayer {kind: PapileDermique,    dz_in_mm: 0.10, indice_refrac: 1.39, v_b: 0.04, v_w: 0.50, },
    SkinLayer {kind: DermeSuperieur,    dz_in_mm: 0.08, indice_refrac: 1.39, v_b: 0.30, v_w: 0.60, },
    SkinLayer {kind: DermeReticulaire,  dz_in_mm: 0.20, indice_refrac: 1.39, v_b: 0.04, v_w: 0.70, },
    SkinLayer {kind: DermeProfond,      dz_in_mm: 0.30, indice_refrac: 1.39, v_b: 0.01, v_w: 0.70, },
];

fn vxl_z_ranges() -> [Range<i32>; SkinLayerKind::COUNT] {
    let mut dz_in_vxl_acc = 0;
    SKIN_LAYERS.map(|skin_layer| { 
        let dz_in_vxl_start = dz_in_vxl_acc; 
        dz_in_vxl_acc += mm_to_vxl(skin_layer.dz_in_mm);
        dz_in_vxl_start..dz_in_vxl_acc
    })
}
fn indice_refrac_ratios() -> [f64; SkinLayerKind::COUNT] {
    let mut indice_refrac_denominateur = INDICE_REFRAC_AIR;
    SKIN_LAYERS.map(|skin_layer| { 
        let indice_refrac_numerateur = indice_refrac_denominateur /* precedent*/ ;
        indice_refrac_denominateur = skin_layer.indice_refrac;
        indice_refrac_numerateur / indice_refrac_denominateur
    })
}

static SKIN: OnceLock<Skin> = OnceLock::new();
#[derive(Debug, )]
struct Skin {
    layers: [SkinLayer; SkinLayerKind::COUNT],
    vxl_z_ranges: [Range<i32>; SkinLayerKind::COUNT],
    indice_refrac_ratios: [f64; SkinLayerKind::COUNT],
}

impl Skin {
    fn new() -> Self {

        Self {
            layers: SKIN_LAYERS,
            vxl_z_ranges : vxl_z_ranges(),
            indice_refrac_ratios : indice_refrac_ratios(),
        }
        
    }
    fn vxl_z_ranges() -> &'static [Range<i32>; SkinLayerKind::COUNT] {
        &SKIN.get_or_init(Skin::new).vxl_z_ranges
    }
    fn indice_refrac_ratios() -> &'static [f64; SkinLayerKind::COUNT] {
        &SKIN.get_or_init(Skin::new).indice_refrac_ratios
    }
    fn vxl_z_range(skin_layer_kind: SkinLayerKind) -> Range<i32> {
        Self::vxl_z_ranges()[skin_layer_kind as usize].clone()
    }
    fn indice_refrac_ratio(skin_layer_kind: SkinLayerKind) -> f64 {
        Self::indice_refrac_ratios()[skin_layer_kind as usize]
    }
    fn mu_t(mu_a:f64) -> f64 {
        mu_a + MU_S
    }
    fn skin_layer_kind(photon_pos_z: i32) -> Option<(i32, SkinLayerKind)> {
        SKIN_LAYERS.iter()
            .find_map(|skin_layer| {
                let vxl_z_range = Self::vxl_z_range(skin_layer.kind);
                if vxl_z_range.contains(&photon_pos_z) {
                    Some((vxl_z_range.start, skin_layer.kind))
                } else {None}
            })
    }

    fn absorption(photon: &Photon) -> f64 {
        let Photon { pos, path_seg, .. } = photon; 
     
        let src_pos_z = pos.z;
        if let Some((_, src_skin_layer_kind)) = Self::skin_layer_kind(pos.z) {
            
            let dst_photon_pos_z = src_pos_z + path_seg.dz();
            if let Some((dst_layer_first_vxl_z, dst_skin_layer_kind)) =  Self::skin_layer_kind(dst_photon_pos_z) {

                //  Trans skin layers with different indice_refrac.
                if src_skin_layer_kind != dst_skin_layer_kind {
                    let refrac_ratio = Self::indice_refrac_ratio(dst_skin_layer_kind);
                    if 1.0 != refrac_ratio {
                        debug!("src: {src_skin_layer_kind:?}, dst: {dst_skin_layer_kind:?} refrac_ratio: {refrac_ratio:?}, dst_layer_first_vxl_z: {dst_layer_first_vxl_z:?}");
                    }
                }
            }
        }
     
        0.
    }
}

#[derive(Debug, )]
struct PhotonProfile {
    wavelength: f64,
    mu_ts : [f64;  SkinLayerKind::COUNT],
    mu_a_on_mu_ts : [f64;  SkinLayerKind::COUNT],
    rng: ThreadRng,
}

impl PhotonProfile {
    fn try_new(wavelength_i:i32) -> Option<Self> {
        let wavelength = wavelength_i as f64;
        if let Some((val_hb_0, val_hb_1)) = hb_at_index(wavelength_i)  &&
           let Some(val_water) = water_hashmap(wavelength_i) {

            //  "J'ai assumé que l'output était en cm^-1, donc divisé par 10 pour mm^-1"
            let mu_a_baseline=  7.84e-7 * wavelength.powf(-3.255) / UNIT_SCALE;   //  reciprocal_mm
            
            let mu_absorption_oxyhemoglobine = val_hb_0 * MAGIG_NUM_0; // reciprocal_mm
            let mu_a_hemoglobine = val_hb_1 * MAGIG_NUM_0;  //  reciprocal_mm
        
            let mu_w = val_water / UNIT_SCALE;  //  reciprocal_mm

            //  Coeff artere & veine
            let (mu_art, mu_vei) = (
                SA_O2 * mu_absorption_oxyhemoglobine + (1. - SA_O2) * mu_a_hemoglobine,
                SV_O2 * mu_absorption_oxyhemoglobine + (1. - SV_O2) * mu_a_hemoglobine );
                   
            let mu_as = SKIN_LAYERS.map(|skin_layer| skin_layer.mu_a(wavelength, mu_w, mu_art, mu_vei, mu_a_baseline));
            
            Some(Self {
                wavelength,
                mu_ts: mu_as.map(|mu_a| Skin::mu_t(mu_a)),
                mu_a_on_mu_ts: mu_as.map(|mu_a| mu_a / Skin::mu_t(mu_a)),
                rng:  rand::rng(),
            })
        } else { None }
    }
    
    fn mu_t(&self, skin_layer_kind: SkinLayerKind) -> f64 {
        self.mu_ts[skin_layer_kind as usize]
    }
    
    fn new_path_seg_len(&mut self, src_skin_layer_kind: SkinLayerKind) -> f64 {
        //  parcours du groupe de photon pour cette itération
        let mut xsi_1:f64 = self.rng.random();
        while xsi_1 == 0. {             //  Exclure 0 de l'interval parce que ln(0) = ∞, ou NaN.
            xsi_1 = self.rng.random();
        }
        //  -xsi_1.ln() runs from 0. to very large (f64::MAX ?).
        -xsi_1.ln() / self.mu_t(src_skin_layer_kind)
    }

}

#[derive(Debug, Clone, Copy)]
struct VoxelPos {
    x: i32,
    y: i32,
    z: i32,
}
#[derive(Debug, Clone, Copy)]
struct VecUnit{
    ux: f64,
    uy: f64,
    uz: f64,
}
impl Default for VecUnit {
    fn default() -> Self {
        Self {ux: 0., uy: 0., uz: 1.,}
    }
}
#[derive(Debug, Clone, Copy, Default)]
struct PhotonPathSeg {
    len : f64,
    dir : VecUnit,
}

impl PhotonPathSeg {
    
    fn set_dir(&mut self, ux_uy_uz:&(f64, f64, f64)) {
        self.dir.ux = ux_uy_uz.0;
        self.dir.uy = ux_uy_uz.1;
        self.dir.uz = ux_uy_uz.2;
    }
    fn dz(&self) -> i32 {
        (self.len * self.dir.uz) as i32
    }
    fn generate_next(self, len: f64) -> Self {
        Self {
            len,
            dir: self.dir
        }
    }
}
#[derive(Debug)]
struct Photon<'a> {
    profile: &'a mut PhotonProfile,
    poids: f64,
    pos: VoxelPos,
    path_seg: PhotonPathSeg,
    skin_layer_kind: SkinLayerKind,
}

impl<'a> Photon<'a> {
    fn new(profile: &'a mut PhotonProfile) -> Self {
        Self {
            profile, 
            poids: 1.0,
            pos: VoxelPos { 
                x: VXLS_X_SIZE / 2,
                y: VXLS_Y_SIZE / 2,
                z: 0,
            },
            path_seg: PhotonPathSeg::default(),
            skin_layer_kind: SkinLayerKind::default(),
        }
    }
    
    fn generate_next_path_seg(&mut self) {
        self.path_seg.generate_next(self.profile.new_path_seg_len(self.skin_layer_kind));
    }
}

/// Return sin if cos and cos if sin
fn co_sin(sin_or_cos:f64) -> f64 {
    (1. - sin_or_cos.powi(2)).sqrt()
}


/// wavelength in nm 
fn monte_carlo1(wavelength_i:i32, _a:&mut VoxelsXZ) {

    if let Some(mut photon_profile) = PhotonProfile::try_new(wavelength_i) {
        for _n_th in 0..NB_PHOTONS_1_2 {
            let mut  photon = Photon::new(&mut photon_profile);
            dbg!(&photon);
            
            while SEUIL_DE_POIDS_CRITIQUE < photon.poids {
                photon.generate_next_path_seg();
                dbg!(&photon);
                break;
            }
            
        }
    }
}


fn plot(vxl:&VoxelsXZ) -> Result<()> {
    
    Plot::new()
        .heatmap(vxl, Some(HeatmapConfig::new()
            .colorbar(true)
            .colorbar_label("Énergie absorbée")
            //  No way to scale the colorbar axis in log yet
            ))
        .ylim(VXLS_Z_SIZE as f64 * CM_PER_VXL, 0.)    
        .xlim( 0., VXLS_X_SIZE as f64 * CM_PER_VXL)    
        .ylabel("Profondeur (mm)")
        .xlabel("Position en x(mm)")
        .save(".heatmap.png")?;

    Ok(())
}
// region:      --- Absorption funtions

// type Cm = uom::si::Quantity<(dyn uom::si::Dimension<I = uom::typenum::Z0, J = uom::typenum::Z0, Kind = (dyn uom::Kind + 'static), L = uom::typenum::PInt<uom::typenum::UInt<uom::typenum::UTerm, uom::typenum::B1>>, M = uom::typenum::Z0, N = uom::typenum::Z0, T = uom::typenum::Z0, Th = uom::typenum::Z0> + 'static), (dyn uom::si::Units<f64, amount_of_substance = uom::si::amount_of_substance::mole, electric_current = uom::si::electric_current::ampere, length = uom::si::length::meter, luminous_intensity = uom::si::luminous_intensity::candela, mass = uom::si::mass::kilogram, thermodynamic_temperature = uom::si::thermodynamic_temperature::kelvin, time = uom::si::time::second> + 'static), f64>;
// fn cm(cm :f64) -> Cm {
//     Length::new::<centimeter>(cm)
// }

/// Returns mu_a_baseline
fn absorption_baseline(wavelength:f64) -> f64 {
    //  "J'ai assumé que l'output était en cm^-1, donc divisé par 10 pour mm^-1"
    7.84e-7 * wavelength.powf(-3.255) / UNIT_SCALE
}

/// Returns (mu_art, mu_vei) le mu des artériel et veineux
fn coeff_artere_veine( mu_absorption_oxyhemoglobine: f64, mu_a_hemoglobine: f64) -> (f64, f64){
    (   SA_O2 * mu_absorption_oxyhemoglobine + (1. - SA_O2) * mu_a_hemoglobine,
        SV_O2 * mu_absorption_oxyhemoglobine + (1. - SV_O2) * mu_a_hemoglobine )
}

/// Return mu_a
fn absorption_couche(v_a:f64, mu_art:f64, v_v:f64, mu_vei:f64, v_w:f64, mu_w:f64, mu_a_baseline:f64) -> f64 {
    v_a * mu_art  +  v_v * mu_vei  +  v_w * mu_w  +  (1.- (v_a + v_v + v_w)) * mu_a_baseline
}

/// Return mu_a_epi
fn absorption_epiderme(wavelength:f64, v_w_epi:f64, mu_w_epi:f64, mu_a_baseline:f64) -> f64 {
            //  mu_a_mel : "J'ai assumé que l'output était en cm^-1, donc divisé par 10 pour mm^-1"
    V_MEL * (6.6e+10 * wavelength.powf(-3.33) / UNIT_SCALE)  +  v_w_epi * mu_w_epi  +  (1. - (V_MEL + v_w_epi)) * mu_a_baseline
}
// endregion:   --- Absorption functions


/// wavelength in nm 
fn monte_carlo(wavelength_i:i32, a:&mut VoxelsXZ) {

    let wavelength_f = wavelength_i as f64;
    let g_square = G.powi(2);
            
    if let Some((val_hb_0, val_hb_1)) = hb_at_index(wavelength_i)  &&
       let Some(val_water) = water_hashmap(wavelength_i) {

        let _dz = CM_PER_VXL; 
        let dz_scaling = MM_PER_VXL;

        let mu_a_baseline = absorption_baseline(wavelength_f); //  #mm^-1
        let mu_absorption_oxyhemoglobine = val_hb_0 * MAGIG_NUM_0; // mm^-1
        let mu_a_hemoglobine = val_hb_1 * MAGIG_NUM_0;  //  mm^-1
        let mu_w = val_water / UNIT_SCALE;  //  mm^-1
        let (mu_art, mu_vei) = coeff_artere_veine(mu_absorption_oxyhemoglobine, mu_a_hemoglobine);
        let mu_s = 20.1* dz_scaling; 
        
        let add_to = |prev_z: i32, increment: f64| prev_z + (increment as i32);
        let z_stratum = add_to(0, cm(0.02) / dz_scaling);
        let z_epiderme = add_to( z_stratum, cm(0.25) / dz_scaling);
        let z_papille = add_to( z_epiderme, cm(0.1) / dz_scaling);
        let z_derme_superieur = add_to( z_papille, cm(0.08) / dz_scaling);
        let z_derme_reticulaire = add_to( z_derme_superieur, cm(0.2) / dz_scaling);
        let z_derme_profond = add_to( z_derme_reticulaire, cm(0.2) / dz_scaling);


        let mut rng = rand::rng();  //  random number generator

        
        for n_th in 0..NB_PHOTONS_1_2 {

            let mut x = VXLS_X_SIZE / 2;
            let mut y = VXLS_Y_SIZE / 2;
            let mut z: i32 = 0;

            let hors_jeu =  |x, y, z|  !(0..VXLS_X_SIZE).contains(&x) || !(0..VXLS_Y_SIZE).contains(&y)   ||  !(0..VXLS_Z_SIZE).contains(&z);

            let mut poids = 1.;     //  poids du photon
            //  parcours du groupe de photon pour cette itération
            let mut prcr = PhotonPathSeg {   
                len: 0.,
                dir: VecUnit { ux: 0., uy: 0., uz: 1. }
            };


            while SEUIL_DE_POIDS_CRITIQUE < poids {

                let mut mu_a = if z <= z_stratum {
                    absorption_epiderme(wavelength_f, 0.05, mu_w, mu_a_baseline)
                } else if z_stratum < z && z <= z_epiderme {
                    absorption_epiderme(wavelength_f, 0.20, mu_w, mu_a_baseline)
                } else if z_epiderme < z && z <= z_papille {
                    let v_b = 0.04;
                    let half_v_b = v_b / 2.;
                    absorption_couche(half_v_b, mu_art, half_v_b, mu_vei, 0.5, mu_w, mu_a_baseline)                           
                } else if z_papille < z && z <= z_derme_superieur {
                    let v_b = 0.3;
                    let half_v_b = v_b / 2.;
                    absorption_couche(half_v_b, mu_art, half_v_b, mu_vei, 0.6, mu_w, mu_a_baseline)                           
                } else if z_derme_superieur < z && z <= z_derme_reticulaire {
                    let v_b = 0.04;
                    let half_v_b = v_b / 2.;
                    absorption_couche(half_v_b, mu_art, half_v_b, mu_vei, 0.7, mu_w, mu_a_baseline)                           
                } else if z_derme_reticulaire < z && z <= z_derme_profond {
                    let v_b = 0.1;
                    let half_v_b = v_b / 2.;
                    absorption_couche(half_v_b, mu_art, half_v_b, mu_vei, 0.7, mu_w, mu_a_baseline)                           
                } else {
                    break;
                };

                mu_a *= dz_scaling;     //  Conversion de mm^-1 à voxels^-1
                let mu_t = mu_a + mu_s;

                //  parcours du groupe de photon pour cette itération
                let mut xsi_1:f64 = rng.random();
                while xsi_1 == 0. {             //  Exclure 0 de l'interval parce que ln(0) = ∞, ou NaN.
                    xsi_1 = rng.random();
                }
                //  Distance de parcours
                prcr.len = -xsi_1.ln() / mu_t;    //  toujours positif

                //  On est à la frontière (entre l'épiderme et la papille dermique) et on doit vérifier s'il y a réflexion ou transmission
                if z <= z_epiderme  &&  z_epiderme < z + (prcr.dir.uz * prcr.len) as i32 {

                    let distance_frontiere = (z_epiderme - z) as f64 / prcr.dir.uz;
                    let distance_restante = prcr.len - distance_frontiere;

                    //  le photon est rendu à la frontière 
                    z = z_epiderme;
                    x += (prcr.dir.ux * distance_frontiere) as i32;  
                    y += (prcr.dir.uy * distance_frontiere) as i32;

                    if hors_jeu(x, y, z) {
                        break;
                    }

                    let mut cos_theta_i = prcr.dir.uz.abs();
                    if 1.0 < cos_theta_i {                              //  cos_theta_i in [0., 1.] range
                        cos_theta_i = 1.0; 
                    }                                                   
                    let sin_theta_i = (1. - cos_theta_i.powi(2)).sqrt();   // sin_theta_i in [0., 1.] range too.
                    let sin_theta_t = N1_SUR_N2 * sin_theta_i;  //  Loi de Snell-Descartes

                    //  Calcul de réflexion
                    let (r, cos_theta_t)  = if sin_theta_t >= 1.0 {     //  Réflexion totale interne
                        (1.0, f64::NAN)     //  100% réflexion, cos_theta_t non utilisé plus bas
                    } else {
                        let cos_theta_t = co_sin(sin_theta_t);
                        ((N1 * cos_theta_i  -  N2 * cos_theta_t) / (N1 * cos_theta_i  +  N2 * cos_theta_t), cos_theta_t)
                    };

                    let r_square = r.powi(2);

                    let delta_poids = if rng.random::<f64>() < r {   //  xsi_5 <  r : Reflexion
                        //  Seul dir.uz est réfléchi les autres directions unitaires ne changent pas.
                        prcr.dir.uz = - prcr.dir.uz;

                        (1. - r_square) * poids
                    } else {                        //  xsi_5 >= r : Transmission
                        // Calcul de la nouvelle direction après réfraction
                        prcr.dir.ux *=  N1_SUR_N2;
                        prcr.dir.uy *=  N1_SUR_N2;
                        prcr.dir.uz = cos_theta_t;

                            r_square * poids
                    };
                    poids -= delta_poids;
                    a.vxls[x as usize][z as usize] += delta_poids;

                    x += (prcr.dir.ux * distance_restante) as i32;
                    y += (prcr.dir.uy * distance_restante) as i32;
                    z += (prcr.dir.uz * distance_restante) as i32;

                    if hors_jeu(x, y, z) {
                        break;
                    }
                } 
                //  Déplacement normal, pas de frontières et réflexion de Fresnel
                else {
                    x += (prcr.dir.ux * prcr.len) as i32;
                    y += (prcr.dir.uy * prcr.len) as i32;
                    z += (prcr.dir.uz * prcr.len) as i32;

                    if hors_jeu(x, y, z) {
                        break;
                    }

                    //  Absorption
                    let delta_poids = mu_a / mu_t * poids;
                    poids -= delta_poids;
                    a.vxls[x as usize][z as usize] += delta_poids;

                }

                let phi = 2.0 * PI * rng.random::<f64>();                                          //  xsi_2
                let cos_phi = phi.cos();
                let sin_phi = phi.sin();
                                                                                                        //  xsi_3
                let mut cos_theta = 1. / (2.* G) * (1. + g_square - ((1. - g_square)/(1. - G + 2. * G * rng.random::<f64>())).powi(2)); 
                //  Encore une fois, pour éviter les erreurs arithmétiques flottantes (cos doit être entre -1, 1)
                cos_theta = cos_theta.clamp(-1., 1.);
                let sin_theta = (1. - cos_theta.powi(2)).sqrt();
                // let theta = sin_theta.asin();

                prcr.set_dir(&( if prcr.dir.uz.abs() == 1. {  // On regarde si mu_z != 1
                    (   sin_theta * sin_phi,    //  prcr.dir.ux
                        sin_theta * cos_phi,    //  prcr.dir.uy
                        cos_theta,              //  prcr.dir.uz
                    )
                } else {
                    let not_drctn_uz = (1. - prcr.dir.uz.powi(2)).sqrt();
                    let drctn_uz_cos_phi = prcr.dir.uz * cos_phi;

                    (   sin_theta * (prcr.dir.ux * drctn_uz_cos_phi  -  prcr.dir.uy * sin_phi) / not_drctn_uz  +  cos_theta * prcr.dir.ux,
                        sin_theta * (prcr.dir.uy * drctn_uz_cos_phi  -  prcr.dir.ux * sin_phi) / not_drctn_uz  +  cos_theta * prcr.dir.uy,
                        prcr.dir.uz * cos_theta  -  not_drctn_uz * sin_theta * cos_phi
                    )
                }));

                //  Si le poids du paquet de photon soit en dessous d’un seuil critique (Wc) prédéterminé, on joue
                //  à la roulette: Le paquet a une chance sur M de continuer à être simulé, sinon "il disparaı̂t". 
                //  Pour éviter que de l’énergie disparaisse pris dans l'ensemble, lorsque le photon survie, 
                //  on multiplie son énergie par M afin de compenser pour les fois ou l’énergie disparait.
                if poids < SEUIL_DE_POIDS_CRITIQUE  &&  rng.random::<f64>() < M_INV {   //  Roulette !  (xsi_4)
                    poids *= M; 
                }
            }
            if n_th % 5 == 0 {
                println!("{n_th}th photon : a[x:{x}][z:{z}]{:?}", a[x as usize][z as usize]);
            }

        }
    }
}
