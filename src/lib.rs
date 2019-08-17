// Dual-Quaternion Libraly
// f64 only

extern crate dual_number;
extern crate quaternion as quat;
use dual_number::DualNumber;
use quat::{Quaternion, Vector3};

// (primary part, dual part)
pub type DualQuaternion<T> = (Quaternion<T>, Quaternion<T>);
// Identity dual quaternion
pub const IDENTITY: DualQuaternion<f64> = ( quat::IDENTITY, (0.0, [0.0; 3]) );


/// qの回転とrの移動を表す二重四元数を生成
#[inline(always)]
pub fn from_quat_vector(q: Quaternion<f64>, r: Vector3<f64>) -> DualQuaternion<f64> {
    let tmp  = quat::mul( (0.0, r), q );
    let dual = quat::scale(0.5, tmp);
    (q, dual)
}

/// 回転・移動を表す二重四元数から移動量を取り出す
#[inline(always)]
pub fn get_translation(a: DualQuaternion<f64>) -> Vector3<f64> {
    let r = quat::scale( 2.0, quat::mul( quat::conj(a.0), a.1 ) );
    r.1
}

#[inline(always)]
pub fn add(a: DualQuaternion<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = quat::add(a.0, b.0);
    let dual = quat::add(a.1, b.1);
    (prim, dual)
}

/// Add "DualNumber" and "DualQuaternion"
#[inline(always)]
pub fn add_num_quat(a: DualNumber<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = ( a.0 + primary(b).0, primary(b).1 );
    let dual = ( a.1 + dual(b).0, dual(b).1 );
    (prim, dual)
}

#[inline(always)]
pub fn sub(a: DualQuaternion<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = quat::sub(a.0, b.0);
    let dual = quat::sub(a.1, b.1);
    (prim, dual)
}

#[inline(always)]
pub fn scale(s: f64, q: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = quat::scale(s, q.0);
    let dual = quat::scale(s, q.1);
    (prim, dual)
}

#[inline(always)]
pub fn mul(a: DualQuaternion<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = quat::mul(a.0, b.0);
    let tmp0 = quat::mul(a.0, b.1);
    let tmp1 = quat::mul(a.1, b.0);
    let dual = quat::add(tmp0, tmp1);
    (prim, dual)
}

/// Multiplication of "Dual Number" and "Dual Quaternion"
/// 二重数と二重四元数の積は可換．
#[inline(always)]
pub fn mul_num_quat(a: DualNumber<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = quat::scale(a.0, b.0);
    let tmp0 = quat::scale(a.0, b.1);
    let tmp1 = quat::scale(a.1, b.0);
    let dual = quat::add(tmp0, tmp1);
    (prim, dual)
}

/// Conjugate of DualQuaternion
/// (q0 + εq1)  -->  (q0* + εq1*)
#[inline(always)]
pub fn conj(a: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = quat::conj(a.0);
    let dual = quat::conj(a.1);
    (prim, dual)
}

/// Conjugate of DualNumber
/// (q0 + εq1)  -->  (q0 - εq1)
#[inline(always)]
pub fn conj_dual_num(a: DualQuaternion<f64>) -> DualQuaternion<f64> {
    ( a.0, quat::negate(a.1) )
}

#[inline(always)]
pub fn norm(a: DualQuaternion<f64>) -> DualNumber<f64> {
    let prim_norm = quat::norm(a.0);
    let dual = quat::dot(a.0, a.1) / prim_norm;
    (prim_norm, dual)
}

#[inline(always)]
pub fn inverse(a: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let norm = norm(a);
    let norm_square = dual_number::mul(norm, norm);
    let prim = quat::scale(norm_square.0, a.0);
    let tmp0 = quat::scale(norm_square.0, a.1);
    let tmp1 = quat::scale(norm_square.1, a.0);
    let dual = quat::sub(tmp0, tmp1);
    scale( (norm_square.0 * norm_square.0).recip(), (prim, dual) )
}

/// 位置ベクトルの回転・移動を行う
#[inline(always)]
pub fn vector_translation(a: DualQuaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let tmp0 = quat::scale_vec( primary(a).0, dual(a).1 );
    let tmp1 = quat::scale_vec( dual(a).0, primary(a).1 );
    let tmp2 = quat::cross_vec( primary(a).1, dual(a).1 );
    let term1 = quat::add_vec( quat::sub_vec(tmp0, tmp1), tmp2 );
    let term2 = quat::vector_rotation( primary(a), r );
    quat::scale_add_vec(2.0, term1, term2)
}

/// 座標系の回転・移動を行う（位置ベクトルの逆）
#[inline(always)]
pub fn frame_translation(a: DualQuaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let tmp0 = quat::scale_vec( primary(a).0, dual(a).1 );
    let tmp1 = quat::scale_vec( dual(a).0, primary(a).1 );
    let tmp2 = quat::cross_vec( dual(a).1, primary(a).1 );
    let term1 = quat::add_vec( quat::sub_vec(tmp0, tmp1), tmp2 );
    let term2 = quat::frame_rotation( primary(a), r );
    quat::scale_add_vec(2.0, term1, term2)
}

/// 二重四元数のPrimary partを取り出す．
#[inline(always)]
fn primary(d: DualQuaternion<f64>) -> Quaternion<f64> {
    d.0
}

/// 二重四元数のDual partを取り出す．
#[inline(always)]
fn dual(d: DualQuaternion<f64>) -> Quaternion<f64> {
    d.1
}