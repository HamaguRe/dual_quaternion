// Dual-Quaternion Libraly
// f64 only

extern crate dual_number;
extern crate quaternion;

use dual_number::DualNumber;
use quaternion as quat;
use quaternion::{Quaternion, Vector3};

// (primary part, dual part)
pub type DualQuaternion<T> = (Quaternion<T>, Quaternion<T>);


#[inline(always)]
pub fn id() -> DualQuaternion<f64> {
    ( quat::id(), (0.0, [0.0; 3]) )
}

/// aの回転とbの移動を表すデュアルクォータニオンを生成
#[inline(always)]
pub fn from_quat_and_vec(a: Quaternion<f64>, b: Vector3<f64>) -> DualQuaternion<f64> {
    let dual = quat::mul( (0.0, b), a );
    let dual = quat::scale(0.5, dual);
    (a, dual)
}

/// 位置ベクトルの回転・移動を行う
#[inline(always)]
pub fn vector_translation(a: DualQuaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let r_dq = ( quat::id(), (0.0, r) );
    let a_conj_dq = conj(a);
    let a_conj_dq_dn = conj_dual_num(a_conj_dq);
    let result = mul( a, mul(r_dq, a_conj_dq_dn) );
    (result.1).1
}

/// 座標系の回転・移動を行う
#[inline(always)]
pub fn coordinate_translation(a: DualQuaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let r_dq = ( quat::id(), (0.0, r) );
    let a_conj_dq = conj(a);
    let a_conj_dq_dn = conj_dual_num(a_conj_dq);
    let result = mul( a_conj_dq_dn, mul(r_dq, a) );
    (result.1).1
}

/// デュアルクォータニオンから移動量を取り出す
#[inline(always)]
pub fn get_translation(a: DualQuaternion<f64>) -> Vector3<f64> {
    let primary_conj = quat::conj(a.0);
    let r = quat::mul(primary_conj, a.1);
    let r = quat::scale(2.0, r);
    r.1
}

#[inline(always)]
pub fn get_rotation(a: DualQuaternion<f64>) -> Quaternion<f64> {
    a.0
}

#[inline(always)]
pub fn add(a: DualQuaternion<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let primary = quat::add(a.0, b.0);
    let dual    = quat::add(a.1, b.1);
    (primary, dual)
}

/// Add "DualNumber" and "DualQuaternion"
#[inline(always)]
fn add_num_quat(a: DualNumber<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let primary = quat::add( (a.0, [0.0; 3]), b.0 );
    let dual    = quat::add( (a.1, [0.0; 3]), b.1 );
    (primary, dual)
}

#[inline(always)]
pub fn mul(a: DualQuaternion<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let primary = quat::mul(a.0, b.0);
    let dual_0  = quat::mul(a.0, b.1);
    let dual_1  = quat::mul(a.1, b.0);
    let dual = quat::add(dual_0, dual_1);
    (primary, dual)
}

/// Multiplication of "Dual Number" and "Dual Quaternion"
#[inline(always)]
fn mul_num_quat(a: DualNumber<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let primary = quat::scale(a.0, b.0);
    let dual_0 = quat::scale(a.0, b.1);
    let dual_1 = quat::scale(a.1, b.0);
    let dual = quat::add(dual_0, dual_1);
    (primary, dual)
}

/// Conjugate of DualQuaternion
/// (q0 + εq1)  -->  (q0* + εq1*)
#[inline(always)]
pub fn conj(a: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let primary = quat::conj(a.0);
    let dual = quat::conj(a.1);
    (primary, dual)
}

/// Conjugate of DualNumber
/// (q0 + εq1)  -->  (q0 - εq1)
#[inline(always)]
pub fn conj_dual_num(a: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let dual = quat::negate(a.1);
    (a.0, dual)
}

#[inline(always)]
pub fn inverse(a: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let a_conj = conj(a);
    let a_norm = norm(a);
    let b = dual_number::mul(a_norm, a_norm);
    let tmp = b.0.recip();
    let primary = quat::scale( tmp, a_conj.0 );
    let dual_0 =  quat::scale( tmp, a_conj.1 );
    let dual_1 = quat::scale(-b.1 / (b.0 * b.0), a_conj.0);
    let dual = quat::add(dual_0, dual_1);
    (primary, dual)
}


// ---------- これ以降実装が合っているか怪しい ---------- //

#[inline(always)]
pub fn norm(a: DualQuaternion<f64>) -> DualNumber<f64> {
    let primary_norm = quat::norm(a.0);
    let dual = quat::dot(a.0, a.1) / primary_norm;
    (primary_norm, dual)
}

#[inline(always)]
pub fn normalize(a: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let tmp = dual_number::inverse( norm(a) );
    mul_num_quat(tmp, a)
}

#[inline(always)]
pub fn dot(a: DualQuaternion<f64>, b: DualQuaternion<f64>) -> f64 {
    let primary = quat::dot(a.0, b.0);
    let dual = quat::dot(a.1, b.1);
    dual + primary
}

#[inline(always)]
pub fn exp(a: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let s = normalize(a);
    let arg = norm(a);
    let dual_num = dual_number::cos(arg);
    let tmp = dual_number::sin(arg);
    let dual_quat = mul_num_quat(tmp, s);
    add_num_quat(dual_num, dual_quat)
}
