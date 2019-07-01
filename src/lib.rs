// Dual-Quaternion Libraly
// f64 only

extern crate dual_number;
extern crate quaternion as quat;
use dual_number::DualNumber;
use quat::{Quaternion, Vector3};

// (primary part, dual part)
pub type DualQuaternion<T> = (Quaternion<T>, Quaternion<T>);


#[inline(always)]
pub fn id() -> DualQuaternion<f64> {
    ( quat::id(), (0.0, [0.0; 3]) )
}

/// aの回転とbの移動を表すデュアルクォータニオンを生成
#[inline(always)]
pub fn from_quat_vector(a: Quaternion<f64>, b: Vector3<f64>) -> DualQuaternion<f64> {
    let tmp  = quat::mul( (0.0, b), a );
    let dual = quat::scale(0.5, tmp);
    (a, dual)
}

/// 位置ベクトルの回転・移動を行う
#[inline(always)]
pub fn vector_translation(a: DualQuaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let term1 = quat_mul_only_vec( a.1, quat::conj(a.0) );
    let term2 = quat_mul_only_vec( a.0, quat::conj(a.1) );
    let term3 = quat::vector_rotation(a.0, r);
    quat::add_vec( quat::sub_vec(term1, term2), term3 )
}

/// 座標系の回転・移動を行う（位置ベクトルの逆）
#[inline(always)]
pub fn coordinate_translation(a: DualQuaternion<f64>, r: Vector3<f64>) -> Vector3<f64> {
    let term1 = quat_mul_only_vec( quat::conj(a.0), a.1 );
    let term2 = quat_mul_only_vec( quat::conj(a.1), a.0 );
    let term3 = quat::coordinate_rotation(a.0, r);
    quat::add_vec( quat::sub_vec(term1, term2), term3 )
}

/// デュアルクォータニオンから移動量を取り出す
#[inline(always)]
pub fn get_translation(a: DualQuaternion<f64>) -> Vector3<f64> {
    let tmp = quat::mul( quat::conj(a.0), a.1 );
    let r = quat::scale(2.0, tmp);
    r.1
}

#[inline(always)]
pub fn add(a: DualQuaternion<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = quat::add(a.0, b.0);
    let dual = quat::add(a.1, b.1);
    (prim, dual)
}

#[inline(always)]
pub fn sub(a: DualQuaternion<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = quat::sub(a.0, b.0);
    let dual = quat::sub(a.1, b.1);
    (prim, dual)
}

/// Add "DualNumber" and "DualQuaternion"
#[inline(always)]
fn add_num_quat(a: DualNumber<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let prim = ( a.0 + (b.0).0, (b.0).1 );
    let dual = ( a.1 + (b.1).0, (b.1).1 );
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
#[inline(always)]
fn mul_num_quat(a: DualNumber<f64>, b: DualQuaternion<f64>) -> DualQuaternion<f64> {
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
    let denom = norm.0 * norm.0;
    let prim = quat::scale(norm.0, a.0);
    let tmp0 = quat::scale(norm.0, a.1);
    let tmp1 = quat::scale(norm.1, a.0);
    let dual = quat::sub(tmp0, tmp1);
    scale( denom.recip(), (prim, dual) )
}


// ---------- これ以降実装が合っているか怪しい ---------- //
#[inline(always)]
pub fn normalize(a: DualQuaternion<f64>) -> DualQuaternion<f64> {
    let tmp = dual_number::inverse( norm(a) );
    mul_num_quat(tmp, a)
}

#[inline(always)]
pub fn dot(a: DualQuaternion<f64>, b: DualQuaternion<f64>) -> f64 {
    let prim = quat::dot(a.0, b.0);
    let dual = quat::dot(a.1, b.1);
    dual + prim
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

/// ハミルトン積の結果でベクトル部のみ返す．
/// ベクトル移動などで使用．
#[inline(always)]
fn quat_mul_only_vec(a: Quaternion<f64>, b: Quaternion<f64>) -> Vector3<f64> {
    let tmp = quat::scale_add_vec( b.0, a.1, quat::cross_vec(a.1, b.1) );
    quat::scale_add_vec( a.0, b.1, tmp )
}