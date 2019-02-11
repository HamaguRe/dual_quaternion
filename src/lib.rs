extern crate dual_number;
extern crate quaternion;

pub type Vector3 = [f64; 3];
pub type Quaternion = (f64, Vector3);
pub type DualNumber = (f64, f64);  // (real part, dual part)
pub type DualQuaternion = (Quaternion, Quaternion);  // (primary part, dual part)


#[inline(always)]
pub fn id() -> DualQuaternion {
    ( (1.0, [0.0; 3]), (0.0, [0.0; 3]) )
}

/// aの回転とbの移動を表すデュアルクォータニオンを生成
#[inline(always)]
pub fn from_quat_and_vec(a: Quaternion, b: Vector3) -> DualQuaternion {
    let dual = quaternion::mul( (0.0, b), a );
    let dual = quaternion::mul_scalar_quat(1.0 / 2.0, dual);
    (a, dual)
}

#[inline(always)]
pub fn vector_translation(a: DualQuaternion, r: Vector3) -> Vector3 {
    let r_dq: DualQuaternion = ( (1.0, [0.0; 3]), (0.0, r) );
    let a_conj_dq = conj(a);
    let a_conj_dq_dn = conj_dual_num(a_conj_dq);
    let result = mul( a, mul(r_dq, a_conj_dq_dn) );
    (result.1).1
}

#[inline(always)]
pub fn get_vec(a: DualQuaternion) -> Vector3 {
    let primary_conj = quaternion::conj(a.0);
    let r = quaternion::mul(primary_conj, a.1);
    let r = quaternion::mul_scalar_quat(2.0, r);
    r.1
}

#[inline(always)]
pub fn add(a: DualQuaternion, b: DualQuaternion) -> DualQuaternion {
    let primary = quaternion::add(a.0, b.0);
    let dual    = quaternion::add(a.1, b.1);
    (primary, dual)
}

/// Add "DualNumber" and "DualQuaternion"
#[inline(always)]
pub fn add_num_quat(a: DualNumber, b: DualQuaternion) -> DualQuaternion {
    let primary = quaternion::add( (a.0, [0.0; 3]), b.0 );
    let dual    = quaternion::add( (a.1, [0.0; 3]), b.1 );
    (primary, dual)
}

#[inline(always)]
pub fn norm(a: DualQuaternion) -> DualNumber {
    let primary_norm = quaternion::norm(a.0);
    let dual = quaternion::dot(a.0, a.1) / primary_norm;
    (primary_norm, dual)
}

#[inline(always)]
pub fn normalize(a: DualQuaternion) -> DualQuaternion {
    let tmp = dual_number::inverse( norm(a) );
    mul_num_quat(tmp, a)
}

// あってるかわからん
#[inline(always)]
pub fn dot(a: DualQuaternion, b: DualQuaternion) -> f64 {
    let primary = quaternion::dot(a.0, b.0);
    let dual = quaternion::dot(a.1, b.1);
    dual + primary
}

#[inline(always)]
pub fn mul(a: DualQuaternion, b: DualQuaternion) -> DualQuaternion {
    let primary = quaternion::mul(a.0, b.0);
    let dual_0  = quaternion::mul(a.0, b.1);
    let dual_1  = quaternion::mul(a.1, b.0);
    let dual = quaternion::add(dual_0, dual_1);
    (primary, dual)
}

/// Multiplication of "Dual Number" and "Dual Quaternion"
#[inline(always)]
fn mul_num_quat(a: DualNumber, b: DualQuaternion) -> DualQuaternion {
    let primary = quaternion::mul_scalar_quat(a.0, b.0);
    let dual_0 = quaternion::mul_scalar_quat(a.0, b.1);
    let dual_1 = quaternion::mul_scalar_quat(a.1, b.0);
    let dual = quaternion::add(dual_0, dual_1);
    (primary, dual)
}

/// Conjugate of DualQuaternion
/// (q0 + εq1)  -->  (q0* + εq1*)
#[inline(always)]
pub fn conj(a: DualQuaternion) -> DualQuaternion {
    let primary = quaternion::conj(a.0);
    let dual = quaternion::conj(a.1);
    (primary, dual)
}

/// Conjugate of DualNumber
/// (q0 + εq1)  -->  (q0 - εq1)
#[inline(always)]
pub fn conj_dual_num(a: DualQuaternion) -> DualQuaternion {
    let dual = quaternion::mul_scalar_quat(-1.0, a.1);
    (a.0, dual)
}

#[inline(always)]
pub fn inverse(a: DualQuaternion) -> DualQuaternion {
    let a_conj = conj(a);
    let a_norm = norm(a);
    let b = dual_number::mul(a_norm, a_norm);
    let primary = quaternion::mul_scalar_quat(1.0 / b.0, a_conj.0);
    let dual_0 = quaternion::mul_scalar_quat(1.0 / b.0, a_conj.1);
    let dual_1 = quaternion::mul_scalar_quat(-b.1 / (b.0 * b.0), a_conj.0);
    let dual = quaternion::add(dual_0, dual_1);
    (primary, dual)
}

// あってるか分からない
#[inline(always)]
pub fn exp(a: DualQuaternion) -> DualQuaternion {
    let s = normalize(a);
    let arg = norm(a);
    let dual_num = dual_number::cos(arg);
    let tmp = dual_number::sin(arg);
    let dual_quat = mul_num_quat(tmp, s);
    add_num_quat(dual_num, dual_quat)
}
