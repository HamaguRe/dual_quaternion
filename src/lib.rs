extern crate num_traits;
extern crate dual_number;
extern crate quaternion;

use num_traits::float::Float;
use dual_number::DualNumber;
use quaternion as quat;
use quaternion::{Quaternion, Vector3};

pub type DualQuaternion<T> = (Quaternion<T>, Quaternion<T>);  // (primary part, dual part)


#[inline(always)]
pub fn id<T>() -> DualQuaternion<T> 
where T: Float {
    let zero = T::zero();

    ( quat::id(), (zero, [zero; 3]) )
}

/// aの回転とbの移動を表すデュアルクォータニオンを生成
#[inline(always)]
pub fn from_quat_and_vec<T>(a: Quaternion<T>, b: Vector3<T>) -> DualQuaternion<T> 
where T: Float {
    let zero = T::zero();
    let one  = T::one();
    let two  = one + one;

    let dual = quat::mul( (zero, b), a );
    let dual = quat::mul_scalar_quat(one / two, dual);
    (a, dual)
}

/// 位置ベクトルの回転・移動を行う
#[inline(always)]
pub fn vector_translation<T>(a: DualQuaternion<T>, r: Vector3<T>) -> Vector3<T> 
where T: Float {
    let zero = T::zero();
    let one  = T::one();

    let r_dq = ( (one, [zero; 3]), (zero, r) );
    let a_conj_dq = conj(a);
    let a_conj_dq_dn = conj_dual_num(a_conj_dq);
    let result = mul( a, mul(r_dq, a_conj_dq_dn) );
    (result.1).1
}

/// 座標系の回転・移動を行う
#[inline(always)]
pub fn coordinate_translation<T>(a: DualQuaternion<T>, r: Vector3<T>) -> Vector3<T> 
where T: Float {
    let zero = T::zero();
    let one  = T::one();

    let r_dq = ( (one, [zero; 3]), (zero, r) );
    let a_conj_dq = conj(a);
    let a_conj_dq_dn = conj_dual_num(a_conj_dq);
    let result = mul( a_conj_dq_dn, mul(r_dq, a) );
    (result.1).1
}

/// 角速度と加速度からdt間のデュアルクォータニオンの変化量を求める．
#[inline(always)]
pub fn integration<T>(omega: Vector3<T>, accel: Vector3<T>, dt: T) -> DualQuaternion<T> 
where T: Float {
    let dq = quat::integration(omega, dt);
    let dr = quat::mul_scalar_vec(dt * dt, accel);  // a[m/s]^2 * t[s]^2
    from_quat_and_vec(dq, dr)
}

/// ボディ角速度を用いた場合の積分
#[inline(always)]
pub fn vector_integration<T>(d: DualQuaternion<T>, omega: Vector3<T>, accel: Vector3<T>, dt: T) -> DualQuaternion<T> 
where T: Float {
    let ddq = integration(omega, accel, dt);
    mul(ddq, d)
}

/// 空間角速度を使った積分
#[inline(always)]
pub fn coordinate_integration<T>(d: DualQuaternion<T>, omega: Vector3<T>, accel: Vector3<T>, dt: T) -> DualQuaternion<T> 
where T: Float {
    let ddq = integration(omega, accel, dt);
    mul(d, ddq)
}

#[inline(always)]
pub fn get_translation<T>(a: DualQuaternion<T>) -> Vector3<T> 
where T: Float {
    let two = T::one() + T::one();

    let primary_conj = quat::conj(a.0);
    let r = quat::mul(primary_conj, a.1);
    let r = quat::mul_scalar_quat(two, r);
    r.1
}

#[inline(always)]
pub fn get_rotation<T>(a: DualQuaternion<T>) -> Quaternion<T> {
    a.0
}

#[inline(always)]
pub fn add<T>(a: DualQuaternion<T>, b: DualQuaternion<T>) -> DualQuaternion<T> 
where T: Float {
    let primary = quat::add(a.0, b.0);
    let dual    = quat::add(a.1, b.1);
    (primary, dual)
}

/// Add "DualNumber" and "DualQuaternion"
#[inline(always)]
fn add_num_quat<T>(a: DualNumber<T>, b: DualQuaternion<T>) -> DualQuaternion<T> 
where T: Float {
    let zero = T::zero();

    let primary = quat::add( (a.0, [zero; 3]), b.0 );
    let dual    = quat::add( (a.1, [zero; 3]), b.1 );
    (primary, dual)
}

#[inline(always)]
pub fn mul<T>(a: DualQuaternion<T>, b: DualQuaternion<T>) -> DualQuaternion<T> 
where T: Float {
    let primary = quat::mul(a.0, b.0);
    let dual_0  = quat::mul(a.0, b.1);
    let dual_1  = quat::mul(a.1, b.0);
    let dual = quat::add(dual_0, dual_1);
    (primary, dual)
}

/// Multiplication of "Dual Number" and "Dual Quaternion"
#[inline(always)]
fn mul_num_quat<T>(a: DualNumber<T>, b: DualQuaternion<T>) -> DualQuaternion<T> 
where T: Float {
    let primary = quat::mul_scalar_quat(a.0, b.0);
    let dual_0 = quat::mul_scalar_quat(a.0, b.1);
    let dual_1 = quat::mul_scalar_quat(a.1, b.0);
    let dual = quat::add(dual_0, dual_1);
    (primary, dual)
}

/// Conjugate of DualQuaternion
/// (q0 + εq1)  -->  (q0* + εq1*)
#[inline(always)]
pub fn conj<T>(a: DualQuaternion<T>) -> DualQuaternion<T> 
where T: Float {
    let primary = quat::conj(a.0);
    let dual = quat::conj(a.1);
    (primary, dual)
}

/// Conjugate of DualNumber
/// (q0 + εq1)  -->  (q0 - εq1)
#[inline(always)]
pub fn conj_dual_num<T>(a: DualQuaternion<T>) -> DualQuaternion<T> 
where T: Float {
    let one = T::one();

    let dual = quat::mul_scalar_quat(-one, a.1);
    (a.0, dual)
}

#[inline(always)]
pub fn inverse<T>(a: DualQuaternion<T>) -> DualQuaternion<T> 
where T: Float {
    let one = T::one();

    let a_conj = conj(a);
    let a_norm = norm(a);
    let b = dual_number::mul(a_norm, a_norm);
    let primary = quat::mul_scalar_quat(one / b.0, a_conj.0);
    let dual_0 = quat::mul_scalar_quat(one / b.0, a_conj.1);
    let dual_1 = quat::mul_scalar_quat(-b.1 / (b.0 * b.0), a_conj.0);
    let dual = quat::add(dual_0, dual_1);
    (primary, dual)
}


// ---------- これ以降実装が合っているか怪しい ---------- //

#[inline(always)]
pub fn norm<T>(a: DualQuaternion<T>) -> DualNumber<T> 
where T: Float {
    let primary_norm = quat::norm(a.0);
    let dual = quat::dot(a.0, a.1) / primary_norm;
    (primary_norm, dual)
}

#[inline(always)]
pub fn normalize<T>(a: DualQuaternion<T>) -> DualQuaternion<T> 
where T: Float {
    let tmp = dual_number::inverse( norm(a) );
    mul_num_quat(tmp, a)
}

#[inline(always)]
pub fn dot<T>(a: DualQuaternion<T>, b: DualQuaternion<T>) -> T 
where T: Float {
    let primary = quat::dot(a.0, b.0);
    let dual = quat::dot(a.1, b.1);
    dual + primary
}

#[inline(always)]
pub fn exp<T>(a: DualQuaternion<T>) -> DualQuaternion<T> 
where T: Float {
    let s = normalize(a);
    let arg = norm(a);
    let dual_num = dual_number::cos(arg);
    let tmp = dual_number::sin(arg);
    let dual_quat = mul_num_quat(tmp, s);
    add_num_quat(dual_num, dual_quat)
}
