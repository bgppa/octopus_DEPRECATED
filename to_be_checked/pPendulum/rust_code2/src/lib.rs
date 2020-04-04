#![feature(clamp)]
use std::cmp::Ordering::Equal;

use nalgebra as na;

mod acrobot;

/// This function get the parameters of a neural network.
///
/// It uses this parameters to approximate a value function
/// for the classic inverse pendulum problem in reinforcement problem.
/// The return value is the number of steps needed to swing up or a maximum
/// value.
#[allow(non_snake_case)]
#[no_mangle]
pub extern "C" fn G(
    x: *const libc::c_double,
    len_x: libc::c_int,
    y: *mut libc::c_double,
    len_y: libc::c_int,
) {
    let x = unsafe { std::slice::from_raw_parts(x, len_x as usize) };
    let y = unsafe { std::slice::from_raw_parts_mut(y, len_y as usize) };

    let (m1, x) = x.split_at(8 * 6);
    let (b1, x) = x.split_at(8);
    let m1 = na::MatrixSliceMN::<f64, na::U8, na::U6>::from_slice(m1);
    let b1 = na::VectorSliceN::<f64, na::U8>::from_slice(b1);

    let (m2, x) = x.split_at(16 * 8);
    let (b2, x) = x.split_at(16);
    let m2 = na::MatrixSliceMN::<f64, na::U16, na::U8>::from_slice(m2);
    let b2 = na::VectorSliceN::<f64, na::U16>::from_slice(b2);

    let (m3, _x) = x.split_at(3 * 16);
    let m3 = na::MatrixSliceMN::<f64, na::U3, na::U16>::from_slice(m3);

    let mut env = acrobot::Acrobot::new();
    let mut reward = 0.0;
    let mut obs = env.reset();
    for _ in 0..500 {
        let o1 = (m1 * obs + b1).map(|i| i.max(0.0));
        let o2 = (m2 * o1 + b2).map(|i| i.max(0.0));
        let o3 = m3 * o2;

        let action = match o3
            .iter()
            .enumerate()
            .max_by(|(_i, a), (_j, b)| a.partial_cmp(b).unwrap_or(Equal))
            .unwrap()
            .0
        {
            0 => acrobot::AcrobotAction::Left,
            1 => acrobot::AcrobotAction::No,
            2 => acrobot::AcrobotAction::Right,
            _ => unreachable!(),
        };
        let (new_obs, r, terminal) = env.step(action);
        obs = new_obs;
        reward += r;
        if terminal {
            break;
        }
    }
    // println!("reward = {}", reward);
    y[0] = reward;
}
