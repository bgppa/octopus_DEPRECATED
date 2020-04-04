//! This module implements the acrobot environment.
//! The code was ported from [gym](https://github.com/openai/gym/blob/master/gym/envs/classic_control/acrobot.py).
use std::f64::consts::PI;
const LINK_LENGTH_1: f64 = 1.;
const LINK_MASS_1: f64 = 1.;
const LINK_MASS_2: f64 = 1.;
/// Position of the center of mass of link 1
const LINK_COM_POS_1: f64 = 0.5;
/// Position of the center of mass of link 2
const LINK_COM_POS_2: f64 = 0.5;
/// Moments of intertia for both links
const LINK_MOI: f64 = 1.;
const MAX_VEL_1: f64 = 4. * PI;
const MAX_VEL_2: f64 = 9. * PI;
const AVAIL_TORQUE: [f64; 3] = [-1., 0., 1.];

use nalgebra as na;

type State = na::Vector5<f64>;

pub(crate) type Observation = na::Vector6<f64>;

pub(crate) enum AcrobotAction {
    Left,
    No,
    Right,
}

pub(crate) struct Acrobot {
    state: [f64; 4],
    dt: f64,
}

impl Acrobot {
    pub(crate) fn new() -> Acrobot {
        Acrobot {
            // deterministic start for now
            state: [0.; 4],
            dt: 0.2,
        }
    }

    pub(crate) fn reset(&mut self) -> Observation {
        self.state = [0.; 4];
        self.observation()
    }

    pub(crate) fn step(&mut self, action: AcrobotAction) -> (Observation, f64, bool) {
        let s_augmented: State = {
            let torque = match action {
                AcrobotAction::Left => AVAIL_TORQUE[0],
                AcrobotAction::No => AVAIL_TORQUE[1],
                AcrobotAction::Right => AVAIL_TORQUE[2],
            };
            let mut tmp = [0f64; 5];
            tmp[0..4].clone_from_slice(&self.state);
            tmp[4] = torque;
            tmp.into()
        };
        let ns = {
            let mut ns = rk4(dsdt, s_augmented, self.dt);
            ns[0] = wrap(ns[0], -PI, PI);
            ns[1] = wrap(ns[1], -PI, PI);
            ns[2] = ns[2].clamp(-MAX_VEL_1, MAX_VEL_1);
            ns[3] = ns[3].clamp(-MAX_VEL_2, MAX_VEL_2);
            ns
        };
        self.state.clone_from_slice(&ns.as_slice()[..4]);
        let terminal = -ns[0].cos() - (ns[1] + ns[0]).cos() > 1.;
        let reward = if terminal { 0. } else { -0.1 };
        (self.observation().into(), reward, terminal)
    }

    fn observation(&self) -> Observation {
        let s = self.state;
        let obs = [s[0].cos(), s[0].sin(), s[1].cos(), s[1].sin(), s[2], s[3]];
        obs.into()
    }
}

fn dsdt(s_augmented: State) -> State {
    let m1 = LINK_MASS_1;
    let m2 = LINK_MASS_2;
    let l1 = LINK_LENGTH_1;
    let lc1 = LINK_COM_POS_1;
    let lc2 = LINK_COM_POS_2;
    let i1 = LINK_MOI;
    let i2 = LINK_MOI;
    let g = 9.8;
    let theta1 = s_augmented[0];
    let theta2 = s_augmented[1];
    let dtheta1 = s_augmented[2];
    let dtheta2 = s_augmented[3];
    let a = s_augmented[4];

    let d1 = m1 * lc1.powf(2.)
        + m2 * (l1.powf(2.) + lc2.powf(2.) + 2. * l1 * lc2 * theta2.cos())
        + i1
        + i2;
    let d2 = m2 * (lc2.powf(2.) + l1 * lc2 * theta2.cos()) + i2;
    let phi2 = m2 * lc2 * g * (theta1 + theta2 - PI / 2.).cos();
    let phi1 = -m2 * l1 * lc2 * dtheta2.powf(2.) * theta2.sin()
        - 2. * m2 * l1 * lc2 * dtheta2 * dtheta1 * theta2.sin()
        + (m1 * lc1 + m2 * l1) * g * (theta1 - PI / 2.).cos()
        + phi2;
    let ddtheta2 = (a + d2 / d1 * phi1 - m2 * l1 * lc2 * dtheta1.powf(2.) * theta2.sin() - phi2)
        / (m2 * lc2.powf(2.) + i2 - d2.powf(2.) / d1);
    let ddtheta1 = -(d2 * ddtheta2 + phi1) / d1;
    [dtheta1, dtheta2, ddtheta1, ddtheta2, 0.].into()
}

fn wrap(mut x: f64, low: f64, high: f64) -> f64 {
    let diff = high - low;
    while x > high {
        x -= diff;
    }
    while x < low {
        x += diff;
    }
    x
}

fn rk4<F>(derivs: F, y0: State, dt: f64) -> State
where
    F: Fn(State) -> State,
{
    let dt2 = dt / 2.;
    let k1 = derivs(y0);
    let k2 = derivs(y0 + dt2 * k1);
    let k3 = derivs(y0 + dt2 * k2);
    let k4 = derivs(y0 + dt * k3);
    y0 + dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
}

#[test]
fn test_dsdt() {
    let aug = [-0.05195153, 0.08536712, -0.09221591, -0.08041345, -1.].into();
    let res: State = [
        -0.09221590518472411,
        -0.08041345407137711,
        1.0863252141806476,
        -2.4505276983082602,
        0.0,
    ]
    .into();
    assert!((dsdt(aug) - res).norm() < 1e-7);
}

#[test]
fn test_rk4() {
    let y0 = [-0.05195153, 0.08536712, -0.09221591, -0.08041345, -1.].into();
    let y1: State = [-0.04859367, 0.0214236, 0.12288546, -0.54704153, -1.].into();
    assert!((rk4(dsdt, y0, 0.2) - y1).norm() < 1e-7);
}

#[test]
fn step() {
    let mut ac = Acrobot::new();
    let (obs, r, t) = ac.step(AcrobotAction::Left);
    let obs_ref: Observation = [
        0.99991205,
        0.01326258,
        0.99941225,
        -0.03428051,
        0.12866185,
        -0.33450109,
    ]
    .into();
    assert!((obs - obs_ref).norm() < 1e-7);
    assert_eq!(r, -0.1);
    assert!(!t);

    let (obs, r, t) = ac.step(AcrobotAction::No);
    let obs_ref: Observation = [
        0.99937865,
        0.03524653,
        0.99565478,
        -0.09312125,
        0.08504501,
        -0.24181155,
    ]
    .into();
    assert!((obs - obs_ref).norm() < 1e-7);
    assert_eq!(r, -0.1);
    assert!(!t);

    let (obs, r, t) = ac.step(AcrobotAction::Right);
    let obs_ref: Observation = [
        0.99948691,
        0.03202992,
        0.99573886,
        -0.09221776,
        -0.11624583,
        0.24966734,
    ]
    .into();
    assert!((obs - obs_ref).norm() < 1e-7);
    assert_eq!(r, -0.1);
    assert!(!t);
}
