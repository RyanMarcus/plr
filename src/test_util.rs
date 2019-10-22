// < begin copyright > 
// Copyright Ryan Marcus 2019
// 
// This file is part of plr.
// 
// plr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// plr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with plr.  If not, see <http://www.gnu.org/licenses/>.
// 
// < end copyright > 
use std::collections::VecDeque;
use crate::util::*;
use superslice::*;
use crate::data::*;
use approx::*;

pub fn sin_data() -> Vec<(f64, f64)> {
    let mut data = Vec::new();

    for i in 0..1000 {
        let x = (i as f64) / 1000.0 * 7.0;
        let y = f64::sin(x);
        data.push((x, y));
    }

    return data;
}

pub fn linear_data(slope: f64, intercept: f64) -> Vec<(f64, f64)> {
    let mut data = Vec::new();

    for i in 0..10 {
        let x = (i as f64) / 1000.0;
        let y = x * slope + intercept;
        data.push((x, y));
    }

    return data;
}

pub fn precision_data() -> Vec<(f64, f64)> {
    let mut data = Vec::new();

    for i in 0..1000 {
        let x = ((i as f64) / 1000.0) * f64::powi(2.0, 60);
        
        let y = i as f64;

        data.push((x, y));
    }

    return data;
}

pub fn osm_data() -> Vec<(f64, f64)> {
    return OSM_DATA.iter()
        .map(|&(x, y)| (x as f64, y as f64))
        .collect();
}

pub fn fb_data() -> Vec<(f64, f64)> {
    return FB_DATA.iter()
        .map(|&(x, y)| (x as f64, y as f64))
        .collect();
}

pub fn verify_gamma(gamma: f64, data: &[(f64, f64)], segments: &[Segment]) {
    let mut seg_q = VecDeque::new();

    for segment in segments {
        seg_q.push_back(segment);
    }
    
    for &(x, y) in data {
        while seg_q.front().unwrap().stop <= x {
            seg_q.pop_front();
        }

        let seg = seg_q.front().unwrap();

        assert!(seg.start <= x);
        assert!(seg.stop >= x);

        let line = Line::new(seg.slope, seg.intercept);
        let pred = line.at(x).y;

        if x as u64 == 42285439947654605 {
            println!("{} to {}", x, pred);
            println!("with slope: {} and intercept: {}", seg.slope, seg.intercept);
            println!("segment start: {}", seg.start);
        }

        assert!(f64::abs(pred - y) <= gamma,
                "Prediction of {} was not within gamma ({}) of true value {}",
                pred, gamma, y);
    }
}

fn spline_interpolate(pt: f64, knots: &[(f64, f64)]) -> (f64, f64) {
    let upper_idx = usize::min(knots.len() - 1,
                               knots.upper_bound_by(|x| x.0.partial_cmp(&pt).unwrap()));
    let lower_idx = upper_idx - 1;

    println!("{:?} <= {} < {:?}", knots[lower_idx], pt, knots[upper_idx]);
    return Point::from_tuple(knots[lower_idx])
        .line_to(&Point::from_tuple(knots[upper_idx]))
        .at(pt).to_tuple();
}

pub fn verify_gamma_splines(gamma: f64, data: &[(f64, f64)], pts: &[(f64, f64)]) {
    println!("{:?}", pts);
    for &(x, y) in data {
        let pred = spline_interpolate(x, pts);
        assert!(f64::abs(pred.1 - y) <= gamma,
                "Prediction of {} was not within gamma ({}) of true value {}",
                pred.1, gamma, y);
    }
}



