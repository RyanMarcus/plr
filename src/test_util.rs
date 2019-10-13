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

        println!("{} {} with {:?}", x, y, seg);
        assert!(f64::abs(pred - y) <= gamma,
                "Prediction of {} was not within gamma ({}) of true value {}",
                pred, gamma, y);
    }

}
