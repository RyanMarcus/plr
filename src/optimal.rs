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
use crate::util::{Point, Line, Segment};
use std::collections::VecDeque;

struct Hull {
    upper: bool,
    data: VecDeque<Point>
}

impl Hull {
    fn new_upper() -> Hull {
        return Hull { upper: true, data: VecDeque::new() };
    }

    fn new_lower() -> Hull {
        return Hull { upper: false, data: VecDeque::new() };
    }

    fn remove_front(&mut self, count: usize) {
        for _ in 0..count { self.data.pop_front().unwrap(); }
    }

    fn push(&mut self, pt: Point) {
        self.data.push_back(pt);

        let mut keep_points = VecDeque::new();

        while self.data.len() > 2 {
            let pt1 = self.data.pop_back().unwrap();
            let pt2 = self.data.pop_back().unwrap();
            let pt3 = self.data.pop_back().unwrap();

            let line = pt1.line_to(&pt3).unwrap();

            if (self.upper && pt2.above(&line)) || (!self.upper && pt2.below(&line)) {
                // remove p2
                self.data.push_back(pt3);
                self.data.push_back(pt1);
            } else {
                keep_points.push_front(pt1);
                self.data.push_back(pt3);
                self.data.push_back(pt2);
            }
        }

        while !self.data.is_empty() {
            keep_points.push_front(self.data.pop_back().unwrap());
        }

        self.data = keep_points;
    }

    fn items(&self) -> &VecDeque<Point> {
        return &self.data;
    }
}

#[derive(PartialEq)]
enum OptimalState {
    Need2, Need1, Ready
}

/// Performs an optimal piecewise linear regression (PLR) in an online fashion. This approach uses constant
/// time for each call to [`process`](#method.process) and potentially linear space. In practice,
/// the space used is generally a small fraction of the data. If constant space is required for your
/// applications, please see [`GreedyPLR`](struct.GreedyPLR.html).
///
/// Each call to [`process`](#method.process) consumes a single point. Each time it is called,
/// [`process`](#method.process) returns either a [`Segment`](struct.Segment.html) representing
/// a piece of the final regression, or `None`. If your stream of points terminates, you can call
/// [`finish`](#method.finish) to flush the buffer and return the final segment.
///
/// # Example
/// ```
/// use plr::OptimalPLR;
/// // first, generate some data points...
/// let mut data = Vec::new();
///  
/// for i in 0..1000 {
///     let x = (i as f64) / 1000.0 * 7.0;
///     let y = f64::sin(x);
///     data.push((x, y));
/// }
///  
/// let mut plr = OptimalPLR::new(0.0005); // gamma = 0.0005, the maximum regression error
///  
/// let mut segments = Vec::new();
///  
/// for (x, y) in data {
///     // when `process` returns a segment, we should add it to our list
///     if let Some(segment) = plr.process(x, y) {
///         segments.push(segment);
///     }
/// }
///
/// // because we have a finite amount of data, we flush the buffer and get the potential
/// // last segment.
/// if let Some(segment) = plr.finish() {
///     segments.push(segment);
/// }
///
/// // the `segments` vector now contains all segments for this regression.
/// ```
pub struct OptimalPLR {
    state: OptimalState,
    gamma: f64,
    s0: Option<Point>,
    s1: Option<Point>,
    s_last: Option<Point>,
    rho_lower: Option<Line>,
    rho_upper: Option<Line>,
    upper_hull: Option<Hull>,
    lower_hull: Option<Hull>
        
}

impl OptimalPLR {
    /// Enables performing PLR using an optimal algorithm with a fixed gamma (maximum error).
    ///
    /// # Examples
    ///
    /// To perform an optimal PLR with a maximum error of `0.05`:
    /// ```
    /// use plr::OptimalPLR;
    /// let plr = OptimalPLR::new(0.05);
    /// ```
    pub fn new(gamma: f64) -> OptimalPLR {
        return OptimalPLR {
            state: OptimalState::Need2,
            gamma,
            s0: None, s1: None, s_last: None,
            rho_lower: None, rho_upper: None,
            upper_hull: None, lower_hull: None
        };
    }

    fn setup(&mut self) {
        // we have two points, initialize rho lower and rho upper.
        let gamma = self.gamma;
        let s0 = self.s0.as_ref().unwrap();
        let s1 = self.s1.as_ref().unwrap();

        self.rho_lower = s0.upper_bound(gamma).line_to(&s1.lower_bound(gamma));
        self.rho_upper = s0.lower_bound(gamma).line_to(&s1.upper_bound(gamma));

        let mut upper_hull = Hull::new_upper();
        let mut lower_hull = Hull::new_lower();

        upper_hull.push(s0.upper_bound(gamma));
        upper_hull.push(s1.upper_bound(gamma));
        lower_hull.push(s0.lower_bound(gamma));
        lower_hull.push(s1.lower_bound(gamma));
        
        self.upper_hull = Some(upper_hull);
        self.lower_hull = Some(lower_hull);
    }

    fn current_segment(&self, end: f64) -> Segment {
        assert!(self.state == OptimalState::Ready);
        let sint = Line::intersection(self.rho_lower.as_ref().unwrap(),
                                      self.rho_upper.as_ref().unwrap());
        let segment_start = self.s0.as_ref().unwrap().to_tuple().0;
        let segment_stop = end;

        let avg_slope = Line::average_slope(self.rho_lower.as_ref().unwrap(),
                                            self.rho_upper.as_ref().unwrap());
        
        let (sint_x, sint_y) = sint.unwrap().to_tuple();
        let intercept = -avg_slope * sint_x + sint_y;
        return Segment {
            start: segment_start, stop: segment_stop,
            slope: avg_slope, intercept
        };
    }

    fn process_pt(&mut self, pt: Point) -> Option<Segment> {
        assert!(self.state == OptimalState::Ready);
        if !(pt.above(self.rho_lower.as_ref().unwrap()) && pt.below(self.rho_upper.as_ref().unwrap())) {
            // we cannot adjust either extreme slope to fit this point, we have to
            // start a new segment.
            let current_segment = self.current_segment(pt.to_tuple().0);

            self.s0 = Some(pt);
            self.state = OptimalState::Need1;
            return Some(current_segment);
        }

        // otherwise, we can adjust the extreme slopes to fit the point.
        let s_upper = pt.upper_bound(self.gamma);
        let s_lower = pt.lower_bound(self.gamma);
        if s_upper.below(self.rho_upper.as_ref().unwrap()) {
            let lower_hull: &mut Hull = self.lower_hull.as_mut().unwrap();
            
            // find the point in the lower hull that would minimize the slope
            // between that point and s_upper
            let it = lower_hull.items().iter().enumerate()
                .map(|(idx, pt)| (idx, pt.line_to(&s_upper).unwrap()));

            // get the min, because we can't use .min() on f64
            let mut curr_best_val = std::f64::INFINITY;
            let mut curr_best_idx = 0;
            for (idx, line) in it {
                if line.slope() < curr_best_val {
                    curr_best_idx = idx;
                    curr_best_val = line.slope();
                }
            }

            self.rho_upper = s_upper.line_to(
                &lower_hull.items()[curr_best_idx]
            );
            
            lower_hull.remove_front(curr_best_idx);
            lower_hull.push(s_lower.clone());
        }

        if s_lower.above(self.rho_lower.as_ref().unwrap()) {
            let upper_hull: &mut Hull = self.upper_hull.as_mut().unwrap();
            
            // find the point in the upper hull that would maximize the slope
            // between that point and s_upper
            let it = upper_hull.items().iter().enumerate()
                .map(|(idx, pt)| (idx, pt.line_to(&s_lower).unwrap()));

            // get the max, because we can't use .max() on f64
            let mut curr_best_val = -std::f64::INFINITY;
            let mut curr_best_idx = 0;
            for (idx, line) in it {
                if line.slope() > curr_best_val {
                    curr_best_idx = idx;
                    curr_best_val = line.slope();
                }
            }

            self.rho_lower = s_lower.line_to(
                &upper_hull.items()[curr_best_idx]
            );

            upper_hull.remove_front(curr_best_idx);
            upper_hull.push(s_upper.clone());
        }

        return None
    }

    /// Processes a single point using the optimal PLR algorithm. This function returns
    /// a new [`Segment`](struct.Segment.html) when the current segment cannot accommodate
    /// the passed point, and returns None if the current segment could be (greedily) adjusted to
    /// fit the point.
    pub fn process(&mut self, x: f64, y: f64) -> Option<Segment> {
        let pt = Point::new(x, y);
        self.s_last = Some(pt.clone());
        
        let mut returned_segment: Option<Segment> = None;

        let new_state = match self.state {
            OptimalState::Need2 => { self.s0 = Some(pt); OptimalState::Need1 },
            OptimalState::Need1 => { self.s1 = Some(pt); self.setup(); OptimalState::Ready }
            OptimalState::Ready => {
                returned_segment = self.process_pt(pt);

                if returned_segment.is_some() {
                    OptimalState::Need1
                } else {
                    OptimalState::Ready
                }
            }
        };

        self.state = new_state;
        return returned_segment;
    }

    /// Terminates the PLR process, returning a final segment if required.
    pub fn finish(self) -> Option<Segment> {
        return match self.state {
            OptimalState::Need2 => None,
            OptimalState::Need1 => {
                let s0 = self.s0.unwrap().to_tuple();
                Some(Segment { start: s0.0, stop: s0.0 + 1.0, slope: 0.0, intercept: s0.1 })
            },
            OptimalState::Ready => Some(self.current_segment(
                self.s_last.as_ref().unwrap().to_tuple().0 + 1.0))
        };
    }
}


#[cfg(test)]
mod test {
    use crate::test_util::*;
    use approx::*;
    use super::*;

    #[test]
    fn test_upper_hull() {
        let mut hull = Hull::new_upper();

        hull.push(Point::new(1.0, 1.0));
        hull.push(Point::new(2.0, 1.0));
        hull.push(Point::new(3.0, 3.0));
        hull.push(Point::new(4.0, 3.0));

        assert_eq!(hull.items().len(), 3);

        let items = hull.items();

        // (1.0, 1.0)
        assert_relative_eq!(items[0].to_tuple().0, 1.0);
        assert_relative_eq!(items[0].to_tuple().1, 1.0);

        // (2.0, 1.0)
        assert_relative_eq!(items[1].to_tuple().0, 2.0);
        assert_relative_eq!(items[1].to_tuple().1, 1.0);

        // (4.0, 3.0)
        assert_relative_eq!(items[2].to_tuple().0, 4.0);
        assert_relative_eq!(items[2].to_tuple().1, 3.0);
    }

    #[test]
    fn test_lower_hull() {
        let mut hull = Hull::new_lower();

        hull.push(Point::new(1.0, 1.0));
        hull.push(Point::new(2.0, 1.0));
        hull.push(Point::new(3.0, 3.0));
        hull.push(Point::new(4.0, 3.0));

        assert_eq!(hull.items().len(), 3);

        let items = hull.items();

        // (1.0, 1.0)
        assert_relative_eq!(items[0].to_tuple().0, 1.0);
        assert_relative_eq!(items[0].to_tuple().1, 1.0);

        // (3.0, 3.0)
        assert_relative_eq!(items[1].to_tuple().0, 3.0);
        assert_relative_eq!(items[1].to_tuple().1, 3.0);

        // (4.0, 3.0)
        assert_relative_eq!(items[2].to_tuple().0, 4.0);
        assert_relative_eq!(items[2].to_tuple().1, 3.0);
    }
        
    #[test]
    fn test_sin() {
        let mut plr = OptimalPLR::new(0.0005);
        let data = sin_data();

        let mut segments = Vec::new();
        
        for &(x, y) in data.iter() {
            if let Some(segment) = plr.process(x, y) {
                segments.push(segment);
            }
        }

        if let Some(segment) = plr.finish() {
            segments.push(segment);
        }
        
        assert_eq!(segments.len(), 66);
        verify_gamma(0.0005, &data, &segments);
    }

    #[test]
    fn test_linear() {
        let mut plr = OptimalPLR::new(0.0005);
        let data = linear_data(10.0, 25.0);

        let mut segments = Vec::new();
        
        for &(x, y) in data.iter() {
            if let Some(segment) = plr.process(x, y) {
                segments.push(segment);
            }
        }

        if let Some(segment) = plr.finish() {
            segments.push(segment);
        }

        assert_eq!(segments.len(), 1);
        verify_gamma(0.0005, &data, &segments);
    }

}
