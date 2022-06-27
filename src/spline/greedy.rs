use crate::util::Point;

/// A greedy spline regressor using the greedy spline corridor algorithm.
/// Each call to `process` will return a new spline point if one is needed
/// to maintain the error bound, or `None` if the passed point works with
/// the current spline.
///
/// Note that the first point (passed to `new`) won't be returned by
/// `process`, but still needs to be included in the spline knots.
///
/// Call `finish` to get the final point of the spline.
#[derive(Clone)]
pub struct GreedySpline {
    error: f64,
    pt1: Point,
    pt2: Option<Point>,
    slope_ub: f64,
    slope_lb: f64,
}

impl GreedySpline {
    /// Constructs a new online greedy spline regression. The first point,
    /// passed to this function, must be used as the first spline point.
    pub fn new(pt1: Point, err: f64) -> GreedySpline {
        GreedySpline {
            error: err,
            pt1,
            pt2: None,
            slope_lb: std::f64::NEG_INFINITY,
            slope_ub: std::f64::INFINITY,
        }
    }

    /// Add another point the spline regression, potentially returning a new knot.
    pub fn process(&mut self, pt: &Point) -> Option<Point> {
        match self.pt2 {
            None => {
                self.pt2 = Some(*pt);
                return None;
            }

            Some(pt2) => {
                let potential_upper = self.pt1.line_to(&pt.upper_bound(self.error)).slope();
                let potential_lower = self.pt1.line_to(&pt.lower_bound(self.error)).slope();

                if potential_upper > self.slope_ub || potential_lower < self.slope_lb {
                    // this point is outside our bound, have to create a new spline point
                    self.pt1 = pt2;
                    self.slope_lb = std::f64::NEG_INFINITY;
                    self.slope_ub = std::f64::INFINITY;
                    self.pt2 = None;
                    return Some(pt2);
                }

                // otherwise, we have a new pt2 and need to update the slope bounds
                self.pt2 = Some(*pt);
                self.slope_lb = f64::max(self.slope_lb, potential_lower);
                self.slope_ub = f64::min(self.slope_ub, potential_upper);

                return None;
            }
        };
    }

    /// Finish the spline regression, returning the final knot.
    pub fn finish(self) -> Point {
        match self.pt2 {
            Some(pt) => pt,
            None => self.pt1,
        }
    }
}

/// Learns a spline regression over `data` with an error bound of `err`.
/// Each spline point is returned in (x, y) coordinates.
pub fn greedy_spline(data: &[(f64, f64)], err: f64) -> Vec<(f64, f64)> {
    assert!(data.len() >= 2);
    if data.len() == 2 {
        return vec![data[0], data[1]];
    }

    let mut pts = vec![data[0]];
    let mut gs = GreedySpline::new(Point::from_tuple(data[0]), err);

    for pt in &data[1..] {
        if let Some(knot) = gs.process(&Point::from_tuple(*pt)) {
            pts.push(knot.as_tuple());
        }
    }
    pts.push(gs.finish().as_tuple());
    return pts;
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_util::*;

    #[test]
    fn test_sin() {
        let data = sin_data();
        let pts = greedy_spline(&data, 0.0005);
        assert!(pts.len() < 500);
        verify_gamma_splines(0.0005, &data, &pts);
    }

    #[test]
    fn test_linear() {
        let data = linear_data(10.0, 25.0);
        let pts = greedy_spline(&data, 0.0005);
        assert_eq!(pts.len(), 2);
        verify_gamma_splines(0.0005, &data, &pts);
    }

    #[test]
    fn test_precision() {
        let data = precision_data();
        let pts = greedy_spline(&data, 0.0005);
        assert_eq!(pts.len(), 2);
        verify_gamma_splines(0.0005, &data, &pts);
    }

    #[test]
    fn test_osm() {
        let data = osm_data();
        let pts = greedy_spline(&data, 64.0);
        verify_gamma_splines(64.0, &data, &pts);
    }

    #[test]
    fn test_fb() {
        let data = fb_data();
        let pts = greedy_spline(&data, 64.0);
        verify_gamma_splines(64.0, &data, &pts);
    }
}
