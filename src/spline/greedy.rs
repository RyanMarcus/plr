use crate::util::{Point};


fn greedy_spline_corridor(data: &[(f64, f64)], err: f64) -> usize {
    assert!(data.len() >= 2);

    if data.len() == 2 {
        return 1;
    }

    let pt1 = Point::from_tuple(data[0]);
    let pt2 = Point::from_tuple(data[1]);

    let mut upper = pt1.line_to(&pt2.upper_bound(err)).slope();
    let mut lower = pt1.line_to(&pt2.lower_bound(err)).slope();

    for (idx, pt_tuple) in data[2..].iter().enumerate() {
        let pt = Point::from_tuple(*pt_tuple);

        // check to make sure the line from pt1 to our new point
        // respects the current bounds.
        let l = pt1.line_to(&pt);
        if l.slope() > upper || l.slope() < lower {
            // we cannot fit the bound.
            return idx + 1;
        }

        let potential_upper = pt1.line_to(&pt.upper_bound(err)).slope();
        let potential_lower = pt1.line_to(&pt.lower_bound(err)).slope();
        
        if potential_upper < upper {
            upper = potential_upper;
        }

        if potential_lower > lower {
            lower = potential_lower;
        }
    }
    
    return data.len() - 1;
}

pub fn greedy(data: &[(f64, f64)], err: f64) -> Vec<(f64, f64)> {
    assert!(data.len() >= 2);
    if data.len() == 2 {
        return vec![data[0].clone(), data[1].clone()];
    }
    
    let mut pts = Vec::new();
    let mut last_idx = 0;

    pts.push(data[0].clone());
    while data[last_idx..].len() > 1 {
        let next_point = greedy_spline_corridor(&data[last_idx..], err) + last_idx;
        debug_assert!(next_point > last_idx);
        pts.push(data[next_point].clone());
        last_idx = next_point;
    }

    return pts;
}

#[cfg(test)]
mod test {
    use crate::test_util::*;
    use super::*;
    
    #[test]
    fn test_sin() {
        let data = sin_data();
        let pts = greedy(&data, 0.0005);
        verify_gamma_splines(0.0005, &data, &pts);
    }

    #[test]
    fn test_linear() {
        let data = linear_data(10.0, 25.0);
        let pts = greedy(&data, 0.0005);
        assert_eq!(pts.len(), 2);
        verify_gamma_splines(0.0005, &data, &pts);
    }

    #[test]
    fn test_precision() {
        let data = precision_data();
        let pts = greedy(&data, 0.0005);
        assert_eq!(pts.len(), 2);
        verify_gamma_splines(0.0005, &data, &pts);
    }
}
