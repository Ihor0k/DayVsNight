use std::f32::consts::PI;
use std::ops::{Add, AddAssign, Mul, Sub};

use crate::FIELD_SIZE;

const BALL_RADIUS: f32 = 0.5;

pub struct State {
    pub field: Field,
    pub day_ball: Ball,
    pub night_ball: Ball,
}

#[derive(Debug)]
pub struct Field {
    pub grid: [[bool; FIELD_SIZE]; FIELD_SIZE],
    pub day_count: u32,
    pub night_count: u32,
}

#[derive(Debug)]
pub struct Ball {
    pub pos: PointF,
    pub vec: PointF,
    pub side: bool,
}

#[derive(Debug, Clone, Copy)]
struct PointI {
    x: i32,
    y: i32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointF {
    pub x: f32,
    pub y: f32,
}

impl Default for Field {
    fn default() -> Self {
        let half_size = FIELD_SIZE / 2;
        let grid = (0..FIELD_SIZE).map(|_| {
            (0..FIELD_SIZE).map(|x| x < half_size).collect::<Vec<bool>>().try_into().unwrap()
        }).collect::<Vec<_>>().try_into().unwrap();

        let half_cell_num = (FIELD_SIZE * FIELD_SIZE / 2) as u32;
        Field {
            grid,
            day_count: half_cell_num,
            night_count: half_cell_num,
        }
    }
}

impl Field {
    fn invert(&mut self, cell: PointI) {
        if self.is_in_boundaries(cell) {
            let (x, y) = (cell.x as usize, cell.y as usize);
            self.grid[y][x] = !self.grid[y][x];
            if self.grid[y][x] {
                self.day_count += 1;
                self.night_count -= 1;
            } else {
                self.day_count -= 1;
                self.night_count += 1;
            }
        }
    }

    fn has_value(&self, cell: PointI, value: bool) -> bool {
        self.is_in_boundaries(cell) && self.grid[cell.y as usize][cell.x as usize] == value
    }

    fn is_in_boundaries(&self, cell: PointI) -> bool {
        (0..FIELD_SIZE as i32).contains(&cell.x) && (0..FIELD_SIZE as i32).contains(&cell.y)
    }
}

impl Ball {
    fn new(
        pos: PointF,
        vec: PointF,
        side: bool,
    ) -> Self {
        Self {
            pos,
            vec: vec.normalized(),
            side,
        }
    }
}

impl PointI {
    fn new(x: i32, y: i32) -> Self {
        Self { x, y }
    }
}

impl Add for PointF {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self { x: self.x + rhs.x, y: self.y + rhs.y }
    }
}

impl Sub for PointF {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self { x: self.x - rhs.x, y: self.y - rhs.y }
    }
}

impl AddAssign<PointF> for PointF {
    fn add_assign(&mut self, rhs: PointF) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

impl Mul<f32> for PointF {
    type Output = PointF;

    fn mul(self, rhs: f32) -> Self::Output {
        Self { x: self.x * rhs, y: self.y * rhs }
    }
}

impl PointF {
    fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    fn normalized(&self) -> Self {
        let len = self.length();
        Self::new(self.x / len, self.y / len)
    }

    fn length(&self) -> f32 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }

    fn rotated(&self, angle: f32) -> Self {
        let sin = angle.sin();
        let cos = angle.cos();
        Self::new(
            self.x * cos - self.y * sin,
            self.x * sin + self.y * cos,
        )
    }

    fn angle_to(&self, other: &PointF) -> f32 {
        self.rotated(-PI / 2.).dot(other).atan2(self.dot(other))
    }

    fn dot(&self, other: &PointF) -> f32 {
        self.x * other.x + self.y * other.y
    }

    fn signum_i32(&self) -> PointI {
        PointI::new(self.x.signum() as i32, self.y.signum() as i32)
    }

    fn dist_to(&self, other: &PointF) -> f32 {
        ((other.x - self.x).powi(2) + (other.y - self.y).powi(2)).sqrt()
    }
}

impl Default for State {
    fn default() -> Self {
        let field = Field::default();

        let day_x = macroquad::rand::gen_range(BALL_RADIUS, (FIELD_SIZE / 2) as f32 - BALL_RADIUS);
        let day_y = macroquad::rand::gen_range(BALL_RADIUS, FIELD_SIZE as f32 - BALL_RADIUS);
        let day_pos = PointF::new(day_x, day_y);
        let day_vec = PointF::new(macroquad::rand::gen_range(0., 1.), macroquad::rand::gen_range(0., 1.));
        let day_ball = Ball::new(day_pos, day_vec, true);

        let night_x = macroquad::rand::gen_range((FIELD_SIZE / 2) as f32 + BALL_RADIUS, FIELD_SIZE as f32 - BALL_RADIUS);
        let night_y = macroquad::rand::gen_range(BALL_RADIUS, FIELD_SIZE as f32 - BALL_RADIUS);
        let night_pos = PointF::new(night_x, night_y);
        let night_vec = PointF::new(macroquad::rand::gen_range(0., 1.), macroquad::rand::gen_range(0., 1.));
        let night_ball = Ball::new(night_pos, night_vec, false);

        State { field, day_ball, night_ball }
    }
}

impl State {
    pub fn next_step(&mut self, dist: f32) {
        ball_next_step(self, dist, true);
        ball_next_step(self, dist, false);
    }
}

fn ball_next_step(state: &mut State, dist: f32, side: bool) {
    if dist <= 0. { return; }
    let (ball, other_ball) = if side { (&mut state.day_ball, &state.night_ball) } else { (&mut state.night_ball, &state.day_ball) };
    let field = &mut state.field;
    let ball_pos = ball.pos;
    let ball_vec = ball.vec;
    let ball_vec_sign = ball_vec.signum_i32();
    let ball_displacement = ball_vec * dist;
    let ball_side = ball.side;

    let ball_left_pos_vec = ball_vec.rotated(-PI / 2.);
    let ball_right_pos_vec = ball_vec.rotated(PI / 2.);

    let (ball_side_x_vec, ball_side_y_vec) = match ball_vec_sign {
        PointI { x: -1, y: -1 } => (PointF::new(-1., 0.), PointF::new(0., -1.)),
        PointI { x: -1, y: 1 } => (PointF::new(-1., 0.), PointF::new(0., 1.)),
        PointI { x: 1, y: -1 } => (PointF::new(1., 0.), PointF::new(0., -1.)),
        PointI { x: 1, y: 1 } => (PointF::new(1., 0.), PointF::new(0., 1.)),
        _ => panic!()
    };

    let ball_left_pos = ball_pos + ball_left_pos_vec * BALL_RADIUS;
    let ball_side_x_pos = ball_pos + ball_side_x_vec * BALL_RADIUS;
    let ball_side_y_pos = ball_pos + ball_side_y_vec * BALL_RADIUS;
    let ball_right_pos = ball_pos + ball_right_pos_vec * BALL_RADIUS;

    let mut contact_points: Option<(PointF, PointF, PointI)> = None;

    let ball_left_walk_result = walk_grid(ball_left_pos, ball_displacement, ball_side, field);
    let ball_side_x_walk_result = walk_grid(ball_side_x_pos, ball_displacement, ball_side, field);
    let ball_side_y_walk_result = walk_grid(ball_side_y_pos, ball_displacement, ball_side, field);
    let ball_right_walk_result = walk_grid(ball_right_pos, ball_displacement, ball_side, field);
    let min_dist_result_option: Option<(PointF, PointI, f32)> = [
        ball_left_walk_result,
        ball_side_x_walk_result,
        ball_side_y_walk_result,
        ball_right_walk_result,
    ]
        .iter()
        .filter_map(|x| *x)
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
    if let Some((min_dist_result_point, min_dist_result_cell, _)) = min_dist_result_option {
        let ball_side_x_result_point = ball_side_x_walk_result.map(|x| x.0);
        let ball_side_y_result_point = ball_side_y_walk_result.map(|x| x.0);

        let cell_contact_point = find_cell_contact_point(min_dist_result_point, ball_side_x_result_point, ball_side_y_result_point, ball_vec_sign);
        let ball_contact_point = find_ball_contact_point(ball_left_pos, ball_right_pos, ball_vec, cell_contact_point);

        contact_points = Some((cell_contact_point, ball_contact_point, min_dist_result_cell));
    }
    if contact_points.is_none() {
        let target_pos = ball_pos + ball_displacement;
        let cell_pos_x = target_pos.x as i32 + if target_pos.x.fract() >= 0.5 { 1 } else { -1 };
        let cell_pos_y = target_pos.y as i32 + if target_pos.y.fract() >= 0.5 { 1 } else { -1 };
        let cell_pos = PointI::new(cell_pos_x, cell_pos_y);
        if field.has_value(cell_pos, !ball_side) {
            if let Some(cell_contact_point) = detect_collisions(target_pos, cell_pos) {
                let ball_contact_point = find_ball_contact_point(ball_left_pos, ball_right_pos, ball_vec, cell_contact_point);
                contact_points = Some((cell_contact_point, ball_contact_point, cell_pos));
            }
        }
    }
    if let Some((cell_contact_point, ball_contact_point, cell_pos)) = contact_points {
        let ball_pos_delta = cell_contact_point - ball_contact_point;
        let new_ball_pos = ball_pos + ball_pos_delta;

        let ball_center_to_contact_point_vec_normalized = (ball_contact_point - ball_pos).normalized();
        let angle_between_ball_vec_and_ball_contact_point_vec = ball_vec.angle_to(&ball_center_to_contact_point_vec_normalized);
        let ball_vec_rotation_angle = PI - 2. * angle_between_ball_vec_and_ball_contact_point_vec;
        let new_ball_vec = ball_vec.rotated(ball_vec_rotation_angle);

        if detect_collisions(other_ball.pos, cell_pos).is_none() {
            field.invert(cell_pos);
        };

        let dist_moved = ball_pos_delta.length();
        if dist_moved > dist {
            println!("{:?}", ball);
            println!("{:?}", other_ball);
            println!("{:?}", field);
            println!("{:?}", min_dist_result_option);
            println!("cell_cp={:?} ball_cp={:?}", cell_contact_point, ball_contact_point);
            println!("dist={}, dist_moved={}", dist, dist_moved);
        }
        ball.pos = new_ball_pos;
        ball.vec = new_ball_vec;
        ball_next_step(state, dist - dist_moved, side);
    } else {
        ball.pos += ball_displacement;
    }
}

fn detect_collisions(ball_pos: PointF, cell: PointI) -> Option<PointF> {
    let test_x = if ball_pos.x < cell.x as f32 {
        cell.x as f32
    } else if ball_pos.x > cell.x as f32 + 1. {
        cell.x as f32 + 1.
    } else {
        ball_pos.x
    };
    let test_y = if ball_pos.y < cell.y as f32 {
        cell.y as f32
    } else if ball_pos.y > cell.y as f32 + 1. {
        cell.y as f32 + 1.
    } else {
        ball_pos.y
    };
    let test_pos = PointF::new(test_x, test_y);
    let dist = ball_pos.dist_to(&test_pos);
    if dist <= BALL_RADIUS { Some(test_pos) } else { None }
}

fn find_cell_contact_point(
    min_dist_point: PointF,
    ball_side_x_walk_point: Option<PointF>,
    ball_side_y_walk_point: Option<PointF>,
    ball_vec_sign: PointI,
) -> PointF {
    let (round_x_dir, round_y_dir) = if Some(min_dist_point) == ball_side_x_walk_point {
        if min_dist_point.x.fract() == 0. {
            (0, 0)
        } else {
            (-ball_vec_sign.x, 0)
        }
    } else if Some(min_dist_point) == ball_side_y_walk_point {
        if min_dist_point.y.fract() == 0. {
            (0, 0)
        } else {
            (0, -ball_vec_sign.y)
        }
    } else {
        if min_dist_point.x.fract() == 0. {
            (0, ball_vec_sign.y)
        } else {
            (ball_vec_sign.x, 0)
        }
    };
    PointF::new(
        round_coordinate(min_dist_point.x, round_x_dir),
        round_coordinate(min_dist_point.y, round_y_dir),
    )
}

fn round_coordinate(x: f32, dir: i32) -> f32 {
    match dir {
        -1 => x.floor(),
        0 => x,
        1 => x.ceil(),
        _ => panic!()
    }
}

fn find_ball_contact_point(ball_left: PointF, ball_right: PointF, ball_vec: PointF, point: PointF) -> PointF {
    let diameter_vec = ball_right - ball_left;
    let projection_onto_diameter = ball_left + diameter_vec * diameter_vec.dot(&(point - ball_left));
    let ball_left_to_projection_onto_diameter = (projection_onto_diameter - ball_left).length();
    let ball_right_to_projection_onto_diameter = (projection_onto_diameter - ball_right).length();
    let projection_to_ball_edge = (ball_left_to_projection_onto_diameter * ball_right_to_projection_onto_diameter).sqrt();
    projection_onto_diameter + ball_vec * projection_to_ball_edge
}

fn walk_grid(from: PointF, dir: PointF, ball_side: bool, field: &Field) -> Option<(PointF, PointI, f32)> {
    let to = from + dir;
    let dir_sign = dir.signum_i32();
    let mut point = from;
    loop {
        let dx_to_next_x = dist_to_next_cell(point.x, dir_sign.x);
        let dy_to_next_x = dx_to_next_x / dir.x * dir.y;
        let dir_to_next_x = PointF::new(dx_to_next_x, dy_to_next_x);

        let dy_to_next_y = dist_to_next_cell(point.y, dir_sign.y);
        let dx_to_next_y = dy_to_next_y / dir.y * dir.x;
        let dir_to_next_y = PointF::new(dx_to_next_y, dy_to_next_y);

        point += std::cmp::min_by(dir_to_next_x, dir_to_next_y, |a, b| a.length().partial_cmp(&b.length()).unwrap());
        let is_x_reached = if dir_sign.x == 1 { point.x >= to.x } else { point.x <= to.x };
        let is_y_reached = if dir_sign.y == 1 { point.y >= to.y } else { point.y <= to.y };
        let cell = PointI::new(cell_coordinate(point.x, dir_sign.x), cell_coordinate(point.y, dir_sign.y));
        let field_size = FIELD_SIZE as f32;
        let is_border_x_reached = point.x >= field_size || point.x <= 0.;
        let is_border_y_reached = point.y >= field_size || point.y <= 0.;
        if is_x_reached && is_y_reached {
            return None;
        } else if is_border_x_reached || is_border_y_reached || field.has_value(cell, !ball_side) {
            let dist = from.dist_to(&point);
            return Some((point, cell, dist));
        }
    }
}

fn dist_to_next_cell(x: f32, sign: i32) -> f32 {
    let fract = x.fract();
    if sign == 1 {
        1. - fract
    } else if fract == 0. {
        -1.
    } else {
        -fract
    }
}

fn cell_coordinate(x: f32, sign: i32) -> i32 {
    let x_i32 = x as i32;
    if sign == -1 && x.fract() == 0. {
        x_i32 - 1
    } else {
        x_i32
    }
}